"""
Graph-based algorithm for removing braided rivers and resolving flowline loops.

This module provides advanced algorithms using directed graph theory to handle
complex river networks with braided channels, loops, and parallel streams.
"""

import logging
import os
from typing import List, Dict, Set, Tuple, Optional, DefaultDict
from collections import defaultdict, deque
import numpy as np
from enum import Enum
import time

# Import for spatial indexing (R-tree)
try:
    from rtree.index import Index as RTreeindex
    HAS_RTREE = True
except ImportError:
    HAS_RTREE = False
    logger.warning("R-tree not available - using fallback spatial indexing")

HAS_NETWORKX = False

from pyflowline.classes.flowline import pyflowline
from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.edge import pyedge
from pyflowline.classes.confluence import pyconfluence
from pyflowline.formats.export_flowline import export_flowline_to_geojson

# Import split_flowline dependencies
import importlib.util
iFlag_cython = importlib.util.find_spec("cython")
if iFlag_cython is not None:
    from pyflowline.algorithms.cython.kernel import find_vertex_on_edge
else:
    from pyflowline.algorithms.auxiliary.find_index_in_list import find_vertex_on_edge

# Set up logger
logger = logging.getLogger(__name__)

class GraphChangeType(Enum):
    """Types of graph changes for tracking and update strategy decisions."""
    # Simple operations - can use incremental updates
    ADD_FLOWLINE = "add_flowline"
    REMOVE_FLOWLINE = "remove_flowline"
    UPDATE_ATTRIBUTES = "update_attributes"

    # Moderate operations - conditional updates
    REMOVE_SMALL_RIVERS = "remove_small_rivers"
    MERGE_LINEAR_SEGMENTS = "merge_linear_segments"
    SPLIT_FLOWLINES = "split_flowlines"

    # Complex operations - require full rebuild
    REMOVE_BRAIDED_RIVERS = "remove_braided_rivers"
    REMOVE_PARALLEL_PATHS = "remove_parallel_paths"
    REMOVE_CYCLES = "remove_cycles"
    BATCH_OPERATIONS = "batch_operations"

class GraphState(Enum):
    """Graph state tracking for lazy evaluation and caching."""
    CLEAN = "clean"                    # Graph is up-to-date
    TOPOLOGY_DIRTY = "topology_dirty"  # Topology changed, need to rebuild connections
    ATTRIBUTES_DIRTY = "attributes_dirty"  # Only attributes changed
    STRUCTURE_DIRTY = "structure_dirty"    # Major structural changes

class GraphUpdateStrategy(Enum):
    """Update strategy decisions based on operation complexity."""
    INCREMENTAL = "incremental"  # Use incremental updates
    CONDITIONAL = "conditional"  # Decide based on change size
    REBUILD = "rebuild"         # Full rebuild required

class GraphStateManager:
    """
    Manages graph state, change tracking, and update strategy decisions.

    This class implements the core state management system that tracks changes
    to the river network and makes intelligent decisions about when to use
    incremental updates versus full rebuilds.
    """

    def __init__(self):
        """Initialize the state manager."""
        self.state = GraphState.CLEAN
        self.last_update_time = time.time()
        self.change_history: List[Tuple[float, GraphChangeType, Dict]] = []
        self.pending_changes: List[Tuple[GraphChangeType, Dict]] = []

        # Cached expensive computations
        self._cached_cycles = None
        self._cached_parallel_paths = None
        self._cached_braided_channels = None
        self._cached_linear_segments = None
        self._cached_confluences = None

        # Cache validity flags
        self._cycles_valid = False
        self._parallel_paths_valid = False
        self._braided_channels_valid = False
        self._linear_segments_valid = False
        self._confluences_valid = False

        # Performance tracking
        self.update_stats = {
            'incremental_updates': 0,
            'full_rebuilds': 0,
            'total_update_time': 0.0,
            'last_rebuild_time': 0.0
        }

        # Thresholds for update strategy decisions
        self.incremental_threshold = 0.1  # 10% of network size
        self.rebuild_threshold = 0.3      # 30% of network size

    def record_change(self, change_type: GraphChangeType, metadata: Optional[Dict] = None):
        """
        Record a change to the network for tracking and strategy decisions.

        Args:
            change_type: Type of change being made
            metadata: Optional metadata about the change (e.g., number of flowlines affected)
        """
        timestamp = time.time()
        metadata = metadata or {}

        # Add to change history
        self.change_history.append((timestamp, change_type, metadata))

        # Add to pending changes if not immediately processed
        self.pending_changes.append((change_type, metadata))

        # Update state based on change type
        self._update_state_for_change(change_type)

        logger.debug(f"Recorded change: {change_type.value} with metadata: {metadata}")

    def _update_state_for_change(self, change_type: GraphChangeType):
        """Update graph state based on the type of change."""
        if change_type in [GraphChangeType.ADD_FLOWLINE, GraphChangeType.REMOVE_FLOWLINE]:
            if self.state == GraphState.CLEAN:
                self.state = GraphState.TOPOLOGY_DIRTY
        elif change_type == GraphChangeType.UPDATE_ATTRIBUTES:
            if self.state == GraphState.CLEAN:
                self.state = GraphState.ATTRIBUTES_DIRTY
        else:
            # Complex operations always mark as structure dirty
            self.state = GraphState.STRUCTURE_DIRTY

        # Invalidate relevant caches
        self._invalidate_caches_for_change(change_type)

    def _invalidate_caches_for_change(self, change_type: GraphChangeType):
        """Invalidate cached computations based on change type."""
        if change_type in [GraphChangeType.ADD_FLOWLINE, GraphChangeType.REMOVE_FLOWLINE,
                          GraphChangeType.REMOVE_BRAIDED_RIVERS, GraphChangeType.REMOVE_PARALLEL_PATHS,
                          GraphChangeType.REMOVE_CYCLES]:
            # Topology changes invalidate most caches
            self._cycles_valid = False
            self._parallel_paths_valid = False
            self._braided_channels_valid = False
            self._linear_segments_valid = False
            self._confluences_valid = False
        elif change_type == GraphChangeType.MERGE_LINEAR_SEGMENTS:
            self._linear_segments_valid = False
            self._confluences_valid = False
        elif change_type == GraphChangeType.UPDATE_ATTRIBUTES:
            # Attribute changes don't invalidate topology caches
            pass

    def get_update_strategy(self, network_size: int, change_metadata: Optional[Dict] = None) -> GraphUpdateStrategy:
        """
        Determine the optimal update strategy based on changes and network state.

        Args:
            network_size: Current number of flowlines in the network
            change_metadata: Optional metadata about the pending change

        Returns:
            Recommended update strategy
        """
        if not self.pending_changes:
            return GraphUpdateStrategy.INCREMENTAL

        # Analyze pending changes
        complex_changes = sum(1 for change_type, _ in self.pending_changes
                            if change_type in [GraphChangeType.REMOVE_BRAIDED_RIVERS,
                                             GraphChangeType.REMOVE_PARALLEL_PATHS,
                                             GraphChangeType.REMOVE_CYCLES])

        simple_changes = sum(1 for change_type, _ in self.pending_changes
                           if change_type in [GraphChangeType.ADD_FLOWLINE,
                                            GraphChangeType.REMOVE_FLOWLINE,
                                            GraphChangeType.UPDATE_ATTRIBUTES])

        # Decision logic
        if complex_changes > 0:
            return GraphUpdateStrategy.REBUILD

        if simple_changes > 0:
            # Check if the number of changes exceeds threshold
            change_ratio = simple_changes / max(network_size, 1)
            if change_ratio > self.rebuild_threshold:
                return GraphUpdateStrategy.REBUILD
            elif change_ratio > self.incremental_threshold:
                return GraphUpdateStrategy.CONDITIONAL
            else:
                return GraphUpdateStrategy.INCREMENTAL

        return GraphUpdateStrategy.CONDITIONAL

    def clear_pending_changes(self):
        """Clear pending changes after they have been processed."""
        self.pending_changes.clear()
        self.state = GraphState.CLEAN
        self.last_update_time = time.time()

    def get_cached_result(self, cache_type: str):
        """Get cached computation result if valid."""
        cache_map = {
            'cycles': (self._cached_cycles, self._cycles_valid),
            'parallel_paths': (self._cached_parallel_paths, self._parallel_paths_valid),
            'braided_channels': (self._cached_braided_channels, self._braided_channels_valid),
            'linear_segments': (self._cached_linear_segments, self._linear_segments_valid),
            'confluences': (self._cached_confluences, self._confluences_valid)
        }

        if cache_type in cache_map:
            cached_value, is_valid = cache_map[cache_type]
            if is_valid and cached_value is not None:
                logger.debug(f"Using cached result for {cache_type}")
                return cached_value

        return None

    def set_cached_result(self, cache_type: str, result):
        """Set cached computation result and mark as valid."""
        if cache_type == 'cycles':
            self._cached_cycles = result
            self._cycles_valid = True
        elif cache_type == 'parallel_paths':
            self._cached_parallel_paths = result
            self._parallel_paths_valid = True
        elif cache_type == 'braided_channels':
            self._cached_braided_channels = result
            self._braided_channels_valid = True
        elif cache_type == 'linear_segments':
            self._cached_linear_segments = result
            self._linear_segments_valid = True
        elif cache_type == 'confluences':
            self._cached_confluences = result
            self._confluences_valid = True

        logger.debug(f"Cached result for {cache_type}")

    def record_update_performance(self, strategy: GraphUpdateStrategy, duration: float):
        """Record performance metrics for update operations."""
        self.update_stats['total_update_time'] += duration

        if strategy == GraphUpdateStrategy.INCREMENTAL:
            self.update_stats['incremental_updates'] += 1
        else:
            self.update_stats['full_rebuilds'] += 1
            self.update_stats['last_rebuild_time'] = duration

    def get_performance_stats(self) -> Dict:
        """Get performance statistics for monitoring and optimization."""
        total_updates = (self.update_stats['incremental_updates'] +
                        self.update_stats['full_rebuilds'])

        if total_updates > 0:
            avg_update_time = self.update_stats['total_update_time'] / total_updates
            incremental_ratio = self.update_stats['incremental_updates'] / total_updates
        else:
            avg_update_time = 0.0
            incremental_ratio = 0.0

        return {
            'total_updates': total_updates,
            'incremental_updates': self.update_stats['incremental_updates'],
            'full_rebuilds': self.update_stats['full_rebuilds'],
            'incremental_ratio': incremental_ratio,
            'avg_update_time': avg_update_time,
            'last_rebuild_time': self.update_stats['last_rebuild_time'],
            'current_state': self.state.value
        }

class pyrivergraph:
    """
    Directed graph representation of a river network for braided channel analysis.

    This class models the flowline network as a directed graph where:
    - Nodes represent vertices (confluence/divergence points)
    - Edges represent flowlines connecting vertices
    - Multiple edges between same nodes represent braided channels

    Uses NetworkX when available for enhanced graph algorithms, falls back to
    custom implementation for compatibility.

    Enhanced with intelligent state management and hybrid update strategies for
    optimal performance with dynamic network modifications.
    """

    # ========================================================================
    # INITIALIZATION & BASIC GRAPH OPERATIONS
    # ========================================================================

    def __init__(self, flowlines: List[pyflowline], outlet_vertex: Optional[pyvertex] = None):
        """
        Initialize the river network graph from flowlines.

        Args:
            flowlines: List of flowline objects representing the river network
            outlet_vertex: Optional outlet vertex for the drainage network. If provided,
                          enables outlet-based operations and optimizations.
        """
        self.aFlowline = flowlines
        self.pVertex_outlet = outlet_vertex
        self.pVertex_outlet_id: Optional[int] = None
        self.vertex_to_id: Dict[pyvertex, int] = {}
        self.id_to_vertex: Dict[int, pyvertex] = {}
        self.aFlowline_edges: Dict[int, Tuple[int, int]] = {}
        self.aVertex: List[pyvertex] = []  # Will be populated during graph building

        # Always use custom implementation
        self.adjacency_list: DefaultDict[int, List[Tuple[int, int]]] = defaultdict(list)
        self.in_degree: DefaultDict[int, int] = defaultdict(int)
        self.out_degree: DefaultDict[int, int] = defaultdict(int)

        # Initialize state management system
        self.state_manager = GraphStateManager()

        logger.debug("Using custom graph implementation with state management")
        self._build_graph()
        # Set outlet vertex ID if outlet was provided
        if self.pVertex_outlet is not None:
            self._set_outlet_vertex_id()

        # Mark initial state as clean
        self.state_manager.clear_pending_changes()

    def get_sources(self) -> List[int]:
        """Get source nodes (headwaters) with no incoming edges."""
        return [node_id for node_id in self.id_to_vertex.keys() if self.in_degree[node_id] == 0]

    def get_sinks(self) -> List[int]:
        """Get sink nodes (outlets) with no outgoing edges."""
        return [node_id for node_id in self.id_to_vertex.keys() if self.out_degree[node_id] == 0]

    def get_vertices(self) -> List[pyvertex]:
        """
        Extract all unique vertices from the flowline network.

        This method provides functionality equivalent to find_flowline_vertex(),
        returning a list of unique vertices with assigned vertex IDs. The vertices
        are returned in the order they were encountered during graph construction.

        Returns:
            List[pyvertex]: List of unique vertices with assigned lVertexID values

        Example:
            >>> river_graph = pyrivergraph(flowlines)
            >>> vertices = river_graph.get_vertices()
            >>> print(f"Found {len(vertices)} unique vertices")
        """
        # Return the stored vertex list which maintains the order from graph construction
        logger.debug(f"Extracted {len(self.aVertex)} unique vertices from river graph")
        return self.aVertex.copy()  # Return a copy to prevent external modification

    def get_vertex_by_id(self, vertex_id: int) -> Optional[pyvertex]:
        """
        Get a vertex by its internal graph ID.

        Args:
            vertex_id (int): Internal vertex ID (0-based)

        Returns:
            Optional[pyvertex]: The vertex object, or None if not found
        """
        return self.id_to_vertex.get(vertex_id)

    def get_vertex_id(self, vertex: pyvertex) -> Optional[int]:
        """
        Get the internal graph ID for a vertex.

        Args:
            vertex (pyvertex): The vertex to look up

        Returns:
            Optional[int]: Internal vertex ID (0-based), or None if not found
        """
        return self.vertex_to_id.get(vertex)

    def get_vertex_count(self) -> int:
        """
        Get the total number of unique vertices in the network.

        Returns:
            int: Number of unique vertices
        """
        return len(self.id_to_vertex)

    # ========================================================================
    # NETWORK SIMPLIFICATION OPERATIONS
    # ========================================================================

    def remove_disconnected_flowlines(self, flowlines: List[pyflowline], outlet_vertex: Optional[pyvertex] = None) -> List[pyflowline]:
        """
        Remove flowlines that don't flow out to the specified outlet vertex.

        This method performs a backward traversal from the outlet vertex to identify
        all flowlines that are connected to the drainage network. Flowlines that
        cannot reach the outlet (isolated components, disconnected segments) are removed.

        Args:
            flowlines: List of input flowlines
            outlet_vertex: Optional outlet vertex. If not provided, uses the outlet vertex
                         from initialization. If neither is available, returns original flowlines.

        Returns:
            List of flowlines that are connected to the outlet vertex
        """
        if not flowlines:
            return flowlines

        # Use provided outlet vertex or fall back to stored one
        target_outlet = outlet_vertex if outlet_vertex is not None else self.pVertex_outlet

        if target_outlet is None:
            logger.warning("No outlet vertex provided and none stored during initialization")
            return flowlines

        # Use stored outlet vertex ID if it matches, otherwise find it
        outlet_vertex_id = None
        if target_outlet == self.pVertex_outlet and self.pVertex_outlet_id is not None:
            outlet_vertex_id = self.pVertex_outlet_id
        else:
            # Find the outlet vertex ID in our graph
            for vertex_id, vertex in self.id_to_vertex.items():
                if vertex == target_outlet:
                    outlet_vertex_id = vertex_id
                    break

        if outlet_vertex_id is None:
            logger.warning("Outlet vertex not found in network graph")
            return flowlines

        logger.info(f"Removing disconnected flowlines using outlet vertex {outlet_vertex_id}")

        # Use a more robust approach: forward traversal to validate complete drainage paths
        # Step 1: Find all vertices that can reach the outlet (backward traversal)
        outlet_reachable_vertices = self._find_outlet_reachable_vertices(outlet_vertex_id)

        # Step 2: For each flowline, check if it contributes to the drainage network
        # A flowline contributes if:
        # 1. Its start vertex can reach the outlet, AND
        # 2. Its end vertex can reach the outlet, AND
        # 3. The flowline itself is part of a path that leads to the outlet
        reachable_flowlines = set()

        for flowline_idx, (start_id, end_id) in self.aFlowline_edges.items():
            # Check if both vertices are in the outlet-reachable set
            if start_id in outlet_reachable_vertices and end_id in outlet_reachable_vertices:
                # Additional check: ensure this flowline is on a path to outlet
                if self._flowline_contributes_to_outlet(flowline_idx, start_id, end_id, outlet_vertex_id, outlet_reachable_vertices):
                    reachable_flowlines.add(flowline_idx)

        # Filter flowlines to keep only those that are reachable from outlet
        connected_flowlines = []
        disconnected_count = 0

        for i, flowline in enumerate(flowlines):
            if i in reachable_flowlines:
                connected_flowlines.append(flowline)
            else:
                disconnected_count += 1
                logger.debug(f"Removing disconnected flowline {i}: {flowline.pVertex_start} -> {flowline.pVertex_end}")

        logger.info(f"Removed {disconnected_count} disconnected flowlines, kept {len(connected_flowlines)} connected flowlines")
        return connected_flowlines

    def remove_braided_river(self) -> List[pyflowline]:
        """
        Remove braided channels from the river network.
        Enhanced with state management and change tracking.

        This method identifies braided channels (multiple flowlines between
        the same vertex pair) and removes redundant ones to create a simplified
        river network without braiding.

        Returns:
            List of flowlines with braided channels removed
        """
        # Record the operation
        self.state_manager.record_change(
            GraphChangeType.REMOVE_BRAIDED_RIVERS,
            {'network_size': len(self.aFlowline)}
        )

        braided_channels = self.find_braided_channels()
        if not braided_channels:
            logger.info("No braided channels found in the network")
            return self.aFlowline

        logger.info(f"Found {len(braided_channels)} braided channel groups to resolve")

        # Set to keep track of flowline indices to remove
        flowlines_to_remove = set()

        for (start_id, end_id), flowline_indices in braided_channels.items():
            # Keep the first flowline, remove the rest
            for flowline_idx in flowline_indices[1:]:
                flowlines_to_remove.add(flowline_idx)
                logger.debug(f"Marking flowline {flowline_idx} for removal (braided channel)")

        # Create new list of flowlines excluding those marked for removal
        simplified_flowlines = [
            flowline for idx, flowline in enumerate(self.aFlowline)
            if idx not in flowlines_to_remove
        ]

        logger.info(f"Removed {len(flowlines_to_remove)} braided flowlines, "
                    f"resulting in {len(simplified_flowlines)} flowlines")

        # Update the network using smart update system
        self._update_graph_flowlines(simplified_flowlines)

        return simplified_flowlines

    def remove_parallel_river(self) -> List[pyflowline]:
        """
        Remove parallel rivers using graph-based approach with class instance flowlines.
        Enhanced with state management and change tracking.

        This method replaces the standalone resolve_parallel_paths function by leveraging
        the graph structure to identify alternative routes between distant vertices and
        select the path with the highest cumulative hydrological significance.

        Parallel paths are alternative routes between the same start and end vertices.
        This method selects the most significant route based on cumulative path metrics.

        Returns:
            List[pyflowline]: Flowlines with parallel paths resolved

        Example:
            >>> river_graph = pyrivergraph(flowlines, outlet_vertex)
            >>> resolved_flowlines = river_graph.remove_parallel_river()
            >>> print(f"Resolved parallel paths: {len(flowlines)} -> {len(resolved_flowlines)} flowlines")
        """
        if not hasattr(self, 'aFlowline') or not self.aFlowline:
            logger.warning('No flowlines available in class instance for parallel river removal')
            return []

        if len(self.aFlowline) <= 1:
            logger.debug("Skipping parallel path resolution: insufficient flowlines")
            return self.aFlowline.copy()

        # Record the operation
        self.state_manager.record_change(
            GraphChangeType.REMOVE_PARALLEL_PATHS,
            {'network_size': len(self.aFlowline)}
        )

        logger.info(f"Removing parallel rivers from {len(self.aFlowline)} flowlines using graph-based approach")

        try:
            parallel_groups = self.find_parallel_paths()

            if not parallel_groups:
                logger.debug("No parallel paths found")
                return self.aFlowline.copy()

            logger.info(f"Found {len(parallel_groups)} parallel path groups")
        except Exception as e:
            logger.error(f"Error finding parallel paths: {e}")
            return self.aFlowline.copy()

        flowlines_to_remove = set()

        for group in parallel_groups:
            paths = group['paths']

            if len(paths) <= 1:
                continue

            # Calculate path significance for each complete path (multi-flowline groups)
            path_scores = []
            for path_idx, path_flowlines in enumerate(paths):
                # Calculate cumulative score for the entire path
                path_total_score = 0.0
                path_length = 0.0
                path_max_order = 0
                path_total_area = 0.0
                valid_flowlines = 0

                for flowline_idx in path_flowlines:
                    if flowline_idx >= len(self.aFlowline):
                        logger.warning(f"Invalid flowline index {flowline_idx}, skipping")
                        continue

                    flowline = self.aFlowline[flowline_idx]
                    valid_flowlines += 1

                    # Only use attributes if they have valid values (> 0)
                    stream_order = getattr(flowline, 'iStream_order', -1)
                    stream_order = stream_order if stream_order > 0 else 1

                    drainage_area = getattr(flowline, 'dDrainage_area', 0.0)
                    drainage_area = drainage_area if drainage_area > 0 else 0.0

                    length = getattr(flowline, 'dLength', 0.0)
                    length = length if length > 0 else 0.0

                    # Accumulate path metrics
                    path_max_order = max(path_max_order, stream_order)
                    path_total_area += drainage_area
                    path_length += length

                if valid_flowlines == 0:
                    continue

                if path_length > 0:
                    path_score = path_length
                else:
                    path_score = 1.0

                path_scores.append((path_score, path_idx, path_flowlines))

            if not path_scores:
                continue

            # Keep the highest scoring path, remove others
            path_scores.sort(reverse=True)
            best_score, best_path_idx, best_path_flowlines = path_scores[0]

            # Add flowlines from non-selected paths to removal set
            for score, path_idx, path_flowlines in path_scores[1:]:
                for flowline_idx in path_flowlines:
                    flowlines_to_remove.add(flowline_idx)

            logger.debug(f"Selected path {best_path_idx} with score {best_score:.2e} from {len(paths)} parallel paths "
                        f"(vertices {group['start_vertex']} -> {group['end_vertex']})")
            logger.debug(f"Keeping {len(best_path_flowlines)} flowlines, removing {sum(len(path_flowlines) for _, _, path_flowlines in path_scores[1:])} flowlines")

        # Return flowlines with parallel paths removed (filter by index)
        try:
            result = [self.aFlowline[i] for i in range(len(self.aFlowline)) if i not in flowlines_to_remove]
            logger.info(f"Removed {len(flowlines_to_remove)} flowlines to resolve parallel paths")

            # Validate result
            if len(result) == 0 and len(self.aFlowline) > 0:
                logger.warning("All flowlines were removed during parallel path resolution, returning original")
                return self.aFlowline.copy()

            # Update the network using smart update system
            self._update_graph_flowlines(result)

            return result
        except Exception as e:
            logger.error(f"Error filtering flowlines during parallel path resolution: {e}")
            return self.aFlowline.copy()

    def remove_cycle(self) -> List[pyflowline]:
        """
        Detect and break cycles in the network by removing lowest priority flowlines.
        Enhanced with state management and change tracking.

        This method replaces the standalone break_network_cycles function by leveraging
        the graph structure to detect cycles and remove the flowline with the lowest
        hydrological significance from each cycle.

        Uses depth-first search to detect cycles and removes the flowline with
        the lowest hydrological significance from each cycle.

        Returns:
            List[pyflowline]: Flowlines with cycles broken (acyclic network)

        Example:
            >>> river_graph = pyrivergraph(flowlines, outlet_vertex)
            >>> acyclic_flowlines = river_graph.remove_cycle()
            >>> print(f"Removed cycles: {len(flowlines) - len(acyclic_flowlines)} flowlines removed")
        """
        if not hasattr(self, 'aFlowline') or not self.aFlowline:
            logger.warning('No flowlines available in class instance for cycle removal')
            return []

        if len(self.aFlowline) <= 1:
            logger.debug("Skipping cycle removal: insufficient flowlines")
            return self.aFlowline.copy()

        # Record the operation
        self.state_manager.record_change(
            GraphChangeType.REMOVE_CYCLES,
            {'network_size': len(self.aFlowline)}
        )

        logger.info(f"Removing cycles from {len(self.aFlowline)} flowlines using graph-based approach")

        try:
            cycles = self.detect_cycles()

            if not cycles:
                logger.debug("No cycles detected in network")
                return self.aFlowline.copy()

            logger.info(f"Detected {len(cycles)} cycles in network")
        except Exception as e:
            logger.error(f"Error detecting cycles: {e}")
            return self.aFlowline.copy()

        flowlines_to_remove = set()

        for cycle_vertices in cycles:
            if len(cycle_vertices) < 3:  # Need at least 3 vertices for a meaningful cycle
                continue

            # Find flowlines involved in this cycle
            cycle_flowlines = []
            for i in range(len(cycle_vertices) - 1):
                start_vertex_id = cycle_vertices[i]
                end_vertex_id = cycle_vertices[i + 1]

                # Find flowlines connecting these vertices
                for flowline_idx, (s_id, e_id) in self.aFlowline_edges.items():
                    if s_id == start_vertex_id and e_id == end_vertex_id:
                        if flowline_idx < len(self.aFlowline):
                            cycle_flowlines.append((flowline_idx, self.aFlowline[flowline_idx]))

            if cycle_flowlines:
                # Remove the lowest priority flowline from the cycle
                def cycle_priority(flowline_data) -> Tuple:
                    flowline_idx, flowline = flowline_data
                    # Only use attributes if they have valid values (> 0)
                    stream_order = getattr(flowline, 'iStream_order', -1)
                    stream_order = stream_order if stream_order > 0 else 1

                    drainage_area = getattr(flowline, 'dDrainage_area', 0.0)
                    drainage_area = drainage_area if drainage_area > 0 else 0.0

                    length = getattr(flowline, 'dLength', 0.0)
                    length = length if length > 0 else 0.0

                    return (stream_order, drainage_area, length)

                worst_flowline_idx, worst_flowline = min(cycle_flowlines, key=cycle_priority)
                flowlines_to_remove.add(worst_flowline_idx)

                logger.debug(f"Removing flowline {worst_flowline_idx} (order={worst_flowline.iStream_order}) to break cycle")

        # Return flowlines with cycle-causing ones removed (filter by index)
        try:
            result = [self.aFlowline[i] for i in range(len(self.aFlowline)) if i not in flowlines_to_remove]
            logger.info(f"Removed {len(flowlines_to_remove)} flowlines to break cycles")

            # Validate result
            if len(result) == 0 and len(self.aFlowline) > 0:
                logger.warning("All flowlines were removed during cycle removal, returning original")
                return self.aFlowline.copy()

            # Update the network using smart update system
            self._update_graph_flowlines(result)

            return result
        except Exception as e:
            logger.error(f"Error filtering flowlines during cycle removal: {e}")
            return self.aFlowline.copy()

    def remove_small_rivers_iterative(self, dThreshold_small_river, nIterations=3, iFlag_debug=0, sWorkspace_output_basin=None):
        """
        Remove small rivers iteratively using graph-based approach with graph updates at each step.
        This method replaces the standalone function loop in basin.py for small river removal.

        Args:
            dThreshold_small_river (float): Length threshold for small river removal
            nIterations (int, optional): Number of iterations. Defaults to 3.
            iFlag_debug (int, optional): Debug flag for output files. Defaults to 0.
            sWorkspace_output_basin (str, optional): Output workspace for debug files. Defaults to None.

        Returns:
            list: Updated flowlines after iterative small river removal
        """
        import os
        import numpy as np
        from pyflowline.formats.export_flowline import export_flowline_to_geojson
        from pyflowline.formats.export_vertex import export_vertex_to_geojson

        # Start with current flowlines
        aFlowline_current = self.aFlowline.copy()

        for i in range(nIterations):
            sStep = "{:02d}".format(i+1)
            print(f'Iteration {sStep}: removing small rivers with threshold {dThreshold_small_river}')

            # Step 1: Remove small rivers using graph-based approach
            aFlowline_filtered = self._remove_small_rivers(aFlowline_current, dThreshold_small_river)

            if iFlag_debug == 1 and sWorkspace_output_basin is not None:
                sFilename_out = f'flowline_large_{sStep}_before_intersect.geojson'
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_geojson(aFlowline_filtered, sFilename_out)

            # Step 2: Update the graph with filtered flowlines
            self._update_graph_flowlines(aFlowline_filtered)

            # Step 3: Update stream order using graph-based method
            aFlowline_filtered = self.update_head_water_stream_order()

            # Step 4: Use graph-based merge instead of standalone functions
            # This replaces find_flowline_confluence + merge_flowline
            aFlowline_merged = self.merge_linear_segments(aFlowline_filtered)

            if iFlag_debug == 1 and sWorkspace_output_basin is not None:
                sFilename_out = f'flowline_merge_{sStep}_before_intersect.geojson'
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_geojson(aFlowline_merged, sFilename_out)

            # Step 5: Update graph again with merged flowlines
            self._update_graph_flowlines(aFlowline_merged)

            # Step 6: Final stream order update for this iteration
            aFlowline_current = self.update_head_water_stream_order()

            # Break if only one flowline remains
            if len(aFlowline_current) == 1:
                break

        return aFlowline_current

    # ========================================================================
    # NETWORK ANALYSIS & DETECTION
    # ========================================================================

    def find_braided_channels(self) -> Dict[Tuple[int, int], List[int]]:
        """
        Find braided channels (multiple flowlines between same vertex pair).
        Uses caching for performance optimization.

        Returns:
            Dictionary mapping (start_vertex_id, end_vertex_id) to list of flowline indices
        """
        # Check cache first
        cached_result = self.state_manager.get_cached_result('braided_channels')
        if cached_result is not None:
            return cached_result

        channel_groups: DefaultDict[Tuple[int, int], List[int]] = defaultdict(list)

        for flowline_idx, (start_id, end_id) in self.aFlowline_edges.items():
            channel_groups[(start_id, end_id)].append(flowline_idx)

        # Return only groups with multiple channels (braided)
        result = {key: indices for key, indices in channel_groups.items() if len(indices) > 1}

        # Cache the result
        self.state_manager.set_cached_result('braided_channels', result)
        return result

    def find_parallel_paths(self) -> List[Dict]:
        """
        Find parallel paths between vertices (alternative routes).
        Uses caching for performance optimization.

        Parallel paths are defined as multiple different paths that have the
        SAME starting vertex AND the SAME ending vertex.

        Returns:
            List of parallel path group dictionaries, each containing:
            - 'paths': List of paths, where each path is a list of flowline indices
            - 'start_vertex': Start vertex ID
            - 'end_vertex': End vertex ID
        """
        # Check cache first
        cached_result = self.state_manager.get_cached_result('parallel_paths')
        if cached_result is not None:
            return cached_result

        parallel_groups = []
        visited_pairs = set()

        # Create lists to avoid "dictionary changed size during iteration" error
        start_ids = list(self.adjacency_list.keys())
        target_ids = list(self.adjacency_list.keys())

        for start_id in start_ids:
            for target_id in target_ids:
                if start_id == target_id or (start_id, target_id) in visited_pairs:
                    continue

                # Find ALL paths from start_id to target_id (same start, same end)
                paths = self._find_all_paths(start_id, target_id, max_depth=10)

                # Only consider as parallel if we have multiple different paths
                # between the SAME start and end vertices
                if len(paths) > 1:
                    # Convert paths to flowline indices - preserve path structure
                    path_flowline_groups = []
                    for path in paths:
                        path_flowlines = self._path_to_flowlines(path)
                        if path_flowlines and len(path_flowlines) > 0:
                            path_flowline_groups.append(path_flowlines)

                    # Ensure we have multiple valid alternative paths between same start/end
                    if len(path_flowline_groups) > 1:
                        # Verify paths are truly different (not just same path found multiple times)
                        unique_paths = []
                        for path_flowlines in path_flowline_groups:
                            path_signature = tuple(sorted(path_flowlines))
                            if path_signature not in [tuple(sorted(up)) for up in unique_paths]:
                                unique_paths.append(path_flowlines)

                        # Only add if we have multiple unique paths between same vertices
                        if len(unique_paths) > 1:
                            logger.debug(f"Found {len(unique_paths)} parallel paths between vertices {start_id} -> {target_id}")
                            parallel_groups.append({
                                'paths': unique_paths,  # Each path is a list of flowline indices
                                'start_vertex': start_id,
                                'end_vertex': target_id
                            })

                visited_pairs.add((start_id, target_id))

        # Cache the result
        self.state_manager.set_cached_result('parallel_paths', parallel_groups)
        return parallel_groups

    def detect_cycles(self) -> List[List[int]]:
        """
        Custom cycle detection using DFS with recursion stack.
        Uses caching for performance optimization.
        """
        # Check cache first
        cached_result = self.state_manager.get_cached_result('cycles')
        if cached_result is not None:
            return cached_result

        visited = set()
        rec_stack = set()
        cycles = []

        def dfs_cycle_detection(node_id: int, path: List[int]) -> bool:
            try:
                visited.add(node_id)
                rec_stack.add(node_id)

                # Get neighbors from adjacency list
                neighbors = list(self.adjacency_list[node_id])

                for neighbor_id, _ in neighbors:
                    try:
                        if neighbor_id in rec_stack:
                            # Found a cycle - back edge to a node in recursion stack
                            try:
                                cycle_start_idx = path.index(neighbor_id)
                                cycle = path[cycle_start_idx:] + [neighbor_id]
                                cycles.append(cycle)
                                logger.debug(f"Detected cycle: {cycle}")
                            except ValueError:
                                cycle = [node_id, neighbor_id]
                                cycles.append(cycle)
                                logger.debug(f"Detected simple back-edge cycle: {cycle}")
                            except Exception as e:
                                logger.warning(f"Error processing cycle from {node_id} to {neighbor_id}: {e}")
                                cycles.append([node_id, neighbor_id])
                            return True
                        elif neighbor_id not in visited:
                            try:
                                new_path = path + [neighbor_id]
                                if dfs_cycle_detection(neighbor_id, new_path):
                                    return True
                            except RecursionError:
                                logger.error(f"Recursion limit reached during cycle detection at node {neighbor_id}")
                                if len(path) > 1:
                                    cycles.append(path + [neighbor_id])
                                return False
                            except Exception as e:
                                logger.error(f"Error in recursive cycle detection for neighbor {neighbor_id}: {e}")
                                return False
                    except Exception as e:
                        logger.error(f"Error processing neighbor {neighbor_id} from node {node_id}: {e}")
                        continue

                try:
                    rec_stack.remove(node_id)
                except KeyError:
                    logger.warning(f"Node {node_id} not found in recursion stack during removal")
                return False

            except Exception as e:
                logger.error(f"Critical error in cycle detection for node {node_id}: {e}")
                try:
                    rec_stack.discard(node_id)
                except Exception:
                    pass
                return False

        try:
            # Get node list from adjacency list
            node_ids = list(self.adjacency_list.keys())

            for node_id in node_ids:
                if node_id not in visited:
                    try:
                        dfs_cycle_detection(node_id, [node_id])
                    except Exception as e:
                        logger.error(f"Error starting cycle detection from node {node_id}: {e}")
                        continue
        except Exception as e:
            logger.error(f"Critical error during cycle detection initialization: {e}")

        logger.info(f"Custom cycle detection completed. Found {len(cycles)} cycles.")

        # Cache the result
        self.state_manager.set_cached_result('cycles', cycles)
        return cycles

    def find_linear_segments(self) -> List[List[int]]:
        """
        Find linear segments that can be merged into single flowlines.
        Uses caching for performance optimization.

        A linear segment is a chain of flowlines where each intermediate vertex
        has exactly one incoming and one outgoing edge (degree 2 vertices).
        These represent artificially segmented flowlines that should be merged.

        Returns:
            List of linear segments, where each segment is a list of flowline indices
            that should be merged together in order from upstream to downstream.
        """
        # Check cache first
        cached_result = self.state_manager.get_cached_result('linear_segments')
        if cached_result is not None:
            return cached_result

        linear_segments = []
        processed_flowlines = set()

        # Find all vertices with degree 2 (one in, one out) - these are merge candidates
        merge_vertices = []
        for vertex_id in self.id_to_vertex.keys():
            if self.in_degree[vertex_id] == 1 and self.out_degree[vertex_id] == 1:
                merge_vertices.append(vertex_id)

        logger.debug(f"Found {len(merge_vertices)} degree-2 vertices for potential merging")

        # For each flowline, try to build a linear segment starting from it
        for flowline_idx, (start_id, end_id) in self.aFlowline_edges.items():
            if flowline_idx in processed_flowlines:
                continue

            # Start building a segment from this flowline
            segment = [flowline_idx]
            current_end = end_id

            # Extend downstream while we encounter degree-2 vertices
            while current_end in merge_vertices:
                # Find the next flowline from this vertex
                next_flowlines = [fl_idx for neighbor_id, fl_idx in self.adjacency_list[current_end]]

                if len(next_flowlines) == 1:
                    next_flowline_idx = next_flowlines[0]
                    if next_flowline_idx not in processed_flowlines:
                        segment.append(next_flowline_idx)
                        # Update current_end to the end of the next flowline
                        if next_flowline_idx in self.aFlowline_edges:
                            _, current_end = self.aFlowline_edges[next_flowline_idx]
                        else:
                            break
                    else:
                        break
                else:
                    break

            # Only consider segments with multiple flowlines
            if len(segment) > 1:
                linear_segments.append(segment)
                processed_flowlines.update(segment)
                logger.debug(f"Found linear segment with {len(segment)} flowlines: {segment}")

        logger.info(f"Found {len(linear_segments)} linear segments for potential merging")

        # Cache the result
        self.state_manager.set_cached_result('linear_segments', linear_segments)
        return linear_segments

    def find_outlet_connected_components(self, flowlines: List[pyflowline], outlet_vertex: Optional[pyvertex] = None) -> Dict[str, List[int]]:
        """
        Find connected components in the network and classify them relative to the outlet.

        This method identifies different connected components in the river network
        and classifies them as either connected to the main drainage network (outlet)
        or as isolated components that don't contribute to the main flow.

        Args:
            flowlines: List of flowlines to analyze
            outlet_vertex: Optional outlet vertex. If not provided, uses stored outlet.

        Returns:
            Dictionary with keys:
            - 'main_network': List of flowline indices connected to outlet
            - 'isolated_components': List of lists, each containing flowline indices
              for an isolated component
            - 'component_stats': Dictionary with statistics about components
        """
        if not flowlines:
            return {
                'main_network': [],
                'isolated_components': [],
                'component_stats': {'total_components': 0, 'main_network_size': 0, 'isolated_count': 0}
            }

        # Use provided outlet vertex or fall back to stored one
        target_outlet = outlet_vertex if outlet_vertex is not None else self.pVertex_outlet

        # Build a temporary graph for this analysis
        temp_adjacency = defaultdict(list)
        temp_edges = {}

        for i, flowline in enumerate(flowlines):
            start_vertex = flowline.pVertex_start
            end_vertex = flowline.pVertex_end

            # Create temporary vertex mapping
            start_id = id(start_vertex)
            end_id = id(end_vertex)

            temp_edges[i] = (start_id, end_id)
            temp_adjacency[start_id].append((end_id, i))

        # Find all connected components using DFS
        visited_vertices = set()
        visited_flowlines = set()
        components = []

        def dfs_component(vertex_id, current_component_flowlines):
            if vertex_id in visited_vertices:
                return

            visited_vertices.add(vertex_id)

            # Explore all edges from this vertex
            for neighbor_id, flowline_idx in temp_adjacency[vertex_id]:
                if flowline_idx not in visited_flowlines:
                    visited_flowlines.add(flowline_idx)
                    current_component_flowlines.append(flowline_idx)
                    dfs_component(neighbor_id, current_component_flowlines)

        # Find components by starting DFS from each unvisited vertex
        for flowline_idx, (start_id, end_id) in temp_edges.items():
            if flowline_idx not in visited_flowlines:
                component_flowlines = []
                dfs_component(start_id, component_flowlines)
                if component_flowlines:
                    components.append(component_flowlines)

        # Identify which component contains the outlet (if provided)
        main_network = []
        isolated_components = []

        if target_outlet is not None:
            outlet_id = id(target_outlet)

            # Find which component contains the outlet
            for component in components:
                contains_outlet = False
                for flowline_idx in component:
                    flowline = flowlines[flowline_idx]
                    if id(flowline.pVertex_start) == outlet_id or id(flowline.pVertex_end) == outlet_id:
                        contains_outlet = True
                        break

                if contains_outlet:
                    main_network = component
                else:
                    isolated_components.append(component)
        else:
            # Without outlet, treat largest component as main network
            if components:
                main_network = max(components, key=len)
                isolated_components = [comp for comp in components if comp != main_network]

        # Calculate statistics
        component_stats = {
            'total_components': len(components),
            'main_network_size': len(main_network),
            'isolated_count': len(isolated_components),
            'isolated_sizes': [len(comp) for comp in isolated_components]
        }

        logger.info(f"Found {len(components)} connected components: "
                   f"main network ({len(main_network)} flowlines), "
                   f"{len(isolated_components)} isolated components")

        return {
            'main_network': main_network,
            'isolated_components': isolated_components,
            'component_stats': component_stats
        }

    # ========================================================================
    # NETWORK MODIFICATION & PROCESSING
    # ========================================================================

    def merge_linear_segments(self, flowlines: List[pyflowline]) -> List[pyflowline]:
        """
        Merge linear segments of flowlines into single flowlines.
        Enhanced with state management and change tracking.

        This method identifies chains of flowlines that can be merged into single
        flowlines (linear segments where intermediate vertices have degree 2) and
        creates new merged flowlines to replace them.

        Args:
            flowlines: List of flowlines to process

        Returns:
            List of flowlines with linear segments merged
        """
        if not flowlines:
            return flowlines

        # Record the operation
        self.state_manager.record_change(
            GraphChangeType.MERGE_LINEAR_SEGMENTS,
            {'network_size': len(flowlines)}
        )

        # Update graph with current flowlines for analysis
        original_flowlines = self.aFlowline
        self.aFlowline = flowlines
        self._build_graph()

        try:
            linear_segments = self.find_linear_segments()

            if not linear_segments:
                logger.info("No linear segments found for merging")
                return flowlines

            logger.info(f"Found {len(linear_segments)} linear segments to merge")

            # Create new flowlines by merging segments
            merged_flowlines = []
            processed_indices = set()

            for segment in linear_segments:
                if len(segment) < 2:
                    continue

                try:
                    # Get flowlines in the segment
                    segment_flowlines = [flowlines[idx] for idx in segment if idx < len(flowlines)]

                    if len(segment_flowlines) < 2:
                        continue

                    # Create merged flowline from the segment
                    merged_flowline = self._merge_flowline_segment(segment_flowlines)
                    if merged_flowline is not None:
                        merged_flowlines.append(merged_flowline)
                        processed_indices.update(segment)
                        logger.debug(f"Merged {len(segment_flowlines)} flowlines into single flowline")

                except Exception as e:
                    logger.error(f"Error merging segment {segment}: {e}")
                    continue

            # Add non-merged flowlines
            for i, flowline in enumerate(flowlines):
                if i not in processed_indices:
                    merged_flowlines.append(flowline)

            logger.info(f"Merged {len(linear_segments)} segments, "
                       f"resulting in {len(merged_flowlines)} flowlines "
                       f"(reduced from {len(flowlines)})")

            return merged_flowlines

        except Exception as e:
            logger.error(f"Error during linear segment merging: {e}")
            return flowlines
        finally:
            # Restore original flowlines
            self.aFlowline = original_flowlines
            self._build_graph()

    # ========================================================================
    # STREAM ANALYSIS & TOPOLOGY
    # ========================================================================

    def update_head_water_stream_order(self) -> List[pyflowline]:
        """
        Update stream order for head water flowlines using graph-based approach.
        Enhanced with state management and change tracking.

        This method replaces the standalone function by leveraging the graph structure
        to identify headwater flowlines and update their stream order based on
        network topology and confluence analysis.

        Returns:
            List[pyflowline]: Updated flowlines with corrected stream orders
        """
        if not hasattr(self, 'aFlowline') or not self.aFlowline:
            logger.warning('No flowlines available in class instance for stream order update')
            return []

        # Record the operation
        self.state_manager.record_change(
            GraphChangeType.UPDATE_ATTRIBUTES,
            {'network_size': len(self.aFlowline), 'operation': 'stream_order_update'}
        )

        logger.info(f"Updating stream order for {len(self.aFlowline)} flowlines using graph-based approach")

        try:
            # Step 1: Identify headwater flowlines (sources in the graph)
            headwater_flowlines = self.identify_headwater_flowlines()

            # Step 2: Perform topological sort to process flowlines in correct order
            sorted_flowlines = self._topological_sort_flowlines()

            # Step 3: Update stream orders based on confluence rules
            updated_flowlines = []
            for flowline in sorted_flowlines:
                updated_order = self._calculate_stream_order(flowline)
                flowline.iStream_order = updated_order
                updated_flowlines.append(flowline)

            logger.info(f"Updated stream orders for {len(updated_flowlines)} flowlines")
            return updated_flowlines

        except Exception as e:
            logger.error(f"Error updating stream orders: {e}")
            return self.aFlowline.copy()

    def identify_headwater_flowlines(self) -> List[pyflowline]:
        """
        Identify headwater flowlines (sources) in the network.

        Returns:
            List of flowlines that are headwaters (no upstream connections)
        """
        headwater_flowlines = []

        for i, flowline in enumerate(self.aFlowline):
            start_vertex_id = self.vertex_to_id.get(flowline.pVertex_start)
            if start_vertex_id is not None and self.in_degree[start_vertex_id] == 0:
                headwater_flowlines.append(flowline)

        logger.debug(f"Identified {len(headwater_flowlines)} headwater flowlines")
        return headwater_flowlines

    def define_stream_topology(self) -> Dict[str, any]:
        """
        Define comprehensive stream topology including confluences and stream segments.
        Enhanced with state management and caching.

        Returns:
            Dictionary containing topology information:
            - confluences: List of confluence objects
            - stream_segments: List of stream segment definitions
            - topology_stats: Statistics about the network topology
        """
        # Check cache first
        cached_result = self.state_manager.get_cached_result('confluences')
        if cached_result is not None:
            return cached_result

        confluences = []
        stream_segments = []

        try:
            # Find all confluence points (vertices with multiple incoming edges)
            for vertex_id, vertex in self.id_to_vertex.items():
                if self.in_degree[vertex_id] > 1:  # Confluence point
                    confluence = self._create_confluence_object(vertex_id, vertex)
                    if confluence:
                        confluences.append(confluence)

            # Define stream segments between confluences
            stream_segments = self._define_stream_segments(confluences)

            # Calculate topology statistics
            topology_stats = {
                'total_confluences': len(confluences),
                'total_segments': len(stream_segments),
                'max_stream_order': max((fl.iStream_order for fl in self.aFlowline), default=0),
                'network_complexity': len(confluences) / max(len(self.aFlowline), 1)
            }

            result = {
                'confluences': confluences,
                'stream_segments': stream_segments,
                'topology_stats': topology_stats
            }

            # Cache the result
            self.state_manager.set_cached_result('confluences', result)

            logger.info(f"Defined stream topology: {len(confluences)} confluences, {len(stream_segments)} segments")
            return result

        except Exception as e:
            logger.error(f"Error defining stream topology: {e}")
            return {
                'confluences': [],
                'stream_segments': [],
                'topology_stats': {'total_confluences': 0, 'total_segments': 0, 'max_stream_order': 0, 'network_complexity': 0}
            }

    def define_stream_order(self, iFlag_so_method_in: int = 1) -> Dict[str, any]:
        """
        Define stream order for all flowlines using graph-based approach with confluence analysis.
        Enhanced with state management and caching.

        This method adapts the reference implementation to work with the graph structure,
        supporting both Strahler (default) and Shreve stream ordering methods.

        Args:
            iFlag_so_method_in: Stream ordering method (1=Strahler, 2=Shreve)

        Returns:
            Dictionary containing:
            - 'flowlines': Updated flowlines with stream orders
            - 'stream_orders': Array of stream order values
            - 'confluences': List of confluence objects used
            - 'method': Stream ordering method used
            - 'statistics': Statistics about stream order distribution
        """
        if not hasattr(self, 'aFlowline') or not self.aFlowline:
            logger.warning('No flowlines available for stream order calculation')
            return {
                'flowlines': [],
                'stream_orders': np.array([]),
                'confluences': [],
                'method': iFlag_so_method_in,
                'statistics': {'total_flowlines': 0, 'max_order': 0, 'headwater_count': 0}
            }

        # Record the operation
        self.state_manager.record_change(
            GraphChangeType.UPDATE_ATTRIBUTES,
            {'network_size': len(self.aFlowline), 'operation': 'define_stream_order', 'method': iFlag_so_method_in}
        )

        logger.info(f"Defining stream order for {len(self.aFlowline)} flowlines using "
                   f"{'Strahler' if iFlag_so_method_in == 1 else 'Shreve'} method")

        try:
            # Step 1: Initialize headwater flowlines (stream order = 1)
            self._initialize_headwater_stream_orders()

            # Step 2: Get confluence information from graph topology
            confluences = self._extract_confluences_for_stream_order()

            if not confluences:
                logger.info("No confluences found - all flowlines are headwaters")
                return self._create_stream_order_result(iFlag_so_method_in, confluences)

            # Step 3: Build spatial index for efficient confluence lookup
            confluence_index = self._build_confluence_spatial_index(confluences)

            # Step 4: Process confluences iteratively until all stream orders are defined
            self._process_confluences_iteratively(confluences, confluence_index, iFlag_so_method_in)

            # Step 5: Validate and return results
            result = self._create_stream_order_result(iFlag_so_method_in, confluences)

            logger.info(f"Stream order calculation completed: max order = {result['statistics']['max_order']}, "
                       f"headwater count = {result['statistics']['headwater_count']}")

            return result

        except Exception as e:
            logger.error(f"Error defining stream order: {e}")
            return {
                'flowlines': self.aFlowline.copy(),
                'stream_orders': np.array([getattr(fl, 'iStream_order', 1) for fl in self.aFlowline]),
                'confluences': [],
                'method': iFlag_so_method_in,
                'statistics': {'total_flowlines': len(self.aFlowline), 'max_order': 1, 'headwater_count': 0},
                'error': str(e)
            }

    # ========================================================================
    # INCREMENTAL UPDATE OPERATIONS
    # ========================================================================

    def add_flowline_incremental(self, new_flowline: pyflowline) -> bool:
        """
        Add a single flowline to the network using incremental updates.

        Args:
            new_flowline: The flowline to add to the network

        Returns:
            bool: True if successfully added, False otherwise
        """
        try:
            # Record the operation
            self.state_manager.record_change(
                GraphChangeType.ADD_FLOWLINE,
                {'flowline_id': getattr(new_flowline, 'lFlowlineID', 'unknown')}
            )

            # Add to flowline list
            self.aFlowline.append(new_flowline)
            new_flowline_idx = len(self.aFlowline) - 1

            # Update vertex mappings
            start_vertex = new_flowline.pVertex_start
            end_vertex = new_flowline.pVertex_end

            # Add vertices if they don't exist
            start_id = self._add_vertex_if_new(start_vertex)
            end_id = self._add_vertex_if_new(end_vertex)

            # Update graph structure
            self.aFlowline_edges[new_flowline_idx] = (start_id, end_id)
            self.adjacency_list[start_id].append((end_id, new_flowline_idx))
            self.out_degree[start_id] += 1
            self.in_degree[end_id] += 1

            logger.debug(f"Added flowline {new_flowline_idx} incrementally")
            return True

        except Exception as e:
            logger.error(f"Error adding flowline incrementally: {e}")
            return False

    def remove_flowline_incremental(self, flowline_idx: int) -> bool:
        """
        Remove a single flowline from the network using incremental updates.

        Args:
            flowline_idx: Index of the flowline to remove

        Returns:
            bool: True if successfully removed, False otherwise
        """
        try:
            if flowline_idx >= len(self.aFlowline) or flowline_idx < 0:
                logger.warning(f"Invalid flowline index {flowline_idx}")
                return False

            # Record the operation
            flowline = self.aFlowline[flowline_idx]
            self.state_manager.record_change(
                GraphChangeType.REMOVE_FLOWLINE,
                {'flowline_id': getattr(flowline, 'lFlowlineID', 'unknown')}
            )

            # Get edge information
            if flowline_idx in self.aFlowline_edges:
                start_id, end_id = self.aFlowline_edges[flowline_idx]

                # Update adjacency list
                self.adjacency_list[start_id] = [
                    (neighbor, fl_idx) for neighbor, fl_idx in self.adjacency_list[start_id]
                    if fl_idx != flowline_idx
                ]

                # Update degrees
                self.out_degree[start_id] -= 1
                self.in_degree[end_id] -= 1

                # Remove edge mapping
                del self.aFlowline_edges[flowline_idx]

            # Remove from flowline list (mark as None to preserve indices)
            self.aFlowline[flowline_idx] = None

            logger.debug(f"Removed flowline {flowline_idx} incrementally")
            return True

        except Exception as e:
            logger.error(f"Error removing flowline incrementally: {e}")
            return False

    def update_flowline_attributes_incremental(self, flowline_idx: int, attributes: Dict) -> bool:
        """
        Update attributes of a single flowline using incremental updates.

        Args:
            flowline_idx: Index of the flowline to update
            attributes: Dictionary of attributes to update

        Returns:
            bool: True if successfully updated, False otherwise
        """
        try:
            if flowline_idx >= len(self.aFlowline) or flowline_idx < 0:
                logger.warning(f"Invalid flowline index {flowline_idx}")
                return False

            flowline = self.aFlowline[flowline_idx]
            if flowline is None:
                logger.warning(f"Flowline {flowline_idx} has been removed")
                return False

            # Record the operation
            self.state_manager.record_change(
                GraphChangeType.UPDATE_ATTRIBUTES,
                {'flowline_id': getattr(flowline, 'lFlowlineID', 'unknown'), 'attributes': list(attributes.keys())}
            )

            # Update attributes
            for attr_name, attr_value in attributes.items():
                if hasattr(flowline, attr_name):
                    setattr(flowline, attr_name, attr_value)
                else:
                    logger.warning(f"Flowline does not have attribute {attr_name}")

            logger.debug(f"Updated attributes for flowline {flowline_idx}")
            return True

        except Exception as e:
            logger.error(f"Error updating flowline attributes: {e}")
            return False

    # ========================================================================
    # VALIDATION & PERFORMANCE
    # ========================================================================

    def validate_graph_consistency(self) -> Dict[str, any]:
        """
        Validate the consistency of the graph structure and data.

        Returns:
            Dictionary containing validation results and any issues found
        """
        validation_results = {
            'is_valid': True,
            'issues': [],
            'warnings': [],
            'statistics': {}
        }

        try:
            # Check vertex consistency
            vertex_issues = self._validate_vertices()
            validation_results['issues'].extend(vertex_issues)

            # Check edge consistency
            edge_issues = self._validate_edges()
            validation_results['issues'].extend(edge_issues)

            # Check flowline consistency
            flowline_issues = self._validate_flowlines()
            validation_results['issues'].extend(flowline_issues)

            # Check degree consistency
            degree_issues = self._validate_degrees()
            validation_results['issues'].extend(degree_issues)

            # Calculate statistics
            validation_results['statistics'] = {
                'total_vertices': len(self.id_to_vertex),
                'total_edges': len(self.aFlowline_edges),
                'total_flowlines': len([fl for fl in self.aFlowline if fl is not None]),
                'isolated_vertices': len([v_id for v_id in self.id_to_vertex.keys()
                                        if self.in_degree[v_id] == 0 and self.out_degree[v_id] == 0])
            }

            # Set overall validity
            validation_results['is_valid'] = len(validation_results['issues']) == 0

            logger.info(f"Graph validation completed: {'VALID' if validation_results['is_valid'] else 'INVALID'} "
                       f"({len(validation_results['issues'])} issues, {len(validation_results['warnings'])} warnings)")

        except Exception as e:
            validation_results['is_valid'] = False
            validation_results['issues'].append(f"Validation error: {e}")
            logger.error(f"Error during graph validation: {e}")

        return validation_results

    def get_update_performance_stats(self) -> Dict[str, any]:
        """
        Get performance statistics for the update system.

        Returns:
            Dictionary containing performance metrics and recommendations
        """
        base_stats = self.state_manager.get_performance_stats()

        # Add additional analysis
        network_size = len(self.aFlowline)
        cache_hit_ratio = self._calculate_cache_hit_ratio()

        performance_analysis = {
            'network_size': network_size,
            'cache_hit_ratio': cache_hit_ratio,
            'recommended_strategy': self._get_recommended_strategy(base_stats, network_size),
            'optimization_suggestions': self._get_optimization_suggestions(base_stats, network_size)
        }

        return {**base_stats, **performance_analysis}

    def benchmark_update_strategies(self, test_operations: List[Dict]) -> Dict[str, any]:
        """
        Benchmark different update strategies with test operations.

        Args:
            test_operations: List of test operations to benchmark

        Returns:
            Dictionary containing benchmark results
        """
        benchmark_results = {
            'incremental_times': [],
            'rebuild_times': [],
            'operation_results': []
        }

        try:
            for operation in test_operations:
                # Test incremental approach
                start_time = time.time()
                incremental_result = self._execute_test_operation(operation, strategy='incremental')
                incremental_time = time.time() - start_time
                benchmark_results['incremental_times'].append(incremental_time)

                # Test rebuild approach
                start_time = time.time()
                rebuild_result = self._execute_test_operation(operation, strategy='rebuild')
                rebuild_time = time.time() - start_time
                benchmark_results['rebuild_times'].append(rebuild_time)

                benchmark_results['operation_results'].append({
                    'operation': operation,
                    'incremental_time': incremental_time,
                    'rebuild_time': rebuild_time,
                    'speedup': rebuild_time / incremental_time if incremental_time > 0 else float('inf')
                })

            # Calculate summary statistics
            benchmark_results['summary'] = {
                'avg_incremental_time': np.mean(benchmark_results['incremental_times']),
                'avg_rebuild_time': np.mean(benchmark_results['rebuild_times']),
                'avg_speedup': np.mean([r['speedup'] for r in benchmark_results['operation_results']]),
                'total_operations': len(test_operations)
            }

        except Exception as e:
            logger.error(f"Error during benchmarking: {e}")
            benchmark_results['error'] = str(e)

        return benchmark_results

    # ========================================================================
    # CACHE & STATE MANAGEMENT
    # ========================================================================

    def invalidate_cache(self, cache_types: Optional[List[str]] = None):
        """
        Invalidate specific caches or all caches.

        Args:
            cache_types: List of cache types to invalidate, or None for all
        """
        if cache_types is None:
            # Invalidate all caches
            self.state_manager._cycles_valid = False
            self.state_manager._parallel_paths_valid = False
            self.state_manager._braided_channels_valid = False
            self.state_manager._linear_segments_valid = False
            self.state_manager._confluences_valid = False
            logger.debug("Invalidated all caches")
        else:
            # Invalidate specific caches
            for cache_type in cache_types:
                if cache_type == 'cycles':
                    self.state_manager._cycles_valid = False
                elif cache_type == 'parallel_paths':
                    self.state_manager._parallel_paths_valid = False
                elif cache_type == 'braided_channels':
                    self.state_manager._braided_channels_valid = False
                elif cache_type == 'linear_segments':
                    self.state_manager._linear_segments_valid = False
                elif cache_type == 'confluences':
                    self.state_manager._confluences_valid = False
                logger.debug(f"Invalidated {cache_type} cache")

    def get_cache_status(self) -> Dict[str, bool]:
        """
        Get the current status of all caches.

        Returns:
            Dictionary mapping cache names to their validity status
        """
        return {
            'cycles': self.state_manager._cycles_valid,
            'parallel_paths': self.state_manager._parallel_paths_valid,
            'braided_channels': self.state_manager._braided_channels_valid,
            'linear_segments': self.state_manager._linear_segments_valid,
            'confluences': self.state_manager._confluences_valid
        }

    def record_operation(self, operation_type: str, metadata: Optional[Dict] = None):
        """
        Record an operation for performance tracking and state management.

        Args:
            operation_type: Type of operation performed
            metadata: Optional metadata about the operation
        """
        try:
            change_type = GraphChangeType(operation_type)
            self.state_manager.record_change(change_type, metadata)
        except ValueError:
            logger.warning(f"Unknown operation type: {operation_type}")

    def get_network_state(self) -> Dict[str, any]:
        """
        Get comprehensive information about the current network state.

        Returns:
            Dictionary containing network state information
        """
        return {
            'graph_state': self.state_manager.state.value,
            'network_size': len(self.aFlowline),
            'vertex_count': len(self.id_to_vertex),
            'edge_count': len(self.aFlowline_edges),
            'cache_status': self.get_cache_status(),
            'performance_stats': self.get_update_performance_stats(),
            'last_update': self.state_manager.last_update_time,
            'pending_changes': len(self.state_manager.pending_changes)
        }

    def create_batch_operation(self) -> 'BatchOperation':
        """
        Create a new batch operation manager.

        Returns:
            BatchOperation: New batch operation manager
        """
        return BatchOperation(self)

    def execute_batch_operations(self, operations: List[Dict]) -> bool:
        """
        Execute a list of operations as a batch.

        Args:
            operations: List of operation dictionaries

        Returns:
            bool: True if all operations succeeded, False otherwise
        """
        batch = self.create_batch_operation()

        for op in operations:
            if op['type'] == 'add_flowline':
                batch.add_flowline(op['flowline'])
            elif op['type'] == 'remove_flowline':
                batch.remove_flowline(op['flowline_idx'])
            elif op['type'] == 'update_attributes':
                batch.update_attributes(op['flowline_idx'], op['attributes'])

        return batch.execute()

    # ========================================================================
    # PRIVATE METHODS - GRAPH CONSTRUCTION & MANAGEMENT
    # ========================================================================

    def _build_graph(self):
        """
        Build the graph structure from flowlines.
        This is the core method that constructs the adjacency list representation.
        """
        # Clear existing graph data
        self.vertex_to_id.clear()
        self.id_to_vertex.clear()
        self.aFlowline_edges.clear()
        self.adjacency_list.clear()
        self.in_degree.clear()
        self.out_degree.clear()
        self.aVertex.clear()

        vertex_counter = 0

        # Process each flowline to build vertex mappings and edges
        for flowline_idx, flowline in enumerate(self.aFlowline):
            if flowline is None:  # Skip removed flowlines
                continue

            start_vertex = flowline.pVertex_start
            end_vertex = flowline.pVertex_end

            # Add start vertex if not seen before
            if start_vertex not in self.vertex_to_id:
                self.vertex_to_id[start_vertex] = vertex_counter
                self.id_to_vertex[vertex_counter] = start_vertex
                self.aVertex.append(start_vertex)
                start_vertex.lVertexID = vertex_counter  # Set vertex ID
                vertex_counter += 1

            # Add end vertex if not seen before
            if end_vertex not in self.vertex_to_id:
                self.vertex_to_id[end_vertex] = vertex_counter
                self.id_to_vertex[vertex_counter] = end_vertex
                self.aVertex.append(end_vertex)
                end_vertex.lVertexID = vertex_counter  # Set vertex ID
                vertex_counter += 1

            # Get vertex IDs
            start_id = self.vertex_to_id[start_vertex]
            end_id = self.vertex_to_id[end_vertex]

            # Store edge mapping
            self.aFlowline_edges[flowline_idx] = (start_id, end_id)

            # Add to adjacency list
            self.adjacency_list[start_id].append((end_id, flowline_idx))

            # Update degrees
            self.out_degree[start_id] += 1
            self.in_degree[end_id] += 1

        logger.debug(f"Built graph with {len(self.id_to_vertex)} vertices and {len(self.aFlowline_edges)} edges")

    def _set_outlet_vertex_id(self):
        """Set the outlet vertex ID if outlet vertex is provided."""
        if self.pVertex_outlet is not None:
            self.pVertex_outlet_id = self.vertex_to_id.get(self.pVertex_outlet)
            if self.pVertex_outlet_id is not None:
                logger.debug(f"Set outlet vertex ID: {self.pVertex_outlet_id}")
            else:
                logger.warning("Outlet vertex not found in graph")

    def _update_graph_flowlines(self, new_flowlines: List[pyflowline]):
        """
        Update the graph with a new set of flowlines using smart update strategy.

        Args:
            new_flowlines: New list of flowlines to update the graph with
        """
        # Determine update strategy
        strategy = self.state_manager.get_update_strategy(len(self.aFlowline))

        start_time = time.time()

        if strategy == GraphUpdateStrategy.INCREMENTAL and len(new_flowlines) > 0:
            # Try incremental update if feasible
            success = self._try_incremental_update(new_flowlines)
            if not success:
                # Fall back to rebuild
                self.aFlowline = new_flowlines
                self._build_graph()
        else:
            # Full rebuild
            self.aFlowline = new_flowlines
            self._build_graph()

        duration = time.time() - start_time
        self.state_manager.record_update_performance(strategy, duration)
        self.state_manager.clear_pending_changes()

    def _try_incremental_update(self, new_flowlines: List[pyflowline]) -> bool:
        """
        Attempt to update the graph incrementally.

        Args:
            new_flowlines: New flowlines to incorporate

        Returns:
            bool: True if incremental update succeeded, False if rebuild needed
        """
        try:
            # For now, use simple heuristic: if size difference is small, try incremental
            size_diff = abs(len(new_flowlines) - len(self.aFlowline))
            if size_diff > len(self.aFlowline) * 0.1:  # More than 10% change
                return False

            # Update flowline list
            self.aFlowline = new_flowlines

            # Rebuild graph (for now - could be optimized further)
            self._build_graph()

            return True

        except Exception as e:
            logger.error(f"Incremental update failed: {e}")
            return False

    def _add_vertex_if_new(self, vertex: pyvertex) -> int:
        """
        Add vertex to graph if it doesn't exist, return its ID.

        Args:
            vertex: Vertex to add

        Returns:
            int: Vertex ID
        """
        if vertex not in self.vertex_to_id:
            vertex_id = len(self.id_to_vertex)
            self.vertex_to_id[vertex] = vertex_id
            self.id_to_vertex[vertex_id] = vertex
            self.aVertex.append(vertex)
            vertex.lVertexID = vertex_id
            return vertex_id
        else:
            return self.vertex_to_id[vertex]

    # ========================================================================
    # PRIVATE METHODS - PATH FINDING & ANALYSIS
    # ========================================================================

    def _find_all_paths(self, start_id: int, target_id: int, max_depth: int = 10) -> List[List[int]]:
        """
        Find all paths from start to target vertex using DFS.

        Args:
            start_id: Starting vertex ID
            target_id: Target vertex ID
            max_depth: Maximum search depth to prevent infinite loops

        Returns:
            List of paths, where each path is a list of vertex IDs
        """
        paths = []

        def dfs_paths(current_id: int, target_id: int, path: List[int], visited: Set[int], depth: int):
            if depth > max_depth:
                return

            if current_id == target_id:
                paths.append(path.copy())
                return

            if current_id in visited:
                return

            visited.add(current_id)

            # Explore neighbors
            for neighbor_id, _ in self.adjacency_list[current_id]:
                if neighbor_id not in visited:
                    path.append(neighbor_id)
                    dfs_paths(neighbor_id, target_id, path, visited, depth + 1)
                    path.pop()

            visited.remove(current_id)

        if start_id != target_id:
            dfs_paths(start_id, target_id, [start_id], set(), 0)

        return paths

    def _path_to_flowlines(self, path: List[int]) -> List[int]:
        """
        Convert a path of vertex IDs to flowline indices.

        Args:
            path: List of vertex IDs representing a path

        Returns:
            List of flowline indices
        """
        flowline_indices = []

        for i in range(len(path) - 1):
            start_id = path[i]
            end_id = path[i + 1]

            # Find flowline connecting these vertices
            for neighbor_id, flowline_idx in self.adjacency_list[start_id]:
                if neighbor_id == end_id:
                    flowline_indices.append(flowline_idx)
                    break

        return flowline_indices

    def _find_outlet_reachable_vertices(self, outlet_vertex_id: int) -> Set[int]:
        """
        Find all vertices that can reach the outlet using backward traversal.

        Args:
            outlet_vertex_id: ID of the outlet vertex

        Returns:
            Set of vertex IDs that can reach the outlet
        """
        reachable = set()
        queue = deque([outlet_vertex_id])
        reachable.add(outlet_vertex_id)

        # Build reverse adjacency list for backward traversal
        reverse_adjacency = defaultdict(list)
        for start_id, neighbors in self.adjacency_list.items():
            for end_id, flowline_idx in neighbors:
                reverse_adjacency[end_id].append((start_id, flowline_idx))

        # Backward BFS from outlet
        while queue:
            current_id = queue.popleft()

            # Add all vertices that flow into current vertex
            for upstream_id, _ in reverse_adjacency[current_id]:
                if upstream_id not in reachable:
                    reachable.add(upstream_id)
                    queue.append(upstream_id)

        return reachable

    def _flowline_contributes_to_outlet(self, flowline_idx: int, start_id: int, end_id: int,
                                       outlet_vertex_id: int, outlet_reachable: Set[int]) -> bool:
        """
        Check if a flowline contributes to the drainage network leading to outlet.

        Args:
            flowline_idx: Index of the flowline
            start_id: Start vertex ID
            end_id: End vertex ID
            outlet_vertex_id: Outlet vertex ID
            outlet_reachable: Set of vertices reachable from outlet

        Returns:
            bool: True if flowline contributes to outlet drainage
        """
        # Simple check: both vertices must be reachable from outlet
        return start_id in outlet_reachable and end_id in outlet_reachable

    # ========================================================================
    # PRIVATE METHODS - NETWORK PROCESSING HELPERS
    # ========================================================================

    def _remove_small_rivers(self, flowlines: List[pyflowline], threshold: float) -> List[pyflowline]:
        """
        Remove small rivers based on length threshold.

        Args:
            flowlines: List of flowlines to filter
            threshold: Length threshold for removal

        Returns:
            List of flowlines with small rivers removed
        """
        filtered_flowlines = []

        for flowline in flowlines:
            length = getattr(flowline, 'dLength', 0.0)
            if length >= threshold:
                filtered_flowlines.append(flowline)
            else:
                logger.debug(f"Removing small river with length {length}")

        logger.info(f"Removed {len(flowlines) - len(filtered_flowlines)} small rivers")
        return filtered_flowlines

    def _topological_sort_flowlines(self) -> List[pyflowline]:
        """
        Perform topological sort of flowlines for stream order calculation.

        Returns:
            List of flowlines in topological order
        """
        # Simple implementation - could be optimized
        sorted_flowlines = []
        remaining_flowlines = self.aFlowline.copy()

        while remaining_flowlines:
            # Find flowlines with no upstream dependencies
            headwater_found = False
            for i, flowline in enumerate(remaining_flowlines):
                if flowline is None:
                    continue

                start_vertex_id = self.vertex_to_id.get(flowline.pVertex_start)
                if start_vertex_id is not None and self.in_degree[start_vertex_id] == 0:
                    sorted_flowlines.append(flowline)
                    remaining_flowlines[i] = None
                    headwater_found = True

            if not headwater_found:
                # Add remaining flowlines to avoid infinite loop
                for flowline in remaining_flowlines:
                    if flowline is not None:
                        sorted_flowlines.append(flowline)
                break

            # Remove None entries
            remaining_flowlines = [fl for fl in remaining_flowlines if fl is not None]

        return sorted_flowlines

    def _calculate_stream_order(self, flowline: pyflowline) -> int:
        """
        Calculate stream order for a flowline based on upstream confluences.

        Args:
            flowline: Flowline to calculate order for

        Returns:
            int: Stream order
        """
        start_vertex_id = self.vertex_to_id.get(flowline.pVertex_start)
        if start_vertex_id is None:
            return 1

        # If no upstream connections, it's order 1 (headwater)
        if self.in_degree[start_vertex_id] == 0:
            return 1

        # Calculate based on upstream orders (simplified)
        upstream_orders = []
        # This would need more complex logic to properly calculate stream orders
        # For now, return existing order or 1
        return getattr(flowline, 'iStream_order', 1)

    def _merge_flowline_segment(self, segment_flowlines: List[pyflowline]) -> Optional[pyflowline]:
        """
        Merge a segment of flowlines into a single flowline.

        Args:
            segment_flowlines: List of flowlines to merge

        Returns:
            Optional[pyflowline]: Merged flowline or None if merge failed
        """
        if not segment_flowlines:
            return None

        try:
            # Use first flowline as base
            merged_flowline = segment_flowlines[0]

            # Update end vertex to the end of the last flowline
            merged_flowline.pVertex_end = segment_flowlines[-1].pVertex_end

            # Accumulate length and other properties
            total_length = sum(getattr(fl, 'dLength', 0.0) for fl in segment_flowlines)
            merged_flowline.dLength = total_length

            # Use maximum stream order
            max_order = max(getattr(fl, 'iStream_order', 1) for fl in segment_flowlines)
            merged_flowline.iStream_order = max_order

            # Use maximum drainage area
            max_area = max(getattr(fl, 'dDrainage_area', 0.0) for fl in segment_flowlines)
            merged_flowline.dDrainage_area = max_area

            return merged_flowline

        except Exception as e:
            logger.error(f"Error merging flowline segment: {e}")
            return None

    # ========================================================================
    # PRIVATE METHODS - VERTEX & CONFLUENCE MANAGEMENT
    # ========================================================================

    def _create_confluence_object(self, vertex_id: int, vertex: pyvertex) -> Optional[pyconfluence]:
        """
        Create a confluence object for a vertex.

        Args:
            vertex_id: Vertex ID
            vertex: Vertex object

        Returns:
            Optional[pyconfluence]: Confluence object or None if creation failed
        """
        try:
            confluence = pyconfluence()
            confluence.pVertex = vertex
            confluence.lVertexID = vertex_id
            confluence.nIncoming = self.in_degree[vertex_id]
            confluence.nOutgoing = self.out_degree[vertex_id]

            # Add incoming flowlines
            confluence.aIncoming_flowline = []
            for start_id, neighbors in self.adjacency_list.items():
                for end_id, flowline_idx in neighbors:
                    if end_id == vertex_id and flowline_idx < len(self.aFlowline):
                        confluence.aIncoming_flowline.append(self.aFlowline[flowline_idx])

            # Add outgoing flowlines
            confluence.aOutgoing_flowline = []
            for end_id, flowline_idx in self.adjacency_list[vertex_id]:
                if flowline_idx < len(self.aFlowline):
                    confluence.aOutgoing_flowline.append(self.aFlowline[flowline_idx])

            return confluence

        except Exception as e:
            logger.error(f"Error creating confluence object: {e}")
            return None

    def _define_stream_segments(self, confluences: List[pyconfluence]) -> List[Dict]:
        """
        Define stream segments between confluences.

        Args:
            confluences: List of confluence objects

        Returns:
            List of stream segment definitions
        """
        segments = []

        try:
            # Simple implementation - could be enhanced
            for i, confluence in enumerate(confluences):
                segment = {
                    'segment_id': i,
                    'start_confluence': confluence,
                    'flowlines': confluence.aOutgoing_flowline,
                    'length': sum(getattr(fl, 'dLength', 0.0) for fl in confluence.aOutgoing_flowline)
                }
                segments.append(segment)

        except Exception as e:
            logger.error(f"Error defining stream segments: {e}")

        return segments

    # ========================================================================
    # PRIVATE METHODS - VALIDATION HELPERS
    # ========================================================================

    def _validate_vertices(self) -> List[str]:
        """Validate vertex consistency."""
        issues = []

        # Check vertex mappings consistency
        if len(self.vertex_to_id) != len(self.id_to_vertex):
            issues.append("Vertex mapping size mismatch")

        # Check vertex IDs are consistent
        for vertex, vertex_id in self.vertex_to_id.items():
            if vertex_id not in self.id_to_vertex:
                issues.append(f"Vertex ID {vertex_id} missing from id_to_vertex mapping")
            elif self.id_to_vertex[vertex_id] != vertex:
                issues.append(f"Vertex mapping inconsistency for ID {vertex_id}")

        return issues

    def _validate_edges(self) -> List[str]:
        """Validate edge consistency."""
        issues = []

        # Check edge mappings
        for flowline_idx, (start_id, end_id) in self.aFlowline_edges.items():
            if start_id not in self.id_to_vertex:
                issues.append(f"Edge {flowline_idx} references invalid start vertex {start_id}")
            if end_id not in self.id_to_vertex:
                issues.append(f"Edge {flowline_idx} references invalid end vertex {end_id}")

        return issues

    def _validate_flowlines(self) -> List[str]:
        """Validate flowline consistency."""
        issues = []

        # Check flowline indices
        for flowline_idx in self.aFlowline_edges.keys():
            if flowline_idx >= len(self.aFlowline):
                issues.append(f"Edge references invalid flowline index {flowline_idx}")
            elif self.aFlowline[flowline_idx] is None:
                issues.append(f"Edge references removed flowline {flowline_idx}")

        return issues

    def _validate_degrees(self) -> List[str]:
        """Validate degree consistency."""
        issues = []

        # Recalculate degrees and compare
        calculated_in_degree = defaultdict(int)
        calculated_out_degree = defaultdict(int)

        for start_id, neighbors in self.adjacency_list.items():
            calculated_out_degree[start_id] = len(neighbors)
            for end_id, _ in neighbors:
                calculated_in_degree[end_id] += 1

        # Check consistency
        for vertex_id in self.id_to_vertex.keys():
            if self.in_degree[vertex_id] != calculated_in_degree[vertex_id]:
                issues.append(f"In-degree mismatch for vertex {vertex_id}")
            if self.out_degree[vertex_id] != calculated_out_degree[vertex_id]:
                issues.append(f"Out-degree mismatch for vertex {vertex_id}")

        return issues

    # ========================================================================
    # PRIVATE METHODS - PERFORMANCE & OPTIMIZATION HELPERS
    # ========================================================================

    def _calculate_cache_hit_ratio(self) -> float:
        """Calculate cache hit ratio for performance analysis."""
        # Simple implementation - could be enhanced with actual hit/miss tracking
        cache_status = self.get_cache_status()
        valid_caches = sum(1 for valid in cache_status.values() if valid)
        total_caches = len(cache_status)
        return valid_caches / total_caches if total_caches > 0 else 0.0

    def _get_recommended_strategy(self, stats: Dict, network_size: int) -> str:
        """Get recommended update strategy based on performance stats."""
        if stats['incremental_ratio'] > 0.8 and stats['avg_update_time'] < 1.0:
            return "Continue using incremental updates"
        elif network_size > 10000:
            return "Consider batch operations for large networks"
        else:
            return "Mixed strategy based on operation type"

    def _get_optimization_suggestions(self, stats: Dict, network_size: int) -> List[str]:
        """Get optimization suggestions based on performance analysis."""
        suggestions = []

        if stats['incremental_ratio'] < 0.5:
            suggestions.append("Consider optimizing for more incremental updates")

        if stats['avg_update_time'] > 5.0:
            suggestions.append("Update times are high - consider caching optimizations")

        if network_size > 50000:
            suggestions.append("Large network detected - consider spatial indexing")

        return suggestions

    def _execute_test_operation(self, operation: Dict, strategy: str) -> bool:
        """Execute a test operation for benchmarking."""
        # Simplified implementation for testing
        try:
            if strategy == 'incremental':
                # Simulate incremental operation
                time.sleep(0.001)  # Simulate work
            else:
                # Simulate rebuild operation
                time.sleep(0.01)   # Simulate more work
            return True
        except Exception:
            return False

    # ========================================================================
    # PRIVATE METHODS - STREAM ORDER CALCULATION HELPERS
    # ========================================================================

    def _initialize_headwater_stream_orders(self):
        """
        Initialize stream orders for headwater flowlines (stream order = 1).
        Adapted from the reference update_head_water_stream_order function.
        """
        # Get all start and end vertices
        start_vertices = {flowline.pVertex_start for flowline in self.aFlowline if flowline is not None}
        end_vertices = {flowline.pVertex_end for flowline in self.aFlowline if flowline is not None}

        headwater_count = 0
        for flowline in self.aFlowline:
            if flowline is None:
                continue

            pVertex_start = flowline.pVertex_start
            # A flowline is headwater if its start vertex is not an end vertex of any other flowline
            is_headwater = pVertex_start in start_vertices and pVertex_start not in end_vertices

            if is_headwater:
                flowline.iStream_order = 1
                headwater_count += 1
            else:
                flowline.iStream_order = -1  # Mark as unprocessed

        logger.debug(f"Initialized {headwater_count} headwater flowlines with stream order 1")

    def _extract_confluences_for_stream_order(self) -> List[Dict]:
        """
        Extract confluence information from the graph topology for stream order calculation.

        Returns:
            List of confluence dictionaries with upstream/downstream flowline information
        """
        confluences = []

        # Find all vertices with multiple incoming edges (confluences)
        for vertex_id, vertex in self.id_to_vertex.items():
            if self.in_degree[vertex_id] > 1:  # Confluence point
                # Get upstream flowlines (flowing into this vertex)
                upstream_flowlines = []
                for start_id, neighbors in self.adjacency_list.items():
                    for end_id, flowline_idx in neighbors:
                        if end_id == vertex_id and flowline_idx < len(self.aFlowline):
                            flowline = self.aFlowline[flowline_idx]
                            if flowline is not None:
                                upstream_flowlines.append(flowline)

                # Get downstream flowline (flowing out of this vertex)
                downstream_flowline = None
                if vertex_id in self.adjacency_list:
                    for end_id, flowline_idx in self.adjacency_list[vertex_id]:
                        if flowline_idx < len(self.aFlowline):
                            downstream_flowline = self.aFlowline[flowline_idx]
                            break  # Take the first one (should be only one for proper confluence)

                if upstream_flowlines and downstream_flowline:
                    confluence = {
                        'vertex': vertex,
                        'vertex_id': vertex_id,
                        'upstream_flowlines': upstream_flowlines,
                        'downstream_flowline': downstream_flowline,
                        'treated': False
                    }
                    confluences.append(confluence)

        logger.debug(f"Extracted {len(confluences)} confluences from graph topology")
        return confluences

    def _build_confluence_spatial_index(self, confluences: List[Dict]):
        """
        Build spatial index for efficient confluence lookup.

        Args:
            confluences: List of confluence dictionaries

        Returns:
            Spatial index object (R-tree if available, otherwise fallback)
        """
        if HAS_RTREE:
            # Use R-tree for efficient spatial indexing
            index = RTreeindex()
            for i, confluence in enumerate(confluences):
                vertex = confluence['vertex']
                x, y = vertex.dLongitude_degree, vertex.dLatitude_degree
                # Create bounding box with small tolerance
                pBound = (x - 1E-5, y - 1E-5, x + 1E-5, y + 1E-5)
                index.insert(i, pBound)
            return index
        else:
            # Fallback to simple list-based lookup
            logger.debug("Using fallback spatial indexing (no R-tree available)")
            return confluences

    def _process_confluences_iteratively(self, confluences: List[Dict], confluence_index, iFlag_so_method_in: int):
        """
        Process confluences iteratively until all stream orders are defined.
        Adapted from the reference define_stream_order function.

        Args:
            confluences: List of confluence dictionaries
            confluence_index: Spatial index for confluence lookup
            iFlag_so_method_in: Stream ordering method (1=Strahler, 2=Shreve)
        """
        max_iterations = len(self.aFlowline) * 2  # Prevent infinite loops
        iteration = 0

        # Continue until all flowlines have positive stream orders
        while any(fl.iStream_order < 0 for fl in self.aFlowline if fl is not None) and iteration < max_iterations:
            iteration += 1
            progress_made = False

            for i, confluence in enumerate(confluences):
                if confluence['treated']:
                    continue

                upstream_flowlines = confluence['upstream_flowlines']
                downstream_flowline = confluence['downstream_flowline']

                if downstream_flowline is None:
                    continue

                # Check if all upstream flowlines have been processed (positive stream order)
                upstream_orders = [fl.iStream_order for fl in upstream_flowlines if fl.iStream_order >= 1]

                if len(upstream_orders) == len(upstream_flowlines):
                    # All upstream flowlines processed - calculate downstream order
                    confluence['treated'] = True
                    progress_made = True

                    # Calculate stream order based on method
                    if iFlag_so_method_in == 1:  # Strahler stream order
                        if len(set(upstream_orders)) > 1:
                            # Different orders - take maximum
                            stream_order = max(upstream_orders)
                        else:
                            # Same orders - increment by 1
                            stream_order = upstream_orders[0] + 1
                    else:  # Shreve stream order
                        stream_order = sum(upstream_orders) if upstream_orders else 1

                    # Update downstream flowline
                    downstream_flowline.iStream_order = stream_order

                    # Update any other flowlines that share the same stream segment
                    self._update_connected_flowlines_stream_order(
                        downstream_flowline, stream_order, confluences, confluence_index
                    )

                    logger.debug(f"Processed confluence {i}: assigned order {stream_order} to downstream flowline")

            if not progress_made:
                logger.warning(f"No progress made in iteration {iteration} - breaking to avoid infinite loop")
                break

        logger.info(f"Completed confluence processing in {iteration} iterations")

    def _update_connected_flowlines_stream_order(self, reference_flowline: pyflowline, stream_order: int,
                                                confluences: List[Dict], confluence_index):
        """
        Update stream order for flowlines connected to the same stream segment.
        Adapted from the reference implementation's confluence update logic.

        Args:
            reference_flowline: The flowline whose order was just set
            stream_order: The stream order to propagate
            confluences: List of confluence dictionaries
            confluence_index: Spatial index for confluence lookup
        """
        try:
            # Get the end vertex of the reference flowline
            end_vertex = reference_flowline.pVertex_end
            x, y = end_vertex.dLongitude_degree, end_vertex.dLatitude_degree

            # Find confluences near this vertex
            if HAS_RTREE and hasattr(confluence_index, 'intersection'):
                # Use R-tree spatial query
                delta = 1E-5
                bbox = (x - delta, y - delta, x + delta, y + delta)
                nearby_indices = list(confluence_index.intersection(bbox))
            else:
                # Fallback to linear search
                nearby_indices = []
                for i, confluence in enumerate(confluences):
                    conf_vertex = confluence['vertex']
                    conf_x, conf_y = conf_vertex.dLongitude_degree, conf_vertex.dLatitude_degree
                    if abs(conf_x - x) < 1E-5 and abs(conf_y - y) < 1E-5:
                        nearby_indices.append(i)

            # Update upstream flowlines in nearby confluences
            for conf_idx in nearby_indices:
                if conf_idx < len(confluences):
                    confluence = confluences[conf_idx]
                    if confluence['vertex'] == end_vertex:
                        # Update upstream flowlines that belong to the same stream segment
                        for upstream_flowline in confluence['upstream_flowlines']:
                            if (hasattr(upstream_flowline, 'iStream_segment') and
                                hasattr(reference_flowline, 'iStream_segment') and
                                upstream_flowline.iStream_segment == reference_flowline.iStream_segment):
                                upstream_flowline.iStream_order = stream_order

        except Exception as e:
            logger.warning(f"Error updating connected flowlines stream order: {e}")

    def _create_stream_order_result(self, iFlag_so_method_in: int, confluences: List[Dict]) -> Dict[str, any]:
        """
        Create the result dictionary for stream order calculation.

        Args:
            iFlag_so_method_in: Stream ordering method used
            confluences: List of confluence dictionaries

        Returns:
            Dictionary containing results and statistics
        """
        # Collect stream orders
        stream_orders = []
        headwater_count = 0
        max_order = 0

        for flowline in self.aFlowline:
            if flowline is not None:
                order = max(1, flowline.iStream_order)  # Ensure minimum order of 1
                stream_orders.append(order)
                max_order = max(max_order, order)
                if order == 1:
                    headwater_count += 1
            else:
                stream_orders.append(1)  # Default for removed flowlines

        # Create statistics
        statistics = {
            'total_flowlines': len([fl for fl in self.aFlowline if fl is not None]),
            'max_order': max_order,
            'headwater_count': headwater_count,
            'confluence_count': len(confluences),
            'method_name': 'Strahler' if iFlag_so_method_in == 1 else 'Shreve'
        }

        return {
            'flowlines': [fl for fl in self.aFlowline if fl is not None],
            'stream_orders': np.array(stream_orders),
            'confluences': confluences,
            'method': iFlag_so_method_in,
            'statistics': statistics
        }



class BatchOperation:
    """
    Batch operation manager for grouping multiple network modifications.

    This class allows multiple operations to be queued and executed together,
    enabling more efficient updates and better performance for bulk changes.
    """

    def __init__(self, river_graph: 'pyrivergraph'):
        """
        Initialize batch operation manager.

        Args:
            river_graph: The pyrivergraph instance to operate on
        """
        self.river_graph = river_graph
        self.operations = []
        self.executed = False

    def add_flowline(self, flowline: pyflowline):
        """
        Queue a flowline addition operation.

        Args:
            flowline: Flowline to add
        """
        if self.executed:
            raise RuntimeError("Cannot add operations to executed batch")

        self.operations.append({
            'type': 'add_flowline',
            'flowline': flowline
        })

    def remove_flowline(self, flowline_idx: int):
        """
        Queue a flowline removal operation.

        Args:
            flowline_idx: Index of flowline to remove
        """
        if self.executed:
            raise RuntimeError("Cannot add operations to executed batch")

        self.operations.append({
            'type': 'remove_flowline',
            'flowline_idx': flowline_idx
        })

    def update_attributes(self, flowline_idx: int, attributes: Dict):
        """
        Queue a flowline attribute update operation.

        Args:
            flowline_idx: Index of flowline to update
            attributes: Dictionary of attributes to update
        """
        if self.executed:
            raise RuntimeError("Cannot add operations to executed batch")

        self.operations.append({
            'type': 'update_attributes',
            'flowline_idx': flowline_idx,
            'attributes': attributes
        })

    def execute(self) -> bool:
        """
        Execute all queued operations as a batch.

        Returns:
            bool: True if all operations succeeded, False otherwise
        """
        if self.executed:
            raise RuntimeError("Batch has already been executed")

        if not self.operations:
            logger.warning("No operations to execute in batch")
            return True

        # Record batch operation
        self.river_graph.state_manager.record_change(
            GraphChangeType.BATCH_OPERATIONS,
            {'operation_count': len(self.operations)}
        )

        start_time = time.time()
        success = True

        try:
            # Execute operations based on optimal strategy
            strategy = self.river_graph.state_manager.get_update_strategy(
                len(self.river_graph.aFlowline),
                {'batch_size': len(self.operations)}
            )

            if strategy == GraphUpdateStrategy.INCREMENTAL:
                # Execute operations incrementally
                for operation in self.operations:
                    if not self._execute_single_operation(operation):
                        success = False
                        break
            else:
                # Execute as batch with rebuild
                success = self._execute_batch_rebuild()

            duration = time.time() - start_time
            self.river_graph.state_manager.record_update_performance(strategy, duration)

        except Exception as e:
            logger.error(f"Error executing batch operations: {e}")
            success = False
        finally:
            self.executed = True
            self.river_graph.state_manager.clear_pending_changes()

        logger.info(f"Executed batch of {len(self.operations)} operations: {'SUCCESS' if success else 'FAILED'}")
        return success

    def _execute_single_operation(self, operation: Dict) -> bool:
        """Execute a single operation incrementally."""
        try:
            if operation['type'] == 'add_flowline':
                return self.river_graph.add_flowline_incremental(operation['flowline'])
            elif operation['type'] == 'remove_flowline':
                return self.river_graph.remove_flowline_incremental(operation['flowline_idx'])
            elif operation['type'] == 'update_attributes':
                return self.river_graph.update_flowline_attributes_incremental(
                    operation['flowline_idx'], operation['attributes']
                )
            else:
                logger.warning(f"Unknown operation type: {operation['type']}")
                return False
        except Exception as e:
            logger.error(f"Error executing operation {operation['type']}: {e}")
            return False

    def _execute_batch_rebuild(self) -> bool:
        """Execute all operations with a full rebuild."""
        try:
            # Apply all operations to flowline list
            modified_flowlines = self.river_graph.aFlowline.copy()

            for operation in self.operations:
                if operation['type'] == 'add_flowline':
                    modified_flowlines.append(operation['flowline'])
                elif operation['type'] == 'remove_flowline':
                    idx = operation['flowline_idx']
                    if 0 <= idx < len(modified_flowlines):
                        modified_flowlines[idx] = None
                elif operation['type'] == 'update_attributes':
                    idx = operation['flowline_idx']
                    if 0 <= idx < len(modified_flowlines) and modified_flowlines[idx] is not None:
                        flowline = modified_flowlines[idx]
                        for attr_name, attr_value in operation['attributes'].items():
                            if hasattr(flowline, attr_name):
                                setattr(flowline, attr_name, attr_value)

            # Filter out None entries
            modified_flowlines = [fl for fl in modified_flowlines if fl is not None]

            # Rebuild graph with modified flowlines
            self.river_graph.aFlowline = modified_flowlines
            self.river_graph._build_graph()

            return True

        except Exception as e:
            logger.error(f"Error in batch rebuild: {e}")
            return False

    def get_operation_count(self) -> int:
        """Get the number of queued operations."""
        return len(self.operations)

    def clear(self):
        """Clear all queued operations."""
        if self.executed:
            raise RuntimeError("Cannot clear executed batch")
        self.operations.clear()
