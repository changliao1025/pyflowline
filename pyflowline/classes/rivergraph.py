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

# Set up logger
logger = logging.getLogger(__name__)
# Import for spatial indexing (R-tree)
try:
    from rtree.index import Index as RTreeindex
    HAS_RTREE = True
except ImportError:
    HAS_RTREE = False
    logger.warning("R-tree not available - using fallback spatial indexing")

import importlib.util

iFlag_cython = importlib.util.find_spec("cython")
if iFlag_cython is not None:
    from pyflowline.algorithms.cython.kernel import find_vertex_on_edge
else:
    from pyflowline.algorithms.auxiliary.find_index_in_list import find_vertex_on_edge


from pyflowline.classes.flowline import pyflowline
from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.edge import pyedge
from pyflowline.classes.confluence import pyconfluence
from pyflowline.formats.export_flowline import export_flowline_to_geojson
from pyflowline.formats.export_vertex import export_vertex_to_geojson


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

    def __init__(
        self, flowlines: List[pyflowline], pVertex_outlet: Optional[pyvertex] = None
    ):
        """
        Initialize the river network graph from flowlines.

        Args:
            flowlines: List of flowline objects representing the river network
            outlet_vertex: Optional outlet vertex for the drainage network. If provided,
                          enables outlet-based operations and optimizations.
        """
        self.aFlowline = flowlines
        self.pVertex_outlet_id: Optional[int] = None
        self.vertex_to_id: Dict[pyvertex, int] = {}
        self.id_to_vertex: Dict[int, pyvertex] = {}
        self.aFlowline_edges: Dict[int, Tuple[int, int]] = {}
        self.aVertex: List[pyvertex] = []  # Will be populated during graph building

        # Always use custom implementation
        self.adjacency_list: DefaultDict[int, List[Tuple[int, int]]] = defaultdict(list)
        self.in_degree: DefaultDict[int, int] = defaultdict(int)
        self.out_degree: DefaultDict[int, int] = defaultdict(int)
        self.aVertex_confluence: List[pyvertex] = []
        self.aConfluence: List[pyconfluence] = []

        # Initialize state management system

        logger.debug("Using custom graph implementation with state management")
        self._build_graph()
        # Set outlet vertex ID if outlet was provided
        if pVertex_outlet is not None:
            # only preserve flowlines that can flow out to the outlet
            self.pVertex_outlet = pVertex_outlet
            self.aFlowline = self.remove_disconnected_flowlines()
            self._set_outlet_vertex_id()

    def get_sources(self) -> List[int]:
        """Get source nodes (headwaters) with no incoming edges."""
        return [
            node_id
            for node_id in self.id_to_vertex.keys()
            if self.in_degree[node_id] == 0
        ]

    def get_sinks(self) -> List[int]:
        """Get sink nodes (outlets) with no outgoing edges."""
        return [
            node_id
            for node_id in self.id_to_vertex.keys()
            if self.out_degree[node_id] == 0
        ]

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

    def remove_disconnected_flowlines(
        self, pVertex_outlet: Optional[pyvertex] = None
    ) -> List[pyflowline]:
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
        flowlines = self.aFlowline

        # Use provided outlet vertex or fall back to stored one
        pVertex_target = (
            pVertex_outlet if pVertex_outlet is not None else self.pVertex_outlet
        )

        if pVertex_target is None:
            logger.warning(
                "No outlet vertex provided and none stored during initialization"
            )
            return flowlines

        # Use stored outlet vertex ID if it matches, otherwise find it
        outlet_vertex_id = None
        if pVertex_target == self.pVertex_outlet and self.pVertex_outlet_id is not None:
            outlet_vertex_id = self.pVertex_outlet_id
        else:
            # Find the outlet vertex ID in our graph
            for vertex_id, vertex in self.id_to_vertex.items():
                if vertex == pVertex_target:
                    outlet_vertex_id = vertex_id
                    break

        if outlet_vertex_id is None:
            logger.warning("Outlet vertex not found in network graph")
            return flowlines

        logger.info(
            f"Removing disconnected flowlines using outlet vertex {outlet_vertex_id}"
        )

        # Use a more robust approach: forward traversal to validate complete drainage paths
        # Step 1: Find all vertices that can reach the outlet (backward traversal)
        outlet_reachable_vertices = self._find_outlet_reachable_vertices(
            outlet_vertex_id
        )

        # Step 2: For each flowline, check if it contributes to the drainage network
        # A flowline contributes if:
        # 1. Its start vertex can reach the outlet, AND
        # 2. Its end vertex can reach the outlet, AND
        # 3. The flowline itself is part of a path that leads to the outlet
        reachable_flowlines = set()

        for flowline_idx, (start_id, end_id) in self.aFlowline_edges.items():
            # Check if both vertices are in the outlet-reachable set
            if (
                start_id in outlet_reachable_vertices
                and end_id in outlet_reachable_vertices
            ):
                # Additional check: ensure this flowline is on a path to outlet
                reachable_flowlines.add(flowline_idx)

        # Filter flowlines to keep only those that are reachable from outlet
        connected_flowlines = []
        disconnected_count = 0

        for i, flowline in enumerate(flowlines):
            if i in reachable_flowlines:
                connected_flowlines.append(flowline)
            else:
                disconnected_count += 1
                logger.debug(
                    f"Removing disconnected flowline {i}: {flowline.pVertex_start} -> {flowline.pVertex_end}"
                )

        logger.info(
            f"Removed {disconnected_count} disconnected flowlines, kept {len(connected_flowlines)} connected flowlines"
        )

        self._update_graph_flowlines(connected_flowlines)
        return connected_flowlines

    def remove_braided_river(self) -> List[pyflowline]:
        """
        Remove braided channels from the river network.
        Enhanced with state management and change tracking.

        This method identifies braided channels (multiple flowlines between
        the same vertex pair) and removes redundant ones to create a simplified
        river network without braiding. For each braided region, keeps the
        shortest flowline and removes the rest.

        Returns:
            List of flowlines with braided channels removed
        """

        braided_regions = self.find_braided_channels()
        if not braided_regions:
            logger.info("No braided channels found in the network")
            return self.aFlowline

        logger.info(f"Found {len(braided_regions)} braided regions to resolve")

        # Set to keep track of flowline indices to remove
        flowlines_to_remove = set()

        for region_flowlines in braided_regions:
            if len(region_flowlines) <= 1:
                continue

            # Find the shortest flowline in this region (keep it)
            shortest_flowline = min(
                region_flowlines, key=lambda fl: getattr(fl, "dLength", float("inf"))
            )

            # Mark all other flowlines in this region for removal
            for flowline in region_flowlines:
                if flowline != shortest_flowline:
                    # Find the index of this flowline in self.aFlowline
                    flowline_idx = self.aFlowline.index(flowline)
                    flowlines_to_remove.add(flowline_idx)
                    logger.debug(
                        f"Marking flowline {flowline_idx} for removal (braided channel, length={getattr(flowline, 'dLength', 'unknown')})"
                    )

        # Create new list of flowlines excluding those marked for removal
        simplified_flowlines = [
            flowline
            for idx, flowline in enumerate(self.aFlowline)
            if idx not in flowlines_to_remove
        ]

        logger.info(
            f"Removed {len(flowlines_to_remove)} braided flowlines, "
            f"resulting in {len(simplified_flowlines)} flowlines"
        )

        # Update the network using smart update system
        self._update_graph_flowlines(simplified_flowlines)
        self.merge_flowline()

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
        if not hasattr(self, "aFlowline") or not self.aFlowline:
            logger.warning(
                "No flowlines available in class instance for parallel river removal"
            )
            return []

        if len(self.aFlowline) <= 1:
            logger.debug("Skipping parallel path resolution: insufficient flowlines")
            return self.aFlowline.copy()

        logger.info(
            f"Removing parallel rivers from {len(self.aFlowline)} flowlines using graph-based approach"
        )

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
            paths = group["paths"]

            if len(paths) <= 1:
                continue

            # Calculate path significance for each complete path (multi-flowline groups)
            path_scores = []
            for path_idx, path_flowlines in enumerate(paths):
                # Calculate cumulative score for the entire path
                path_length = 0.0
                valid_flowlines = 0

                for flowline_idx in path_flowlines:
                    if flowline_idx >= len(self.aFlowline):
                        logger.warning(
                            f"Invalid flowline index {flowline_idx}, skipping"
                        )
                        continue

                    flowline = self.aFlowline[flowline_idx]
                    valid_flowlines += 1

                    # Only use length attributes if they have valid values (> 0)
                    length = getattr(flowline, "dLength", 0.0)
                    length = length if length > 0 else 0.0

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
            # loop through and remove safe flowlines
            for _, path_idx, path_flowlines in path_scores:
                # Reverse the flowlines in the path
                path_flowlines = list(path_flowlines)
                for flowline_idx in path_flowlines[
                    :-1
                ]:  # Exclude the last flowline in the path
                    # Check if the flowline has confluence points
                    flowline = self.aFlowline[flowline_idx]
                    # Check if the end point of the flowline is a confluence
                    end_vertex_id = self.vertex_to_id.get(flowline.pVertex_end)
                    if end_vertex_id is not None and self.in_degree[end_vertex_id] > 1:
                        flowlines_to_remove.add(flowline_idx)
                        logger.debug(
                            f"Removing flowline {flowline_idx} as its end point is a confluence"
                        )
                        break

        # Return flowlines with parallel paths removed (filter by index)
        try:
            result = [
                self.aFlowline[i]
                for i in range(len(self.aFlowline))
                if i not in flowlines_to_remove
            ]
            logger.info(
                f"Removed {len(flowlines_to_remove)} flowlines to resolve parallel paths"
            )
            # Validate result
            if len(result) == 0 and len(self.aFlowline) > 0:
                logger.warning(
                    "All flowlines were removed during parallel path resolution, returning original"
                )
                return self.aFlowline.copy()

            # Update the network using smart update system
            self._update_graph_flowlines(result)
            self.merge_flowline()

            return result
        except Exception as e:
            logger.error(
                f"Error filtering flowlines during parallel path resolution: {e}"
            )
            return self.aFlowline.copy()

    def remove_cycle(self) -> List[pyflowline]:
        """
        Remove cycles from the directed flow network using a single-pass
        **reverse-BFS from the outlet**.

        Algorithm
        ---------
        A valid flow network is a DAG where every flowline points toward the
        outlet.  We exploit this by growing a set ``known`` of vertices that
        are confirmed to have a path to the outlet, starting from the outlet
        itself and expanding upstream one flowline at a time.

        For each flowline ``start → end``:
        - If ``end ∈ known`` and ``start ∉ known``:
            The flowline correctly flows toward the outlet.
            Add ``start`` to ``known`` and enqueue it for further expansion.
        - If ``end ∈ known`` and ``start ∈ known``:
            Both endpoints are already confirmed outlet-reachable, but this
            flowline creates a redundant/cyclic connection — it flows from a
            vertex that is already upstream-reachable back toward a vertex
            that is also already known.  **Remove this flowline.**
        - If ``end ∉ known``:
            This flowline's end is not yet confirmed; skip for now (it will
            be processed when its end vertex is added to ``known``).

        This is O(V + E) — a single BFS pass — and requires no separate cycle
        detection step.  It naturally preserves all flowlines that contribute
        to the outlet-connected network and removes only those that create
        directed cycles.

        Fallback
        --------
        If no outlet vertex is available (``pVertex_outlet_id`` is None), cycle
        removal is skipped and the original flowlines are returned unchanged,
        because flow direction cannot be determined without an outlet reference.

        Returns
        -------
        List[pyflowline]
            Flowlines with all directed cycles removed.

        Example
        -------
        >>> river_graph = pyrivergraph(flowlines, outlet_vertex)
        >>> acyclic_flowlines = river_graph.remove_cycle()
        >>> print(f"Removed: {len(flowlines) - len(acyclic_flowlines)} flowlines")
        """
        if not hasattr(self, "aFlowline") or not self.aFlowline:
            logger.warning("No flowlines available in class instance for cycle removal")
            return []

        if len(self.aFlowline) <= 1:
            logger.debug("Skipping cycle removal: insufficient flowlines")
            return self.aFlowline.copy()

        logger.info(
            f"Removing cycles from {len(self.aFlowline)} flowlines using reverse-BFS approach"
        )

        outlet_id = getattr(self, "pVertex_outlet_id", None)

        if outlet_id is None:
            # Cycle removal requires an outlet vertex to determine flow direction.
            # Without it we cannot reliably identify which edge in a cycle flows
            # away from the outlet, so we skip removal and return unchanged.
            logger.warning(
                "remove_cycle: no outlet vertex available — "
                "cycle removal skipped (outlet is required to determine flow direction)"
            )
            return self.aFlowline.copy()

        # ----------------------------------------------------------------
        # Main path: reverse-BFS from the outlet.
        #
        # known_vertices: set of vertex IDs confirmed to have a path to the
        #                 outlet (either the outlet itself, or a start-vertex
        #                 of a flowline whose end is already in known_vertices).
        #
        # reverse_adj: maps end_vertex_id → list of (start_vertex_id, flowline_idx)
        #              so we can efficiently find all flowlines that flow INTO
        #              a given vertex (i.e. whose end == that vertex).
        # ----------------------------------------------------------------
        from collections import deque

        # Build reverse adjacency: end_id → [(start_id, flowline_idx), ...]
        reverse_adj: Dict[int, List[Tuple[int, int]]] = {}
        for fl_idx, (s_id, e_id) in self.aFlowline_edges.items():
            reverse_adj.setdefault(e_id, []).append((s_id, fl_idx))

        # known_vertices: vertex IDs that are confirmed outlet-reachable
        known_vertices: set = {outlet_id}
        flowlines_to_remove: set = set()

        queue = deque([outlet_id])

        while queue:
            current_end = queue.popleft()

            # Find all flowlines whose end vertex == current_end
            for start_id, fl_idx in reverse_adj.get(current_end, []):
                if fl_idx >= len(self.aFlowline):
                    continue  # guard against stale indices

                fl = self.aFlowline[fl_idx]

                if start_id not in known_vertices:
                    # Normal case: this flowline correctly flows toward outlet.
                    # Mark start_id as outlet-reachable and expand upstream.
                    known_vertices.add(start_id)
                    queue.append(start_id)
                    logger.debug(
                        f"  [remove_cycle] Accepted flowline idx={fl_idx} "
                        f"({start_id} → {current_end})"
                    )
                else:
                    # start_id is already in known_vertices — this flowline
                    # creates a cycle (flows from an already-known vertex back
                    # toward another already-known vertex).  Remove it.
                    flowlines_to_remove.add(fl_idx)
                    print(
                        f"  [remove_cycle] Cycle detected — removing flowline idx={fl_idx} "
                        f"start=({fl.pVertex_start.dLongitude_degree:.4f},"
                        f"{fl.pVertex_start.dLatitude_degree:.4f}) "
                        f"end=({fl.pVertex_end.dLongitude_degree:.4f},"
                        f"{fl.pVertex_end.dLatitude_degree:.4f}) "
                        f"order={fl.iStream_order} "
                        f"length={getattr(fl,'dLength',0):.1f}"
                    )

        print(
            f"  [remove_cycle] Reverse-BFS complete: "
            f"known={len(known_vertices)} vertices, "
            f"removing {len(flowlines_to_remove)} cyclic flowline(s)"
        )
        print(f"  [remove_cycle] Total flowlines before removal: {len(self.aFlowline)}")

        try:
            result = [
                self.aFlowline[i]
                for i in range(len(self.aFlowline))
                if i not in flowlines_to_remove
            ]
            print(f"  [remove_cycle] Total flowlines after removal: {len(result)}")

            if len(result) == 0 and len(self.aFlowline) > 0:
                print("  [remove_cycle] WARNING: all flowlines removed, returning original")
                return self.aFlowline.copy()

            self._update_graph_flowlines(result)
            return result
        except Exception as e:
            logger.error(f"Error filtering flowlines during cycle removal: {e}")
            return self.aFlowline.copy()

    def remove_small_river(
        self,
        dThreshold_small_river,
        nIterations=3,
        iFlag_debug=0,
        sWorkspace_output_basin=None,
    ):
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

        # Start with current flowlines
        aFlowline_current = self.aFlowline.copy()

        for i in range(nIterations):
            sStep = "{:02d}".format(i + 1)
            print(
                f"Iteration {sStep}: removing small rivers with threshold {dThreshold_small_river}"
            )

            # Step 1: Remove small rivers using graph-based approach
            aFlowline_filtered = self._remove_small_river_step(
                aFlowline_current, dThreshold_small_river
            )

            if iFlag_debug == 1 and sWorkspace_output_basin is not None:
                sFilename_out = f"flowline_large_{sStep}_before_intersect.geojson"
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_geojson(aFlowline_filtered, sFilename_out)

            # check if whether any flowlinw is actually removed
            if len(aFlowline_filtered) == len(aFlowline_current):
                print(f"No small rivers found in iteration {sStep}, stopping early.")
                break
            else:
                # Step 2: Update the graph with filtered flowlines
                self._update_graph_flowlines(aFlowline_filtered)

                # Step 3: Update stream order using graph-based method
                aFlowline_filtered = self.update_headwater_stream_order()

                # Step 4: Use graph-based merge instead of standalone functions
                # This replaces find_flowline_confluence + merge_flowline
                aFlowline_merged = self.merge_flowline()

                if iFlag_debug == 1 and sWorkspace_output_basin is not None:
                    sFilename_out = f"flowline_merge_{sStep}_before_intersect.geojson"
                    sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                    export_flowline_to_geojson(aFlowline_merged, sFilename_out)

                # Step 6: Final stream order update for this iteration
                aFlowline_current = self.update_headwater_stream_order()

                # Break if only one flowline remains
                if len(aFlowline_current) == 1:
                    break

        return aFlowline_current

    def remove_duplicate_flowlines(
        self, iFlag_direction_insensitive: bool = False
    ) -> List[pyflowline]:
        """
        Remove duplicate flowlines from the network.
        Enhanced with state management and change tracking.

        This method identifies and removes duplicate flowlines based on their start and end vertices
        to ensure each flowline in the network is unique.

        Args:
            iFlag_direction_insensitive (bool): If True, flowlines with opposite directions
                                              (A→B and B→A) are considered duplicates and BOTH are removed.
                                              If False (default), direction matters and they
                                              are treated as different flowlines.

        Returns:
            List of flowlines with duplicates removed

        Examples:
            # Default behavior - direction matters
            >>> unique_flowlines = river_graph.remove_duplicate_flowlines()

            # Direction-insensitive - A→B and B→A are considered duplicates and both are removed
            >>> unique_flowlines = river_graph.remove_duplicate_flowlines(iFlag_direction_insensitive=True)
        """

        if iFlag_direction_insensitive:
            # Special handling for direction-insensitive mode: remove both flowlines if they are opposite
            return self._remove_duplicate_flowlines_direction_insensitive()
        else:
            # Standard duplicate removal where direction matters
            return self._remove_duplicate_flowlines_direction_sensitive()

    def _remove_duplicate_flowlines_direction_sensitive(self) -> List[pyflowline]:
        """
        Remove duplicate flowlines where direction matters (A→B ≠ B→A).
        """
        seen_edges = set()
        unique_flowlines = []
        duplicates_count = 0

        for idx, flowline in enumerate(self.aFlowline):
            start_id = self.vertex_to_id.get(flowline.pVertex_start)
            end_id = self.vertex_to_id.get(flowline.pVertex_end)

            if start_id is None or end_id is None:
                logger.warning(
                    f"Flowline {idx} has unknown vertices, skipping duplicate check"
                )
                unique_flowlines.append(flowline)
                continue

            # Use original direction: (start_id, end_id) for direction-sensitive comparison
            edge_key = (start_id, end_id)

            if edge_key not in seen_edges:
                seen_edges.add(edge_key)
                unique_flowlines.append(flowline)
            else:
                duplicates_count += 1
                logger.debug(
                    f"Removing duplicate flowline {idx} (direction-sensitive): "
                    f"{flowline.pVertex_start} -> {flowline.pVertex_end}"
                )

        logger.info(
            f"Removed {duplicates_count} duplicate flowlines (direction-sensitive), "
            f"resulting in {len(unique_flowlines)} unique flowlines"
        )

        # Update the network using smart update system
        self._update_graph_flowlines(unique_flowlines)

        return unique_flowlines

    def _remove_duplicate_flowlines_direction_insensitive(self) -> List[pyflowline]:
        """
        Remove duplicate flowlines where direction doesn't matter.
        If two flowlines are in opposite directions (A→B and B→A), remove BOTH.
        """
        from collections import defaultdict

        # Group flowlines by normalized edge key (undirected)
        edge_groups = defaultdict(list)

        for idx, flowline in enumerate(self.aFlowline):
            start_id = self.vertex_to_id.get(flowline.pVertex_start)
            end_id = self.vertex_to_id.get(flowline.pVertex_end)

            if start_id is None or end_id is None:
                logger.warning(
                    f"Flowline {idx} has unknown vertices, skipping duplicate check"
                )
                continue

            # Normalize edge key: always use (smaller_id, larger_id) for direction-insensitive comparison
            edge_key = tuple(sorted([start_id, end_id]))
            edge_groups[edge_key].append((idx, flowline, (start_id, end_id)))

        unique_flowlines = []
        duplicates_count = 0
        opposite_pairs_count = 0

        for edge_key, flowline_group in edge_groups.items():
            if len(flowline_group) == 1:
                # Single flowline for this edge - keep it
                _, flowline, _ = flowline_group[0]
                unique_flowlines.append(flowline)
            else:
                # Multiple flowlines for this edge
                # Check if any are in opposite directions
                directions = set()
                for idx, flowline, (start_id, end_id) in flowline_group:
                    directions.add((start_id, end_id))

                # Check if we have opposite directions
                has_opposite_directions = False
                direction_list = list(directions)
                for i in range(len(direction_list)):
                    for j in range(i + 1, len(direction_list)):
                        dir1 = direction_list[i]
                        dir2 = direction_list[j]
                        # Check if dir2 is the reverse of dir1
                        if dir1 == (dir2[1], dir2[0]):
                            has_opposite_directions = True
                            break
                    if has_opposite_directions:
                        break

                if has_opposite_directions:
                    # Keep the flowline that has shorter path to outlet
                    if self.pVertex_outlet_id is not None:
                        flowline_path_lengths = []
                        for idx, flowline, (start_id, end_id) in flowline_group:
                            # Find path from end vertex to outlet
                            paths_to_outlet = self._find_all_paths(
                                end_id, self.pVertex_outlet_id, max_depth=1000
                            )
                            if paths_to_outlet:
                                # Use the shortest path length
                                shortest_path_length = min(
                                    len(path) for path in paths_to_outlet
                                )
                                flowline_path_lengths.append(
                                    (
                                        shortest_path_length,
                                        idx,
                                        flowline,
                                        (start_id, end_id),
                                    )
                                )
                            else:
                                # If no path to outlet, assign very large length
                                flowline_path_lengths.append(
                                    (float("inf"), idx, flowline, (start_id, end_id))
                                )

                        # Sort by path length (ascending - shorter path first)
                        flowline_path_lengths.sort(key=lambda x: x[0])

                        # Keep the one with shortest path, remove the rest
                        shortest_path_length, closest_idx, closest_flowline, _ = (
                            flowline_path_lengths[0]
                        )
                        unique_flowlines.append(closest_flowline)

                        # Remove the ones with longer paths
                        for (
                            path_length,
                            idx,
                            flowline,
                            (start_id, end_id),
                        ) in flowline_path_lengths[1:]:
                            opposite_pairs_count += 1
                            logger.debug(
                                f"Removing flowline {idx} (longer path to outlet, path_length={path_length}): "
                                f"{flowline.pVertex_start} -> {flowline.pVertex_end}"
                            )
                    else:
                        # No outlet defined, remove all as before
                        opposite_pairs_count += len(flowline_group)
                        for idx, flowline, (start_id, end_id) in flowline_group:
                            logger.debug(
                                f"Removing flowline {idx} (opposite direction pair, no outlet): "
                                f"{flowline.pVertex_start} -> {flowline.pVertex_end}"
                            )
                else:
                    # Same direction duplicates - keep only the first one
                    _, first_flowline, _ = flowline_group[0]
                    unique_flowlines.append(first_flowline)
                    duplicates_count += len(flowline_group) - 1
                    for idx, flowline, _ in flowline_group[1:]:
                        logger.debug(
                            f"Removing duplicate flowline {idx} (same direction): "
                            f"{flowline.pVertex_start} -> {flowline.pVertex_end}"
                        )

        # Add flowlines with unknown vertices (they were skipped in the grouping)
        for idx, flowline in enumerate(self.aFlowline):
            start_id = self.vertex_to_id.get(flowline.pVertex_start)
            end_id = self.vertex_to_id.get(flowline.pVertex_end)
            if start_id is None or end_id is None:
                unique_flowlines.append(flowline)

        total_removed = duplicates_count + opposite_pairs_count
        logger.info(
            f"Removed {total_removed} flowlines (direction-insensitive): "
            f"{duplicates_count} same-direction duplicates, "
            f"{opposite_pairs_count} opposite-direction pairs, "
            f"resulting in {len(unique_flowlines)} unique flowlines"
        )

        # Update the network using smart update system
        self._update_graph_flowlines(unique_flowlines)

        return unique_flowlines

    # ========================================================================
    # NETWORK ANALYSIS & DETECTION
    # ========================================================================

    def find_braided_channels(self) -> List[List[pyflowline]]:
        """
        Find braided channels (multiple flowlines between same vertex pair).
        Uses caching for performance optimization.

        Returns:
            List of braided regions, where each region is a list of braided flowlines
            Structure: [[region1_flowlines], [region2_flowlines], ...]
            Each region contains multiple flowlines connecting the same vertex pair
        """

        channel_groups: DefaultDict[Tuple[int, int], List[int]] = defaultdict(list)

        for flowline_idx, (start_id, end_id) in self.aFlowline_edges.items():
            channel_groups[(start_id, end_id)].append(flowline_idx)

        # Return only groups with multiple channels (braided) as separate regions
        braided_groups = []
        for (start_id, end_id), indices in channel_groups.items():
            if (
                len(indices) > 1
            ):  # Only braided channels (multiple flowlines between same vertices)
                region_flowlines = [self.aFlowline[idx] for idx in indices]
                braided_groups.append(region_flowlines)

        return braided_groups

    def find_parallel_paths(self) -> List[Dict]:
        """
        Find divergent parallel sections between vertices.

        This function identifies sections where paths diverge and reconverge, returning only
        the truly parallel (divergent) sections, not the shared portions.

        For example, given paths [14,15,13,11,0] and [14,15,27,13,11,0]:
        - Returns parallel section: [15,13] vs [15,27,13]
        - Does NOT return shared sections: [14,15] or [13,11,0]

        Returns:
            List of parallel path group dictionaries, each containing:
            - 'paths': List of paths, where each path is a list of flowline indices
            - 'start_vertex': Start vertex ID (divergence point)
            - 'end_vertex': End vertex ID (reconvergence point)
        """
        # Only find partial parallel sections (divergent sections)
        return self._find_partial_parallel_sections()

    def _find_partial_parallel_sections(self) -> List[Dict]:
        """
        Find partial parallel sections within longer paths.

        This identifies divergence and reconvergence points in the network where
        multiple paths split and rejoin, returning only the divergent sections.

        For example, given paths [14,15,13,11,0] and [14,15,27,13,11,0]:
        - Identifies divergence at vertex 15 and reconvergence at vertex 13
        - Returns parallel sections: [15,13] vs [15,27,13]
        """
        parallel_sections = []

        # Find all vertices that have multiple outgoing paths (divergence points)
        divergence_vertices = []
        for vertex_id in self.adjacency_list.keys():
            if self.out_degree[vertex_id] > 1:
                divergence_vertices.append(vertex_id)

        logger.debug(f"Found {len(divergence_vertices)} potential divergence vertices")

        # For each divergence vertex, find where the paths reconverge
        for divergence_id in divergence_vertices:
            reconvergence_id = self._find_reconvergence_points(divergence_id)

            if reconvergence_id is not None:
                # Find all paths between divergence and reconvergence points
                divergent_paths = self._find_all_paths(
                    divergence_id, reconvergence_id, max_depth=8
                )

                if len(divergent_paths) > 1:
                    # Convert paths to flowline indices
                    path_flowline_groups = []
                    for path in divergent_paths:
                        path_flowlines = self._path_to_flowlines(path)
                        if path_flowlines and len(path_flowlines) > 0:
                            path_flowline_groups.append(path_flowlines)

                    # Ensure we have multiple unique divergent paths
                    if len(path_flowline_groups) > 1:
                        # Verify paths are truly different
                        unique_paths = []
                        for path_flowlines in path_flowline_groups:
                            path_signature = tuple(sorted(path_flowlines))
                            if path_signature not in [
                                tuple(sorted(up)) for up in unique_paths
                            ]:
                                unique_paths.append(path_flowlines)

                        # Only add if we have multiple unique divergent paths
                        if len(unique_paths) > 1:
                            logger.debug(
                                f"Found {len(unique_paths)} partial parallel paths between vertices {divergence_id} -> {reconvergence_id}"
                            )
                            parallel_sections.append(
                                {
                                    "paths": unique_paths,
                                    "start_vertex": divergence_id,
                                    "end_vertex": reconvergence_id,
                                }
                            )

        return parallel_sections

    def _find_reconvergence_points(
        self, divergence_id: int, max_depth: int = 8
    ) -> List[int]:
        """
        Find vertices where paths from a divergence point reconverge.

        Args:
            divergence_id: The vertex ID where paths diverge
            max_depth: Maximum search depth to prevent infinite loops

        Returns:
            List of vertex IDs where paths reconverge
        """

        # Get all immediate downstream vertices from the divergence point
        downstream_vertices = [
            neighbor_id for neighbor_id, _ in self.adjacency_list[divergence_id]
        ]

        if len(downstream_vertices) < 2:
            return None  # Need at least 2 paths to have reconvergence

        # For each pair of downstream paths, find where they reconverge
        for i in range(len(downstream_vertices)):
            for j in range(i + 1, len(downstream_vertices)):
                path1_start = downstream_vertices[i]
                path2_start = downstream_vertices[j]

                # Find common downstream vertices reachable from both paths
                path1_reachable = self._get_reachable_vertices(path1_start, max_depth)
                path2_reachable = self._get_reachable_vertices(path2_start, max_depth)

                # find the first identical vertex in both reachable lists
                common_vertex = self.find_first_common_vertex(
                    path1_reachable, path2_reachable
                )

        return common_vertex

    def find_first_common_vertex(self, reachable_list1, reachable_list2):
        """
        Find the first vertex that appears in both lists, respecting topological order.
        Returns the vertex that appears earliest in the flow sequence.
        """
        # Convert second list to set for O(1) lookup
        set2 = set(reachable_list2)

        # Find first vertex in list1 that also exists in list2
        for vertex in reachable_list1:
            if vertex in set2:
                return vertex

        return None  # No common vertex found

    def _get_reachable_vertices(self, start_id: int, max_depth: int) -> Set[int]:
        """
        Get all vertices reachable from a starting vertex within max_depth.

        Args:
            start_id: Starting vertex ID
            max_depth: Maximum search depth

        Returns:
            Set of reachable vertex IDs
        """
        reachable = list()
        queue = deque([(start_id, 0)])  # (vertex_id, depth)
        visited = set([start_id])

        while queue:
            current_id, depth = queue.popleft()

            if depth >= max_depth:
                continue

            reachable.append(current_id)
            # Add downstream vertices
            for neighbor_id, _ in self.adjacency_list[current_id]:
                if neighbor_id not in visited:
                    visited.add(neighbor_id)
                    queue.append((neighbor_id, depth + 1))

        return reachable

    def _find_linear_segments(self) -> List[List[int]]:
        """
        Find linear segments that can be merged into single flowlines.
        Uses caching for performance optimization.

        A linear segment is a chain of flowlines where each intermediate vertex
        has exactly one incoming and one outgoing edge (degree 2 vertices).
        These represent artificially segmented flowlines that should be merged.

        Fixed to handle bidirectional segment detection - extends both upstream
        and downstream to capture complete linear segments, including those that
        start from confluences.

        Returns:
            List of linear segments, where each segment is a list of flowline indices
            that should be merged together in order from upstream to downstream.
        """

        linear_segments = []
        processed_flowlines = set()

        # Find all vertices with degree 2 (one in, one out) - these are merge candidates
        merge_vertices = set()
        for vertex_id in self.id_to_vertex.keys():
            if self.in_degree[vertex_id] == 1 and self.out_degree[vertex_id] == 1:
                merge_vertices.add(vertex_id)

        logger.debug(
            f"Found {len(merge_vertices)} degree-2 vertices for potential merging"
        )

        # For each flowline, try to build a complete linear segment using bidirectional search
        for flowline_idx, (start_id, end_id) in self.aFlowline_edges.items():
            if flowline_idx in processed_flowlines:
                continue

            # Build segment by extending both upstream and downstream
            segment = self._build_bidirectional_segment(
                flowline_idx, start_id, end_id, merge_vertices, processed_flowlines
            )

            # Only consider segments with multiple flowlines
            if len(segment) > 1:
                linear_segments.append(segment)
                processed_flowlines.update(segment)
                logger.debug(
                    f"Found linear segment with {len(segment)} flowlines: {segment}"
                )

        logger.info(
            f"Found {len(linear_segments)} linear segments for potential merging"
        )

        return linear_segments

    def _build_bidirectional_segment(
        self,
        start_flowline_idx: int,
        start_vertex_id: int,
        end_vertex_id: int,
        merge_vertices: set,
        processed_flowlines: set,
    ) -> List[int]:
        """
        Build a complete linear segment by extending both upstream and downstream
        from the starting flowline.

        Args:
            start_flowline_idx: Index of the flowline to start building from
            start_vertex_id: Start vertex ID of the starting flowline
            end_vertex_id: End vertex ID of the starting flowline
            merge_vertices: Set of degree-2 vertices that can be merged
            processed_flowlines: Set of already processed flowline indices

        Returns:
            List of flowline indices in upstream-to-downstream order
        """
        # Start with the current flowline
        upstream_segment = []
        downstream_segment = [start_flowline_idx]

        # Extend upstream from start vertex
        current_start = start_vertex_id
        while current_start in merge_vertices:
            # Find upstream flowlines that end at this vertex
            upstream_flowlines = []
            for fl_idx, (fl_start, fl_end) in self.aFlowline_edges.items():
                if fl_end == current_start and fl_idx not in processed_flowlines:
                    upstream_flowlines.append(fl_idx)

            if len(upstream_flowlines) == 1:
                upstream_flowline_idx = upstream_flowlines[0]
                upstream_segment.insert(0, upstream_flowline_idx)  # Insert at beginning
                # Update current_start to the start of the upstream flowline
                current_start, _ = self.aFlowline_edges[upstream_flowline_idx]
            else:
                break

        # Extend downstream from end vertex
        current_end = end_vertex_id
        while current_end in merge_vertices:
            # Find downstream flowlines that start at this vertex
            downstream_flowlines = [
                fl_idx for neighbor_id, fl_idx in self.adjacency_list[current_end]
            ]

            if len(downstream_flowlines) == 1:
                downstream_flowline_idx = downstream_flowlines[0]
                if downstream_flowline_idx not in processed_flowlines:
                    downstream_segment.append(downstream_flowline_idx)
                    # Update current_end to the end of the downstream flowline
                    if downstream_flowline_idx in self.aFlowline_edges:
                        _, current_end = self.aFlowline_edges[downstream_flowline_idx]
                    else:
                        break
                else:
                    break
            else:
                break

        # Combine upstream and downstream segments
        complete_segment = upstream_segment + downstream_segment

        return complete_segment

    # ========================================================================
    # NETWORK MODIFICATION & PROCESSING
    # ========================================================================

    def split_flowline(
        self,
        aVertex_in: Optional[List[pyvertex]] = None,
        iFlag_intersect=None,
        iFlag_use_id=None,
    ) -> List[pyflowline]:
        """
        Split flowline based on the intersection with a list of vertex
        Input:

            aVertex_in: list of vertex (optional)
            iFlag_intersect: 1: a vertex maybe on a line, but it is not a vertex of the line
            iFlag_use_id: 1:
        Output:
            aFlowline_out: list of flowline
        """
        iFlag_graph_update = 0
        aFlowline_in = self.aFlowline
        aFlowline_out = list()
        nFlowline = len(aFlowline_in)

        if aVertex_in is None:
            aVertex_in = self.aVertex

        aVertex_in_set = set(aVertex_in)

        for i in range(nFlowline):
            pFlowline = aFlowline_in[i]
            iStream_order = pFlowline.iStream_order
            # iStream_segment = pFlowline.iStream_segment
            iFlag_dam = pFlowline.iFlag_dam
            # nVertex = pFlowline.nVertex
            nEdge = pFlowline.nEdge
            iPart = 0
            aVertex = list()  # the actual vertex of ROI
            aVertex_all = (
                list()
            )  # include vertex that is not ROI, but we need them to subset
            for j in range(nEdge):
                pEdge = pFlowline.aEdge[j]
                pVertex = pEdge.pVertex_start
                aVertex_all.append(pVertex)

                if aVertex_in and pVertex in aVertex_in_set:
                    iPart += 1
                    aVertex.append(pVertex)

                if iFlag_intersect is not None:
                    if iFlag_use_id is not None:
                        aDistance = list()
                        iFlag_exist = 0
                        aVertex_dummy = list()
                        aIndex = list()
                        npoint = 0
                        for k in range(len(aVertex_in)) if aVertex_in else []:
                            pVertex0 = aVertex_in[k]
                            if pVertex0.lFlowlineID == pFlowline.lFlowlineID:
                                iFlag_exist = 1
                                distance = pEdge.pVertex_start.calculate_distance(
                                    pVertex0
                                )
                                if distance == 0:
                                    continue
                                iPart = iPart + 1
                                aDistance.append(distance)
                                aVertex_dummy.append(pVertex0)
                                aIndex.append(k)
                                npoint = npoint + 1
                            else:
                                pass

                        if aDistance:
                            aIndex_order = np.argsort(aDistance)
                            for k in aIndex_order:
                                pVertex_dummy = aVertex_in[aIndex[k]]
                                aVertex.append(pVertex_dummy)
                                aVertex_all.append(pVertex_dummy)

                    else:
                        iFlag_exist, npoint, aIndex = (
                            find_vertex_on_edge(aVertex_in, pEdge)
                            if aVertex_in
                            else (0, 0, [])
                        )
                        # they must be ordered
                        if iFlag_exist == 1:
                            for m in range(npoint):
                                pVertex_dummy = aVertex_in[aIndex[m]]
                                iPart = iPart + 1
                                aVertex.append(pVertex_dummy)
                                aVertex_all.append(pVertex_dummy)

            # the last ending vertex
            pVertex = pFlowline.pVertex_end
            aVertex_all.append(pVertex)

            if aVertex_in and pVertex in aVertex_in_set:
                iPart = iPart + 1
                aVertex.append(pVertex)
            if iPart == 0:
                print("Something is wrong")
            else:
                if iPart == 1:
                    if iFlag_use_id is not None:
                        pass
                    pass
                else:
                    if iPart >= 2:
                        nLine = iPart - 1
                        aVertex_index = [
                            aVertex_all.index(pVertex)
                            for pVertex in aVertex
                            if pVertex in aVertex_all
                        ]

                        # find duplicate
                        for k in range(nLine):
                            t = aVertex_index[k]
                            s = aVertex_index[k + 1]
                            if s != t:
                                aEdge = [
                                    pyedge.create(aVertex_all[l], aVertex_all[l + 1])
                                    for l in range(t, s)
                                ]
                                # Remove None from the aEdge list
                                aEdge = [edge for edge in aEdge if edge is not None]
                                pFlowline1 = pyflowline(aEdge)
                                pFlowline1.iStream_order = iStream_order
                                pFlowline1.iFlag_dam = iFlag_dam
                                aFlowline_out.append(pFlowline1)
                                nEdge_new = pFlowline1.nEdge
                                if nEdge_new != nEdge:
                                    iFlag_graph_update = 1
                            pass
                        pass
                    pass

        # should we update the graph here?
        if iFlag_graph_update == 1:
            self._update_graph_flowlines(aFlowline_out)
        return aFlowline_out

    def merge_flowline(self) -> List[pyflowline]:
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

        flowlines = self.aFlowline

        try:
            linear_segments = self._find_linear_segments()

            if not linear_segments:
                logger.info("No linear segments found for merging")
                return self.aFlowline

            logger.info(f"Found {len(linear_segments)} linear segments to merge")

            # Create new flowlines by merging segments
            merged_flowlines = []
            processed_indices = set()

            for segment in linear_segments:
                if len(segment) < 2:
                    continue

                try:
                    # Get flowlines in the segment
                    segment_flowlines = [
                        flowlines[idx] for idx in segment if idx < len(flowlines)
                    ]

                    if len(segment_flowlines) < 2:
                        continue

                    # Create merged flowline from the segment
                    merged_flowline = self._merge_flowline_segment(segment_flowlines)
                    if merged_flowline is not None:
                        merged_flowlines.append(merged_flowline)
                        processed_indices.update(segment)
                        logger.debug(
                            f"Merged {len(segment_flowlines)} flowlines into single flowline"
                        )

                except Exception as e:
                    logger.error(f"Error merging segment {segment}: {e}")
                    continue

            # Add non-merged flowlines
            for i, flowline in enumerate(flowlines):
                if i not in processed_indices:
                    merged_flowlines.append(flowline)

            logger.info(
                f"Merged {len(linear_segments)} segments, "
                f"resulting in {len(merged_flowlines)} flowlines "
                f"(reduced from {len(flowlines)})"
            )

            self._update_graph_flowlines(merged_flowlines)
            return merged_flowlines

        except Exception as e:
            logger.error(f"Error during linear segment merging: {e}")
            return flowlines

    # ========================================================================
    # STREAM ANALYSIS & TOPOLOGY
    # ========================================================================

    def update_headwater_stream_order(self) -> List[pyflowline]:
        """
        Update stream order for head water flowlines using graph-based approach.
        Enhanced with state management and change tracking.

        This method replaces the standalone function by leveraging the graph structure
        to identify headwater flowlines and update their stream order based on
        network topology and confluence analysis.

        Returns:
            List[pyflowline]: Updated flowlines with corrected stream orders
        """
        if not hasattr(self, "aFlowline") or not self.aFlowline:
            logger.warning(
                "No flowlines available in class instance for stream order update"
            )
            return []

        logger.info(
            f"Updating stream order for {len(self.aFlowline)} flowlines using graph-based approach"
        )

        try:
            # Step 1: Identify headwater flowlines (sources in the graph)
            aFlowline_headwater = self.identify_headwater_flowlines()
            if not aFlowline_headwater:
                logger.warning(
                    "No headwater flowlines identified for stream order update"
                )
                return self.aFlowline.copy()
            else:
                # set the stream order of headwater to 1
                for pFlowline in aFlowline_headwater:
                    pFlowline.iStream_order = 1

            return self.aFlowline

        except Exception as e:
            logger.error(f"Error updating stream orders: {e}")
            return self.aFlowline.copy()

    def identify_headwater_flowlines(self) -> List[pyflowline]:
        """
        Identify headwater flowlines (sources) in the network.

        Returns:
            List of flowlines that are headwaters (no upstream connections)
        """
        aFlowline_headwater = []
        for i, flowline in enumerate(self.aFlowline):
            start_vertex_id = self.vertex_to_id.get(flowline.pVertex_start)
            if start_vertex_id is not None and self.in_degree[start_vertex_id] == 0:
                aFlowline_headwater.append(flowline)

        logger.debug(f"Identified {len(aFlowline_headwater)} headwater flowlines")
        self.aFlowline_headwater = aFlowline_headwater
        return self.aFlowline_headwater

    def define_stream_segment(self):
        """
        Define stream segments using topological sorting from outlet to headwater.

        This method ensures flowlines are properly ordered from outlet upstream,
        regardless of their input order, using graph topology. Maintains compatibility
        with original numbering: outlet gets largest segment number, decreasing upstream.

        Returns:
            Tuple of (sorted_flowlines, segment_indices)
        """
        from collections import deque

        nFlowline = len(self.aFlowline)

        if nFlowline == 0:
            print("data incomplete")
            return [], []

        try:
            # Step 1: Sort flowlines using BFS from outlet
            sorted_flowlines = self._sort_flowlines_from_outlet()

            # Step 2: Assign segment numbers (outlet gets largest number, decreasing upstream)
            # This maintains compatibility with original method: nFlowline - i + 1
            nFlowline = len(sorted_flowlines)
            for i, flowline in enumerate(sorted_flowlines, start=1):
                flowline.iStream_segment = nFlowline - i + 1

            # Step 3: Update internal flowline list to maintain sorted order
            self.aFlowline = sorted_flowlines
            self._update_graph_flowlines(sorted_flowlines)

            # Step 4: Create segment index list
            aStream_segment = [
                flowline.iStream_segment for flowline in sorted_flowlines
            ]

            logger.info(
                f"Topologically sorted {len(sorted_flowlines)} flowlines from outlet to headwater"
            )

            return sorted_flowlines, aStream_segment

        except Exception as e:
            logger.error(f"Error in topological sorting: {e}, using simple enumeration")
            # Fallback to simple enumeration on error
            for i, pFlowline in enumerate(self.aFlowline, start=1):
                pFlowline.iStream_segment = nFlowline - i + 1
            aStream_segment = [
                pFlowline.iStream_segment for pFlowline in self.aFlowline
            ]
            return self.aFlowline, aStream_segment

    def _sort_flowlines_from_outlet(self):
        """
        Core sorting function using BFS traversal from outlet vertex.

        Returns:
            List of flowlines sorted from outlet to headwater
        """
        from collections import deque

        visited_flowlines = set()
        sorted_flowlines = []
        queue = deque()

        # Start from outlet vertex - ensure it's valid after graph modifications
        outlet_id = self.pVertex_outlet_id

        # Safety check: if outlet_id is invalid, try to refresh it
        if outlet_id is None or outlet_id not in self.id_to_vertex:
            logger.warning("Outlet vertex ID is invalid, attempting to refresh")
            if hasattr(self, "pVertex_outlet") and self.pVertex_outlet is not None:
                self._set_outlet_vertex_id()
                outlet_id = self.pVertex_outlet_id

            if outlet_id is None or outlet_id not in self.id_to_vertex:
                logger.error(
                    "Cannot find valid outlet vertex after refresh, returning original flowline order"
                )
                return self.aFlowline

        # Find all flowlines that end at the outlet (these are the most downstream)
        outlet_flowlines = []
        for flowline_idx, (start_id, end_id) in self.aFlowline_edges.items():
            if end_id == outlet_id:
                outlet_flowlines.append(flowline_idx)

        if not outlet_flowlines:
            logger.warning("No flowlines found ending at outlet vertex")
            return self.aFlowline  # Return original order if no outlet connections

        # Add outlet flowlines to queue
        for flowline_idx in outlet_flowlines:
            queue.append(flowline_idx)

        # BFS traversal upstream from outlet
        while queue:
            flowline_idx = queue.popleft()

            if flowline_idx in visited_flowlines:
                continue

            visited_flowlines.add(flowline_idx)

            # Add flowline to sorted list
            if flowline_idx < len(self.aFlowline):
                sorted_flowlines.append(self.aFlowline[flowline_idx])

            # Find upstream flowlines (those that end where this one starts)
            start_id, _ = self.aFlowline_edges[flowline_idx]

            upstream_flowlines = []
            for upstream_idx, (up_start, up_end) in self.aFlowline_edges.items():
                if up_end == start_id and upstream_idx not in visited_flowlines:
                    upstream_flowlines.append(upstream_idx)

            # Add upstream flowlines to queue (they'll be processed after current level)
            for upstream_idx in upstream_flowlines:
                if upstream_idx not in queue:
                    queue.append(upstream_idx)

        # Log any disconnected flowlines that were excluded
        excluded_count = len(self.aFlowline) - len(sorted_flowlines)
        if excluded_count > 0:
            logger.info(
                f"Excluded {excluded_count} disconnected flowlines that do not flow to outlet"
            )

            # Log details of excluded flowlines for debugging
            for i, flowline in enumerate(self.aFlowline):
                if i not in visited_flowlines:
                    logger.debug(
                        f"Excluded disconnected flowline {i} (ID: {getattr(flowline, 'lFlowlineID', 'unknown')})"
                    )

        logger.info(
            f"Sorted {len(sorted_flowlines)} connected flowlines using BFS from outlet"
        )
        return sorted_flowlines

    def define_river_confluence(self):
        """
        Build the conflence using the in_degree

        Args:
            aFlowline_basin_in (list [pyflowline]): A list of flowlines in this basin
            aVertex_confluence_in (list [pyconfluence]): A list of vertices in this basin

        Returns:
            list [pyconfluence]: A list of confluences in this basin
        """
        # this can only be calculated for confluence
        # Create a dictionary to map each vertex to its upstream and downstream flowlines
        aConfluence = []
        try:
            for vertex_id, vertex in self.id_to_vertex.items():
                if self.in_degree[vertex_id] > 1:  # Confluence point
                    pConfluence = self._create_confluence_object(vertex_id, vertex)
                    if pConfluence:
                        aConfluence.append(pConfluence)
        except:
            logger.error("Error defining river confluences")
            return []

        self.aConfluence = aConfluence

        return aConfluence

    def define_stream_topology(self) -> List[pyflowline]:
        """
        Define comprehensive stream topology including confluences and stream segments.
        Enhanced with state management and caching.

        Returns:
            Dictionary containing topology information:
            - confluences: List of confluence objects
            - stream_segments: List of stream segment definitions
            - topology_stats: Statistics about the network topology
        """
        aFlowline = self.aFlowline
        try:
            # use the connectivity information to set up the topology
            for flowline in aFlowline:
                # get upstream and downstream flowlines
                upstream_flowlines = [
                    self.aFlowline[i] for i in self.get_upstream_indices(flowline)
                ]
                downstream_flowlines = [
                    self.aFlowline[i] for i in self.get_downstream_indices(flowline)
                ]
                flowline.aFlowline_upstream = upstream_flowlines
                flowline.aFlowline_downstream = downstream_flowlines
            pass

        except Exception as e:
            logger.error(f"Error defining stream topology: {e}")

        return aFlowline

    def define_stream_order(self, iFlag_so_method_in: int = 1) -> List[pyflowline]:
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
        if not hasattr(self, "aFlowline") or not self.aFlowline:
            logger.warning("No flowlines available for stream order calculation")
            return None

        logger.info(
            f"Defining stream order for {len(self.aFlowline)} flowlines using "
            f"{'Strahler' if iFlag_so_method_in == 1 else 'Shreve'} method"
        )

        try:
            # Step 1: Initialize headwater flowlines (stream order = 1)
            # which should be already set in previous steps
            self.identify_headwater_flowlines()
            self.update_headwater_stream_order()
            # Step 4: Process confluences iteratively until all stream orders are defined
            self._process_confluences_iteratively(iFlag_so_method_in)

            aStream_order = [flowline.iStream_order for flowline in self.aFlowline]

            return self.aFlowline, aStream_order

        except Exception as e:
            logger.error(f"Error defining stream order: {e}")
            return None, None

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

        logger.debug(
            f"Built graph with {len(self.id_to_vertex)} vertices and {len(self.aFlowline_edges)} edges"
        )

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

        start_time = time.time()
        self.aFlowline = new_flowlines
        self._build_graph()

        # Critical: Update outlet vertex ID after graph rebuild
        # The outlet vertex ID can become invalid after graph modifications
        if hasattr(self, "pVertex_outlet") and self.pVertex_outlet is not None:
            self._set_outlet_vertex_id()
            logger.debug(
                f"Updated outlet vertex ID after graph rebuild: {self.pVertex_outlet_id}"
            )

        return

    # ========================================================================
    # PRIVATE METHODS - PATH FINDING & ANALYSIS
    # ========================================================================

    def _find_all_paths(
        self, start_id: int, target_id: int, max_depth: int = 10
    ) -> List[List[int]]:
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

        def dfs_paths(
            current_id: int,
            target_id: int,
            path: List[int],
            visited: Set[int],
            depth: int,
        ):
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

    def _remove_small_river_step(
        self, flowlines: List[pyflowline], threshold: float
    ) -> List[pyflowline]:
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
            length = getattr(flowline, "dLength", 0.0)
            if length >= threshold:
                filtered_flowlines.append(flowline)
            else:
                logger.debug(f"Removing small river with length {length}")

        logger.info(f"Removed {len(flowlines) - len(filtered_flowlines)} small rivers")
        return filtered_flowlines

    def _merge_flowline_segment(
        self, segment_flowlines: List[pyflowline]
    ) -> Optional[pyflowline]:
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
            # assume that they are connected in order
            # merge using the built-in api in the reverse order
            merged_flowline = segment_flowlines[-1]
            for flowline in segment_flowlines[-2::-1]:
                flowline_up = flowline
                merged_flowline = merged_flowline.merge_upstream(flowline_up)

            return merged_flowline

        except Exception as e:
            logger.error(f"Error merging flowline segment: {e}")
            return None

    # ========================================================================
    # PRIVATE METHODS - VERTEX & CONFLUENCE MANAGEMENT
    # ========================================================================

    def _create_confluence_object(
        self, vertex_id: int, vertex: pyvertex
    ) -> Optional[pyconfluence]:
        """
        Create a confluence object for a vertex.
        Using the right constructor pConfluence = pyconfluence(pVertex, aFlowline_upstream, aFlowline_downstream)

        Args:
            vertex_id: Vertex ID
            vertex: Vertex object

        Returns:
            Optional[pyconfluence]: Confluence object or None if creation failed
        """
        try:
            # Collect upstream flowlines (incoming to this vertex)
            aFlowline_upstream = []
            for start_id, neighbors in self.adjacency_list.items():
                for end_id, flowline_idx in neighbors:
                    if end_id == vertex_id and flowline_idx < len(self.aFlowline):
                        aFlowline_upstream.append(self.aFlowline[flowline_idx])

            # Collect downstream flowlines (outgoing from this vertex)
            aFlowline_downstream = []
            for end_id, flowline_idx in self.adjacency_list[vertex_id]:
                if flowline_idx < len(self.aFlowline):
                    aFlowline_downstream.append(self.aFlowline[flowline_idx])

            # Create confluence using the proper constructor
            confluence = pyconfluence(vertex, aFlowline_upstream, aFlowline_downstream)

            # Set additional attributes
            confluence.lVertexID = vertex_id
            # Note: nUpstream and nDownstream are already set by the constructor

            return confluence

        except Exception as e:
            logger.error(f"Error creating confluence object: {e}")
            return None

    # ========================================================================
    # PRIVATE METHODS - STREAM ORDER CALCULATION HELPERS
    # ========================================================================

    def _process_confluences_iteratively(self, iFlag_so_method_in: int):
        """
        Process confluences iteratively to update stream orders.

        Args:
            confluence_index: Spatial index for confluences
            iFlag_so_method_in: Stream ordering method (1=Strahler, 2=Shreve)
        """
        # Simple iterative approach - can be optimized
        unresolved_confluences = set(range(len(self.aConfluence)))
        iteration = 0

        while unresolved_confluences:
            iteration += 1
            logger.debug(
                f"Stream order iteration {iteration}, "
                f"{len(unresolved_confluences)} confluences remaining"
            )

            resolved_this_iteration = set()

            for idx in unresolved_confluences:
                confluence = self.aConfluence[idx]
                upstream_orders = [
                    fl.iStream_order
                    for fl in confluence.aFlowline_upstream
                    if fl.iStream_order > 0
                ]

                if len(upstream_orders) == len(confluence.aFlowline_upstream):
                    # All upstream orders are known, calculate downstream order
                    if iFlag_so_method_in == 1:  # Strahler
                        if upstream_orders.count(max(upstream_orders)) > 1:
                            downstream_order = max(upstream_orders) + 1
                        else:
                            downstream_order = max(upstream_orders)
                    elif iFlag_so_method_in == 2:  # Shreve
                        downstream_order = sum(upstream_orders)
                    else:
                        logger.warning(
                            f"Unknown stream order method: {iFlag_so_method_in}"
                        )
                        continue

                    # Update downstream flowlines
                    for fl in confluence.aFlowline_downstream:
                        fl.iStream_order = downstream_order

                    resolved_this_iteration.add(idx)

            if not resolved_this_iteration:
                logger.warning(
                    "No confluences resolved in this iteration, stopping to avoid infinite loop"
                )
                break

            unresolved_confluences -= resolved_this_iteration

        return

    def get_upstream_indices(self, flowline: pyflowline) -> List[int]:
        """
        Get indices of upstream flowlines for a given flowline.

        Args:
            flowline: Flowline object

        Returns:
            List of indices of upstream flowlines
        """
        upstream_indices = []
        # Get the start vertex ID of the given flowline
        start_vertex_id = self.vertex_to_id.get(flowline.pVertex_start)
        if start_vertex_id is not None:
            # Find flowlines that end at this flowline's start vertex
            for start_id, neighbors in self.adjacency_list.items():
                for end_id, flowline_idx in neighbors:
                    if end_id == start_vertex_id and flowline_idx < len(self.aFlowline):
                        upstream_indices.append(flowline_idx)
        return upstream_indices

    def get_downstream_indices(self, flowline: pyflowline) -> List[int]:
        """
        Get indices of downstream flowlines for a given flowline.

        Args:
            flowline: Flowline object

        Returns:
            List of indices of downstream flowlines
        """
        downstream_indices = []
        # Get the end vertex ID of the given flowline
        end_vertex_id = self.vertex_to_id.get(flowline.pVertex_end)
        if end_vertex_id is not None:
            # Find flowlines that start at this flowline's end vertex
            for neighbor_id, flowline_idx in self.adjacency_list[end_vertex_id]:
                if flowline_idx < len(self.aFlowline):
                    downstream_indices.append(flowline_idx)
        return downstream_indices
