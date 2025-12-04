
"""
Incremental updates module for pyrivergraph.

This module provides efficient incremental update methods for simple operations
that don't require full graph rebuilds, including single flowline additions/removals,
vertex updates, and localized topology changes.
"""

import time
import logging
from typing import List, Dict, Any, Optional, Set, Tuple, Union
from collections import defaultdict

from .state_management import GraphChangeType, GraphStateManager
from .utils import calculate_distance_2d, find_cycles_dfs

logger = logging.getLogger(__name__)


class IncrementalUpdater:
    """
    Handles incremental updates to the river network graph for simple operations
    that can be performed without full rebuilds.
    """

    def __init__(self, graph_instance: Any, state_manager: Optional[GraphStateManager] = None):
        """
        Initialize the incremental updater.

        Args:
            graph_instance: The pyrivergraph instance to update
            state_manager: Optional state manager for tracking changes
        """
        self.graph = graph_instance
        self.state_manager = state_manager
        self.tolerance = 1e-6  # Distance tolerance for geometric operations

        logger.debug("Initialized IncrementalUpdater")

    def add_single_flowline(self, flowline: Any,
                           update_topology: bool = True,
                           update_stream_order: bool = False) -> bool:
        """
        Add a single flowline to the graph with incremental updates.

        Args:
            flowline: The flowline object to add
            update_topology: Whether to update topology after addition
            update_stream_order: Whether to recalculate stream order

        Returns:
            True if addition was successful
        """
        start_time = time.time()

        try:
            # Get current flowline count for ID assignment
            flowline_id = len(self.graph.aFlowline)
            flowline.lFlowlineIndex = flowline_id

            # Add to flowline list
            self.graph.aFlowline.append(flowline)

            # Update vertices incrementally
            affected_vertices = self._update_vertices_for_flowline(flowline, add=True)

            # Update confluences incrementally
            affected_confluences = self._update_confluences_for_flowline(flowline, add=True)

            # Update topology if requested
            if update_topology:
                self._update_local_topology(flowline_id, affected_vertices, affected_confluences)

            # Update stream order if requested
            if update_stream_order:
                self._update_local_stream_order(flowline_id, affected_vertices)

            # Record the change
            if self.state_manager:
                self.state_manager.record_change(
                    GraphChangeType.ADD_FLOWLINE,
                    {flowline_id},
                    {
                        'flowline_id': flowline_id,
                        'affected_vertices': affected_vertices,
                        'affected_confluences': affected_confluences,
                        'execution_time': time.time() - start_time
                    }
                )

            logger.debug(f"Added flowline {flowline_id} incrementally in {time.time() - start_time:.4f}s")
            return True

        except Exception as e:
            logger.error(f"Failed to add flowline incrementally: {e}")
            return False

    def remove_single_flowline(self, flowline_id: int,
                              update_topology: bool = True,
                              update_stream_order: bool = False) -> bool:
        """
        Remove a single flowline from the graph with incremental updates.

        Args:
            flowline_id: ID of the flowline to remove
            update_topology: Whether to update topology after removal
            update_stream_order: Whether to recalculate stream order

        Returns:
            True if removal was successful
        """
        start_time = time.time()

        try:
            if flowline_id < 0 or flowline_id >= len(self.graph.aFlowline):
                logger.warning(f"Invalid flowline ID: {flowline_id}")
                return False

            flowline = self.graph.aFlowline[flowline_id]

            # Get affected elements before removal
            affected_vertices = self._get_flowline_vertices(flowline)
            affected_confluences = self._get_flowline_confluences(flowline)

            # Remove from flowline list (mark as None to preserve indices)
            self.graph.aFlowline[flowline_id] = None

            # Update vertices incrementally
            self._update_vertices_for_flowline(flowline, add=False)

            # Update confluences incrementally
            self._update_confluences_for_flowline(flowline, add=False)

            # Update topology if requested
            if update_topology:
                self._update_local_topology(flowline_id, affected_vertices, affected_confluences)

            # Update stream order if requested
            if update_stream_order:
                self._update_local_stream_order(flowline_id, affected_vertices)

            # Record the change
            if self.state_manager:
                self.state_manager.record_change(
                    GraphChangeType.REMOVE_FLOWLINE,
                    {flowline_id},
                    {
                        'flowline_id': flowline_id,
                        'affected_vertices': affected_vertices,
                        'affected_confluences': affected_confluences,
                        'execution_time': time.time() - start_time
                    }
                )

            logger.debug(f"Removed flowline {flowline_id} incrementally in {time.time() - start_time:.4f}s")
            return True

        except Exception as e:
            logger.error(f"Failed to remove flowline incrementally: {e}")
            return False

    def modify_single_flowline(self, flowline_id: int, modifications: Dict[str, Any],
                              update_topology: bool = True,
                              update_stream_order: bool = False) -> bool:
        """
        Modify a single flowline with incremental updates.

        Args:
            flowline_id: ID of the flowline to modify
            modifications: Dictionary of attribute changes
            update_topology: Whether to update topology after modification
            update_stream_order: Whether to recalculate stream order

        Returns:
            True if modification was successful
        """
        start_time = time.time()

        try:
            if flowline_id < 0 or flowline_id >= len(self.graph.aFlowline):
                logger.warning(f"Invalid flowline ID: {flowline_id}")
                return False

            flowline = self.graph.aFlowline[flowline_id]
            if flowline is None:
                logger.warning(f"Flowline {flowline_id} has been removed")
                return False

            # Store original values for potential rollback
            original_values = {}
            geometry_changed = False

            # Apply modifications
            for attr, value in modifications.items():
                if hasattr(flowline, attr):
                    original_values[attr] = getattr(flowline, attr)
                    setattr(flowline, attr, value)

                    # Check if geometry-related attributes changed
                    if attr in ['aVertex', 'dLength', 'dSlope']:
                        geometry_changed = True

            affected_vertices = set()
            affected_confluences = set()

            # If geometry changed, update related elements
            if geometry_changed:
                affected_vertices = self._get_flowline_vertices(flowline)
                affected_confluences = self._get_flowline_confluences(flowline)

                # Update topology if requested
                if update_topology:
                    self._update_local_topology(flowline_id, affected_vertices, affected_confluences)

                # Update stream order if requested
                if update_stream_order:
                    self._update_local_stream_order(flowline_id, affected_vertices)

            # Record the change
            if self.state_manager:
                self.state_manager.record_change(
                    GraphChangeType.MODIFY_FLOWLINE,
                    {flowline_id},
                    {
                        'flowline_id': flowline_id,
                        'modifications': modifications,
                        'geometry_changed': geometry_changed,
                        'affected_vertices': affected_vertices,
                        'affected_confluences': affected_confluences,
                        'execution_time': time.time() - start_time
                    }
                )

            logger.debug(f"Modified flowline {flowline_id} incrementally in {time.time() - start_time:.4f}s")
            return True

        except Exception as e:
            logger.error(f"Failed to modify flowline incrementally: {e}")
            return False

    def add_single_vertex(self, vertex: Any) -> bool:
        """
        Add a single vertex to the graph.

        Args:
            vertex: The vertex object to add

        Returns:
            True if addition was successful
        """
        start_time = time.time()

        try:
            # Assign vertex ID
            vertex_id = len(self.graph.aVertex)
            vertex.lVertexIndex = vertex_id

            # Add to vertex list
            self.graph.aVertex.append(vertex)

            # Record the change
            if self.state_manager:
                self.state_manager.record_change(
                    GraphChangeType.ADD_VERTEX,
                    {vertex_id},
                    {
                        'vertex_id': vertex_id,
                        'execution_time': time.time() - start_time
                    }
                )

            logger.debug(f"Added vertex {vertex_id} incrementally")
            return True

        except Exception as e:
            logger.error(f"Failed to add vertex incrementally: {e}")
            return False

    def remove_single_vertex(self, vertex_id: int) -> bool:
        """
        Remove a single vertex from the graph.

        Args:
            vertex_id: ID of the vertex to remove

        Returns:
            True if removal was successful
        """
        start_time = time.time()

        try:
            if vertex_id < 0 or vertex_id >= len(self.graph.aVertex):
                logger.warning(f"Invalid vertex ID: {vertex_id}")
                return False

            # Check if vertex is used by any flowlines
            vertex_in_use = self._is_vertex_in_use(vertex_id)
            if vertex_in_use:
                logger.warning(f"Cannot remove vertex {vertex_id}: still in use by flowlines")
                return False

            # Remove vertex (mark as None to preserve indices)
            self.graph.aVertex[vertex_id] = None

            # Record the change
            if self.state_manager:
                self.state_manager.record_change(
                    GraphChangeType.REMOVE_VERTEX,
                    {vertex_id},
                    {
                        'vertex_id': vertex_id,
                        'execution_time': time.time() - start_time
                    }
                )

            logger.debug(f"Removed vertex {vertex_id} incrementally")
            return True

        except Exception as e:
            logger.error(f"Failed to remove vertex incrementally: {e}")
            return False

    def update_confluence_incrementally(self, confluence_id: int,
                                        modifications: Dict[str, Any]) -> bool:
        """
        Update a confluence incrementally with new modifications.

        Args:
            confluence_id: ID of the confluence to update
            modifications: Dictionary of attribute changes

        Returns:
            True if update was successful
        """
        start_time = time.time()

        try:
            if confluence_id < 0 or confluence_id >= len(self.graph.aConfluence):
                logger.warning(f"Invalid confluence ID: {confluence_id}")
                return False

            confluence = self.graph.aConfluence[confluence_id]
            if confluence is None:
                logger.warning(f"Confluence {confluence_id} has been removed")
                return False

            # Apply modifications
            for attr, value in modifications.items():
                if hasattr(confluence, attr):
                    setattr(confluence, attr, value)

            # Record the change
            if self.state_manager:
                self.state_manager.record_change(
                    GraphChangeType.MODIFY_CONFLUENCE,
                    {confluence_id},
                    {
                        'confluence_id': confluence_id,
                        'modifications': modifications,
                        'execution_time': time.time() - start_time
                    }
                )

            logger.debug(f"Updated confluence {confluence_id} incrementally")
            return True

        except Exception as e:
            logger.error(f"Failed to update confluence incrementally: {e}")
            return False

    # Helper methods for incremental updates

    def _update_vertices_for_flowline(self, flowline: Any, add: bool = True) -> Set[int]:
        """
        Update vertices affected by flowline addition/removal.

        Args:
            flowline: The flowline object
            add: True for addition, False for removal

        Returns:
            Set of affected vertex IDs
        """
        affected_vertices = set()

        try:
            if hasattr(flowline, 'aVertex') and flowline.aVertex:
                for vertex in flowline.aVertex:
                    if hasattr(vertex, 'lVertexIndex'):
                        vertex_id = vertex.lVertexIndex
                        affected_vertices.add(vertex_id)

                        if add:
                            # Add vertex if not already in graph
                            if vertex_id >= len(self.graph.aVertex):
                                # Extend vertex list if needed
                                while len(self.graph.aVertex) <= vertex_id:
                                    self.graph.aVertex.append(None)
                                self.graph.aVertex[vertex_id] = vertex
                        else:
                            # Check if vertex is still used by other flowlines
                            if not self._is_vertex_in_use(vertex_id, exclude_flowline=flowline):
                                # Mark vertex as unused (don't remove to preserve indices)
                                if vertex_id < len(self.graph.aVertex):
                                    self.graph.aVertex[vertex_id] = None

        except Exception as e:
            logger.error(f"Error updating vertices for flowline: {e}")

        return affected_vertices

    def _update_confluences_for_flowline(self, flowline: Any, add: bool = True) -> Set[int]:
        """
        Update confluences affected by flowline addition/removal.

        Args:
            flowline: The flowline object
            add: True for addition, False for removal

        Returns:
            Set of affected confluence IDs
        """
        affected_confluences = set()

        try:
            # Find confluences at flowline endpoints
            if hasattr(flowline, 'aVertex') and flowline.aVertex:
                start_vertex = flowline.aVertex[0] if flowline.aVertex else None
                end_vertex = flowline.aVertex[-1] if len(flowline.aVertex) > 1 else None

                for vertex in [start_vertex, end_vertex]:
                    if vertex and hasattr(vertex, 'lVertexIndex'):
                        confluence_id = self._find_confluence_at_vertex(vertex.lVertexIndex)
                        if confluence_id is not None:
                            affected_confluences.add(confluence_id)

                            if add:
                                # Update confluence to include new flowline
                                self._add_flowline_to_confluence(confluence_id, flowline)
                            else:
                                # Update confluence to remove flowline
                                self._remove_flowline_from_confluence(confluence_id, flowline)

        except Exception as e:
            logger.error(f"Error updating confluences for flowline: {e}")

        return affected_confluences

    def _update_local_topology(self, flowline_id: int,
                              affected_vertices: Set[int],
                              affected_confluences: Set[int]) -> None:
        """
        Update local topology for affected elements.

        Args:
            flowline_id: ID of the modified flowline
            affected_vertices: Set of affected vertex IDs
            affected_confluences: Set of affected confluence IDs
        """
        try:
            # Update topology for affected vertices
            for vertex_id in affected_vertices:
                if vertex_id < len(self.graph.aVertex) and self.graph.aVertex[vertex_id]:
                    self._update_vertex_connections(vertex_id)

            # Update topology for affected confluences
            for confluence_id in affected_confluences:
                if confluence_id < len(self.graph.aConfluence) and self.graph.aConfluence[confluence_id]:
                    self._update_confluence_topology(confluence_id)

            # Update downstream/upstream relationships
            self._update_flow_relationships(flowline_id, affected_vertices)

        except Exception as e:
            logger.error(f"Error updating local topology: {e}")

    def _update_local_stream_order(self, flowline_id: int, affected_vertices: Set[int]) -> None:
        """
        Update stream order for locally affected flowlines.

        Args:
            flowline_id: ID of the modified flowline
            affected_vertices: Set of affected vertex IDs
        """
        try:
            # Find all flowlines connected to affected vertices
            affected_flowlines = set()
            affected_flowlines.add(flowline_id)

            for vertex_id in affected_vertices:
                connected_flowlines = self._get_flowlines_at_vertex(vertex_id)
                affected_flowlines.update(connected_flowlines)

            # Recalculate stream order for affected flowlines
            for fid in affected_flowlines:
                if fid < len(self.graph.aFlowline) and self.graph.aFlowline[fid]:
                    self._calculate_flowline_stream_order(fid)

        except Exception as e:
            logger.error(f"Error updating local stream order: {e}")

    def _get_flowline_vertices(self, flowline: Any) -> Set[int]:
        """Get all vertex IDs associated with a flowline."""
        vertices = set()
        try:
            if hasattr(flowline, 'aVertex') and flowline.aVertex:
                for vertex in flowline.aVertex:
                    if hasattr(vertex, 'lVertexIndex'):
                        vertices.add(vertex.lVertexIndex)
        except Exception as e:
            logger.error(f"Error getting flowline vertices: {e}")
        return vertices

    def _get_flowline_confluences(self, flowline: Any) -> Set[int]:
        """Get all confluence IDs associated with a flowline."""
        confluences = set()
        try:
            if hasattr(flowline, 'aVertex') and flowline.aVertex:
                start_vertex = flowline.aVertex[0] if flowline.aVertex else None
                end_vertex = flowline.aVertex[-1] if len(flowline.aVertex) > 1 else None

                for vertex in [start_vertex, end_vertex]:
                    if vertex and hasattr(vertex, 'lVertexIndex'):
                        confluence_id = self._find_confluence_at_vertex(vertex.lVertexIndex)
                        if confluence_id is not None:
                            confluences.add(confluence_id)
        except Exception as e:
            logger.error(f"Error getting flowline confluences: {e}")
        return confluences

    def _is_vertex_in_use(self, vertex_id: int, exclude_flowline: Any = None) -> bool:
        """Check if a vertex is still used by any flowlines."""
        try:
            for flowline in self.graph.aFlowline:
                if flowline is None or flowline == exclude_flowline:
                    continue

                if hasattr(flowline, 'aVertex') and flowline.aVertex:
                    for vertex in flowline.aVertex:
                        if hasattr(vertex, 'lVertexIndex') and vertex.lVertexIndex == vertex_id:
                            return True
        except Exception as e:
            logger.error(f"Error checking vertex usage: {e}")
        return False

    def _find_confluence_at_vertex(self, vertex_id: int) -> Optional[int]:
        """Find confluence ID at a given vertex."""
        try:
            for i, confluence in enumerate(self.graph.aConfluence):
                if confluence is None:
                    continue

                if hasattr(confluence, 'lVertexIndex') and confluence.lVertexIndex == vertex_id:
                    return i
        except Exception as e:
            logger.error(f"Error finding confluence at vertex: {e}")
        return None

    def _add_flowline_to_confluence(self, confluence_id: int, flowline: Any) -> None:
        """Add a flowline to a confluence."""
        try:
            if confluence_id < len(self.graph.aConfluence):
                confluence = self.graph.aConfluence[confluence_id]
                if confluence and hasattr(confluence, 'aFlowline'):
                    if flowline not in confluence.aFlowline:
                        confluence.aFlowline.append(flowline)
        except Exception as e:
            logger.error(f"Error adding flowline to confluence: {e}")

    def _remove_flowline_from_confluence(self, confluence_id: int, flowline: Any) -> None:
        """Remove a flowline from a confluence."""
        try:
            if confluence_id < len(self.graph.aConfluence):
                confluence = self.graph.aConfluence[confluence_id]
                if confluence and hasattr(confluence, 'aFlowline'):
                    if flowline in confluence.aFlowline:
                        confluence.aFlowline.remove(flowline)
        except Exception as e:
            logger.error(f"Error removing flowline from confluence: {e}")

    def _update_vertex_connections(self, vertex_id: int) -> None:
        """Update connections for a specific vertex."""
        try:
            if vertex_id < len(self.graph.aVertex):
                vertex = self.graph.aVertex[vertex_id]
                if vertex:
                    # Update vertex connections based on current flowlines
                    connected_flowlines = self._get_flowlines_at_vertex(vertex_id)
                    if hasattr(vertex, 'aFlowline'):
                        vertex.aFlowline = [self.graph.aFlowline[fid] for fid in connected_flowlines
                                          if fid < len(self.graph.aFlowline) and self.graph.aFlowline[fid]]
        except Exception as e:
            logger.error(f"Error updating vertex connections: {e}")

    def _update_confluence_topology(self, confluence_id: int) -> None:
        """Update topology for a specific confluence."""
        try:
            if confluence_id < len(self.graph.aConfluence):
                confluence = self.graph.aConfluence[confluence_id]
                if confluence and hasattr(confluence, 'lVertexIndex'):
                    # Update confluence flowline connections
                    vertex_id = confluence.lVertexIndex
                    connected_flowlines = self._get_flowlines_at_vertex(vertex_id)
                    if hasattr(confluence, 'aFlowline'):
                        confluence.aFlowline = [self.graph.aFlowline[fid] for fid in connected_flowlines
                                              if fid < len(self.graph.aFlowline) and self.graph.aFlowline[fid]]
        except Exception as e:
            logger.error(f"Error updating confluence topology: {e}")

    def _update_flow_relationships(self, flowline_id: int, affected_vertices: Set[int]) -> None:
        """Update upstream/downstream relationships for affected flowlines."""
        try:
            # This is a simplified version - full implementation would depend on
            # the specific graph structure and flow direction algorithms
            for vertex_id in affected_vertices:
                connected_flowlines = self._get_flowlines_at_vertex(vertex_id)
                for fid in connected_flowlines:
                    if fid < len(self.graph.aFlowline) and self.graph.aFlowline[fid]:
                        # Update flow relationships based on vertex positions
                        self._update_flowline_relationships(fid)
        except Exception as e:
            logger.error(f"Error updating flow relationships: {e}")

    def _get_flowlines_at_vertex(self, vertex_id: int) -> Set[int]:
        """Get all flowline IDs that use a specific vertex."""
        flowlines = set()
        try:
            for i, flowline in enumerate(self.graph.aFlowline):
                if flowline is None:
                    continue

                if hasattr(flowline, 'aVertex') and flowline.aVertex:
                    for vertex in flowline.aVertex:
                        if hasattr(vertex, 'lVertexIndex') and vertex.lVertexIndex == vertex_id:
                            flowlines.add(i)
                            break
        except Exception as e:
            logger.error(f"Error getting flowlines at vertex: {e}")
        return flowlines

    def _calculate_flowline_stream_order(self, flowline_id: int) -> None:
        """Calculate stream order for a specific flowline."""
        try:
            if flowline_id < len(self.graph.aFlowline):
                flowline = self.graph.aFlowline[flowline_id]
                if flowline:
                    # Simplified stream order calculation
                    # Full implementation would use proper Strahler or Shreve ordering
                    upstream_orders = self._get_upstream_orders(flowline_id)
                    if upstream_orders:
                        max_order = max(upstream_orders)
                        same_max_count = sum(1 for order in upstream_orders if order == max_order)
                        if same_max_count > 1:
                            stream_order = max_order + 1
                        else:
                            stream_order = max_order
                    else:
                        stream_order = 1

                    if hasattr(flowline, 'iStreamOrder'):
                        flowline.iStreamOrder = stream_order
        except Exception as e:
            logger.error(f"Error calculating flowline stream order: {e}")

    def _get_upstream_orders(self, flowline_id: int) -> List[int]:
        """Get stream orders of upstream flowlines."""
        orders = []
        try:
            if flowline_id < len(self.graph.aFlowline):
                flowline = self.graph.aFlowline[flowline_id]
                if flowline and hasattr(flowline, 'aVertex') and flowline.aVertex:
                    # Get upstream flowlines at start vertex
                    start_vertex = flowline.aVertex[0]
                    if hasattr(start_vertex, 'lVertexIndex'):
                        upstream_flowlines = self._get_upstream_flowlines(flowline_id, start_vertex.lVertexIndex)
                        for upstream_id in upstream_flowlines:
                            if upstream_id < len(self.graph.aFlowline):
                                upstream_flowline = self.graph.aFlowline[upstream_id]
                                if upstream_flowline and hasattr(upstream_flowline, 'iStreamOrder'):
                                    orders.append(upstream_flowline.iStreamOrder)
        except Exception as e:
            logger.error(f"Error getting upstream orders: {e}")
        return orders

    def _get_upstream_flowlines(self, flowline_id: int, vertex_id: int) -> Set[int]:
        """Get upstream flowline IDs at a vertex."""
        upstream = set()
        try:
            connected_flowlines = self._get_flowlines_at_vertex(vertex_id)
            for fid in connected_flowlines:
                if fid != flowline_id and fid < len(self.graph.aFlowline):
                    flowline = self.graph.aFlowline[fid]
                    if flowline and hasattr(flowline, 'aVertex') and flowline.aVertex:
                        # Check if this flowline flows into the vertex
                        end_vertex = flowline.aVertex[-1] if flowline.aVertex else None
                        if end_vertex and hasattr(end_vertex, 'lVertexIndex'):
                            if end_vertex.lVertexIndex == vertex_id:
                                upstream.add(fid)
        except Exception as e:
            logger.error(f"Error getting upstream flowlines: {e}")
        return upstream

    def _update_flowline_relationships(self, flowline_id: int) -> None:
        """Update upstream/downstream relationships for a flowline."""
        try:
            if flowline_id < len(self.graph.aFlowline):
                flowline = self.graph.aFlowline[flowline_id]
                if flowline and hasattr(flowline, 'aVertex') and flowline.aVertex:
                    # Update upstream flowlines
                    start_vertex = flowline.aVertex[0]
                    if hasattr(start_vertex, 'lVertexIndex'):
                        upstream_ids = self._get_upstream_flowlines(flowline_id, start_vertex.lVertexIndex)
                        if hasattr(flowline, 'aUpstreamFlowline'):
                            flowline.aUpstreamFlowline = [self.graph.aFlowline[uid] for uid in upstream_ids
                                                        if uid < len(self.graph.aFlowline) and self.graph.aFlowline[uid]]

                    # Update downstream flowlines
                    end_vertex = flowline.aVertex[-1] if len(flowline.aVertex) > 1 else None
                    if end_vertex and hasattr(end_vertex, 'lVertexIndex'):
                        downstream_ids = self._get_downstream_flowlines(flowline_id, end_vertex.lVertexIndex)
                        if hasattr(flowline, 'aDownstreamFlowline'):
                            flowline.aDownstreamFlowline = [self.graph.aFlowline[did] for did in downstream_ids
                                                          if did < len(self.graph.aFlowline) and self.graph.aFlowline[did]]
        except Exception as e:
            logger.error(f"Error updating flowline relationships: {e}")

    def _get_downstream_flowlines(self, flowline_id: int, vertex_id: int) -> Set[int]:
        """Get downstream flowline IDs at a vertex."""
        downstream = set()
        try:
            connected_flowlines = self._get_flowlines_at_vertex(vertex_id)
            for fid in connected_flowlines:
                if fid != flowline_id and fid < len(self.graph.aFlowline):
                    flowline = self.graph.aFlowline[fid]
                    if flowline and hasattr(flowline, 'aVertex') and flowline.aVertex:
                        # Check if this flowline starts from the vertex
                        start_vertex = flowline.aVertex[0] if flowline.aVertex else None
                        if start_vertex and hasattr(start_vertex, 'lVertexIndex'):
                            if start_vertex.lVertexIndex == vertex_id:
                                downstream.add(fid)
        except Exception as e:
            logger.error(f"Error getting downstream flowlines: {e}")
        return downstream
