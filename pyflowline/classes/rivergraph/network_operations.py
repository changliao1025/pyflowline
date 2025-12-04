"""
Network operations module for pyrivergraph.

This module contains functions for network modification and simplification operations
including removal of disconnected flowlines, braided channels, parallel paths, cycles,
small rivers, and merging of linear segments.
"""

import logging
import os
from typing import List, Dict, Set, Tuple, Optional
from collections import defaultdict, deque
import numpy as np

from pyflowline.classes.flowline import pyflowline
from pyflowline.classes.vertex import pyvertex
from pyflowline.formats.export_flowline import export_flowline_to_geojson
from pyflowline.formats.export_vertex import export_vertex_to_geojson

from .state_management import GraphChangeType

logger = logging.getLogger(__name__)


class NetworkOperations:
    """
    Network operations for river network modification and simplification.

    This class provides methods for various network operations that modify
    the structure of the river network, including removal of disconnected
    flowlines, braided channels, parallel paths, cycles, and small rivers.
    """

    def __init__(self, river_graph):
        """
        Initialize network operations with reference to the river graph.

        Args:
            river_graph: The pyrivergraph instance to operate on
        """
        self.river_graph = river_graph

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

        target_outlet = outlet_vertex if outlet_vertex is not None else self.river_graph.pVertex_outlet

        if target_outlet is None:
            logger.warning("No outlet vertex provided and none stored during initialization")
            return flowlines

        outlet_vertex_id = None
        if target_outlet == self.river_graph.pVertex_outlet and self.river_graph.pVertex_outlet_id is not None:
            outlet_vertex_id = self.river_graph.pVertex_outlet_id
        else:
            for vertex_id, vertex in self.river_graph.id_to_vertex.items():
                if vertex == target_outlet:
                    outlet_vertex_id = vertex_id
                    break

        if outlet_vertex_id is None:
            logger.warning("Outlet vertex not found in network graph")
            return flowlines

        logger.info(f"Removing disconnected flowlines using outlet vertex {outlet_vertex_id}")

        outlet_reachable_vertices = self._find_outlet_reachable_vertices(outlet_vertex_id)

        reachable_flowlines = set()

        for flowline_idx, (start_id, end_id) in self.river_graph.aFlowline_edges.items():
            if start_id in outlet_reachable_vertices and end_id in outlet_reachable_vertices:
                if self._flowline_contributes_to_outlet(flowline_idx, start_id, end_id, outlet_vertex_id, outlet_reachable_vertices):
                    reachable_flowlines.add(flowline_idx)

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
        self.river_graph.state_manager.record_change(
            GraphChangeType.REMOVE_BRAIDED_RIVERS,
            {'network_size': len(self.river_graph.aFlowline)}
        )

        braided_channels = self.river_graph.find_braided_channels()
        if not braided_channels:
            logger.info("No braided channels found in the network")
            return self.river_graph.aFlowline

        logger.info(f"Found {len(braided_channels)} braided channel groups to resolve")

        flowlines_to_remove = set()

        for (start_id, end_id), flowline_indices in braided_channels.items():
            for flowline_idx in flowline_indices[1:]:
                flowlines_to_remove.add(flowline_idx)
                logger.debug(f"Marking flowline {flowline_idx} for removal (braided channel)")

        simplified_flowlines = [
            flowline for idx, flowline in enumerate(self.river_graph.aFlowline)
            if idx not in flowlines_to_remove
        ]

        logger.info(f"Removed {len(flowlines_to_remove)} braided flowlines, "
                    f"resulting in {len(simplified_flowlines)} flowlines")

        self.river_graph._update_graph_flowlines(simplified_flowlines)

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
        if not hasattr(self.river_graph, 'aFlowline') or not self.river_graph.aFlowline:
            logger.warning('No flowlines available in class instance for parallel river removal')
            return []

        if len(self.river_graph.aFlowline) <= 1:
            logger.debug("Skipping parallel path resolution: insufficient flowlines")
            return self.river_graph.aFlowline.copy()

        self.river_graph.state_manager.record_change(
            GraphChangeType.REMOVE_PARALLEL_PATHS,
            {'network_size': len(self.river_graph.aFlowline)}
        )

        logger.info(f"Removing parallel rivers from {len(self.river_graph.aFlowline)} flowlines using graph-based approach")

        try:
            parallel_groups = self.river_graph.find_parallel_paths()

            if not parallel_groups:
                logger.debug("No parallel paths found")
                return self.river_graph.aFlowline.copy()

            logger.info(f"Found {len(parallel_groups)} parallel path groups")
        except Exception as e:
            logger.error(f"Error finding parallel paths: {e}")
            return self.river_graph.aFlowline.copy()

        flowlines_to_remove = set()

        for group in parallel_groups:
            paths = group['paths']

            if len(paths) <= 1:
                continue

            path_scores = []
            for path_idx, path_flowlines in enumerate(paths):
                path_total_score = 0.0
                path_length = 0.0
                path_max_order = 0
                path_total_area = 0.0
                valid_flowlines = 0

                for flowline_idx in path_flowlines:
                    if flowline_idx >= len(self.river_graph.aFlowline):
                        logger.warning(f"Invalid flowline index {flowline_idx}, skipping")
                        continue

                    flowline = self.river_graph.aFlowline[flowline_idx]
                    valid_flowlines += 1

                    stream_order = getattr(flowline, 'iStream_order', -1)
                    stream_order = stream_order if stream_order > 0 else 1

                    drainage_area = getattr(flowline, 'dDrainage_area', 0.0)
                    drainage_area = drainage_area if drainage_area > 0 else 0.0

                    length = getattr(flowline, 'dLength', 0.0)
                    length = length if length > 0 else 0.0

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

            path_scores.sort(reverse=True)
            best_score, best_path_idx, best_path_flowlines = path_scores[0]

            for score, path_idx, path_flowlines in path_scores[1:]:
                for flowline_idx in path_flowlines:
                    flowlines_to_remove.add(flowline_idx)

            logger.debug(f"Selected path {best_path_idx} with score {best_score:.2e} from {len(paths)} parallel paths "
                        f"(vertices {group['start_vertex']} -> {group['end_vertex']})")
            logger.debug(f"Keeping {len(best_path_flowlines)} flowlines, removing {sum(len(path_flowlines) for _, _, path_flowlines in path_scores[1:])} flowlines")

        try:
            result = [self.river_graph.aFlowline[i] for i in range(len(self.river_graph.aFlowline)) if i not in flowlines_to_remove]
            logger.info(f"Removed {len(flowlines_to_remove)} flowlines to resolve parallel paths")

            if len(result) == 0 and len(self.river_graph.aFlowline) > 0:
                logger.warning("All flowlines were removed during parallel path resolution, returning original")
                return self.river_graph.aFlowline.copy()

            self.river_graph._update_graph_flowlines(result)

            return result
        except Exception as e:
            logger.error(f"Error filtering flowlines during parallel path resolution: {e}")
            return self.river_graph.aFlowline.copy()

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
        if not hasattr(self.river_graph, 'aFlowline') or not self.river_graph.aFlowline:
            logger.warning('No flowlines available in class instance for cycle removal')
            return []

        if len(self.river_graph.aFlowline) <= 1:
            logger.debug("Skipping cycle removal: insufficient flowlines")
            return self.river_graph.aFlowline.copy()

        self.river_graph.state_manager.record_change(
            GraphChangeType.REMOVE_CYCLES,
            {'network_size': len(self.river_graph.aFlowline)}
        )

        logger.info(f"Removing cycles from {len(self.river_graph.aFlowline)} flowlines using graph-based approach")

        try:
            cycles = self.river_graph.detect_cycles()

            if not cycles:
                logger.debug("No cycles detected in network")
                return self.river_graph.aFlowline.copy()

            logger.info(f"Detected {len(cycles)} cycles in network")
        except Exception as e:
            logger.error(f"Error detecting cycles: {e}")
            return self.river_graph.aFlowline.copy()

        flowlines_to_remove = set()

        for cycle_vertices in cycles:
            if len(cycle_vertices) < 3:
                continue

            cycle_flowlines = []
            for i in range(len(cycle_vertices) - 1):
                start_vertex_id = cycle_vertices[i]
                end_vertex_id = cycle_vertices[i + 1]

                for flowline_idx, (s_id, e_id) in self.river_graph.aFlowline_edges.items():
                    if s_id == start_vertex_id and e_id == end_vertex_id:
                        if flowline_idx < len(self.river_graph.aFlowline):
                            cycle_flowlines.append((flowline_idx, self.river_graph.aFlowline[flowline_idx]))

            if cycle_flowlines:
                def cycle_priority(flowline_data) -> Tuple:
                    flowline_idx, flowline = flowline_data
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

        try:
            result = [self.river_graph.aFlowline[i] for i in range(len(self.river_graph.aFlowline)) if i not in flowlines_to_remove]
            logger.info(f"Removed {len(flowlines_to_remove)} flowlines to break cycles")

            if len(result) == 0 and len(self.river_graph.aFlowline) > 0:
                logger.warning("All flowlines were removed during cycle removal, returning original")
                return self.river_graph.aFlowline.copy()

            self.river_graph._update_graph_flowlines(result)

            return result
        except Exception as e:
            logger.error(f"Error filtering flowlines during cycle removal: {e}")
            return self.river_graph.aFlowline.copy()

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
        aFlowline_current = self.river_graph.aFlowline.copy()

        for i in range(nIterations):
            sStep = "{:02d}".format(i+1)
            print(f'Iteration {sStep}: removing small rivers with threshold {dThreshold_small_river}')

            aFlowline_filtered = self._remove_small_rivers(aFlowline_current, dThreshold_small_river)

            if iFlag_debug == 1 and sWorkspace_output_basin is not None:
                sFilename_out = f'flowline_large_{sStep}_before_intersect.geojson'
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_geojson(aFlowline_filtered, sFilename_out)

            self.river_graph._update_graph_flowlines(aFlowline_filtered)

            aFlowline_filtered = self.river_graph.update_head_water_stream_order()

            aFlowline_merged = self.merge_linear_segments(aFlowline_filtered)

            if iFlag_debug == 1 and sWorkspace_output_basin is not None:
                sFilename_out = f'flowline_merge_{sStep}_before_intersect.geojson'
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_geojson(aFlowline_merged, sFilename_out)

            self.river_graph._update_graph_flowlines(aFlowline_merged)

            aFlowline_current = self.river_graph.update_head_water_stream_order()

            if len(aFlowline_current) == 1:
                break

        return aFlowline_current

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

        self.river_graph.state_manager.record_change(
            GraphChangeType.MERGE_LINEAR_SEGMENTS,
            {'network_size': len(flowlines)}
        )

        original_flowlines = self.river_graph.aFlowline
        self.river_graph.aFlowline = flowlines
        self.river_graph._build_graph()

        try:
            linear_segments = self.river_graph.find_linear_segments()

            if not linear_segments:
                logger.info("No linear segments found for merging")
                return flowlines

            logger.info(f"Found {len(linear_segments)} linear segments to merge")

            merged_flowlines = []
            processed_indices = set()

            for segment in linear_segments:
                if len(segment) < 2:
                    continue

                try:
                    segment_flowlines = [flowlines[idx] for idx in segment if idx < len(flowlines)]

                    if len(segment_flowlines) < 2:
                        continue

                    merged_flowline = self._merge_flowline_segment(segment_flowlines)
                    if merged_flowline is not None:
                        merged_flowlines.append(merged_flowline)
                        processed_indices.update(segment)
                        logger.debug(f"Merged {len(segment_flowlines)} flowlines into single flowline")

                except Exception as e:
                    logger.error(f"Error merging segment {segment}: {e}")
                    continue

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
            self.river_graph.aFlowline = original_flowlines
            self.river_graph._build_graph()

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

        reverse_adjacency = defaultdict(list)
        for start_id, neighbors in self.river_graph.adjacency_list.items():
            for end_id, flowline_idx in neighbors:
                reverse_adjacency[end_id].append((start_id, flowline_idx))

        while queue:
            current_id = queue.popleft()

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
        return start_id in outlet_reachable and end_id in outlet_reachable

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
            merged_flowline = segment_flowlines[0]

            merged_flowline.pVertex_end = segment_flowlines[-1].pVertex_end

            total_length = sum(getattr(fl, 'dLength', 0.0) for fl in segment_flowlines)
            merged_flowline.dLength = total_length

            max_order = max(getattr(fl, 'iStream_order', 1) for fl in segment_flowlines)
            merged_flowline.iStream_order = max_order

            max_area = max(getattr(fl, 'dDrainage_area', 0.0) for fl in segment_flowlines)
            merged_flowline.dDrainage_area = max_area

            return merged_flowline

        except Exception as e:
            logger.error(f"Error merging flowline segment: {e}")
            return None