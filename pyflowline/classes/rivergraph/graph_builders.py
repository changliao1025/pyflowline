"""
Graph builders module for pyrivergraph.

This module contains utilities for graph construction, vertex management,
and graph update operations. It provides the core functionality for building
and maintaining the graph structure from flowline data.
"""

import logging
from typing import List, Dict, Set, Tuple, Optional
from collections import defaultdict
import time

from pyflowline.classes.flowline import pyflowline
from pyflowline.classes.vertex import pyvertex

from .state_management import GraphUpdateStrategy

logger = logging.getLogger(__name__)


class GraphBuilders:
    """
    Graph construction and management utilities for river networks.

    This class provides methods for building graph structures from flowlines,
    managing vertices, and updating graph topology efficiently.
    """

    def __init__(self, river_graph):
        """
        Initialize graph builders with reference to the river graph.

        Args:
            river_graph: The pyrivergraph instance to operate on
        """
        self.river_graph = river_graph

    def build_graph(self):
        """
        Build the graph structure from flowlines.
        This is the core method that constructs the adjacency list representation.
        """
        self.river_graph.vertex_to_id.clear()
        self.river_graph.id_to_vertex.clear()
        self.river_graph.aFlowline_edges.clear()
        self.river_graph.adjacency_list.clear()
        self.river_graph.in_degree.clear()
        self.river_graph.out_degree.clear()
        self.river_graph.aVertex.clear()

        vertex_counter = 0

        for flowline_idx, flowline in enumerate(self.river_graph.aFlowline):
            if flowline is None:
                continue

            start_vertex = flowline.pVertex_start
            end_vertex = flowline.pVertex_end

            if start_vertex not in self.river_graph.vertex_to_id:
                self.river_graph.vertex_to_id[start_vertex] = vertex_counter
                self.river_graph.id_to_vertex[vertex_counter] = start_vertex
                self.river_graph.aVertex.append(start_vertex)
                start_vertex.lVertexID = vertex_counter
                vertex_counter += 1

            if end_vertex not in self.river_graph.vertex_to_id:
                self.river_graph.vertex_to_id[end_vertex] = vertex_counter
                self.river_graph.id_to_vertex[vertex_counter] = end_vertex
                self.river_graph.aVertex.append(end_vertex)
                end_vertex.lVertexID = vertex_counter
                vertex_counter += 1

            start_id = self.river_graph.vertex_to_id[start_vertex]
            end_id = self.river_graph.vertex_to_id[end_vertex]

            self.river_graph.aFlowline_edges[flowline_idx] = (start_id, end_id)

            self.river_graph.adjacency_list[start_id].append((end_id, flowline_idx))

            self.river_graph.out_degree[start_id] += 1
            self.river_graph.in_degree[end_id] += 1

        logger.debug(f"Built graph with {len(self.river_graph.id_to_vertex)} vertices and {len(self.river_graph.aFlowline_edges)} edges")

    def set_outlet_vertex_id(self):
        """Set the outlet vertex ID if outlet vertex is provided."""
        if self.river_graph.pVertex_outlet is not None:
            self.river_graph.pVertex_outlet_id = self.river_graph.vertex_to_id.get(self.river_graph.pVertex_outlet)
            if self.river_graph.pVertex_outlet_id is not None:
                logger.debug(f"Set outlet vertex ID: {self.river_graph.pVertex_outlet_id}")
            else:
                logger.warning("Outlet vertex not found in graph")

    def update_graph_flowlines(self, new_flowlines: List[pyflowline]):
        """
        Update the graph with a new set of flowlines using smart update strategy.

        Args:
            new_flowlines: New list of flowlines to update the graph with
        """
        strategy = self.river_graph.state_manager.get_update_strategy(len(self.river_graph.aFlowline))

        start_time = time.time()

        if strategy == GraphUpdateStrategy.INCREMENTAL and len(new_flowlines) > 0:
            success = self.try_incremental_update(new_flowlines)
            if not success:
                self.river_graph.aFlowline = new_flowlines
                self.build_graph()
        else:
            self.river_graph.aFlowline = new_flowlines
            self.build_graph()

        duration = time.time() - start_time
        self.river_graph.state_manager.record_update_performance(strategy, duration)
        self.river_graph.state_manager.clear_pending_changes()

    def try_incremental_update(self, new_flowlines: List[pyflowline]) -> bool:
        """
        Attempt to update the graph incrementally.

        Args:
            new_flowlines: New flowlines to incorporate

        Returns:
            bool: True if incremental update succeeded, False if rebuild needed
        """
        try:
            size_diff = abs(len(new_flowlines) - len(self.river_graph.aFlowline))
            if size_diff > len(self.river_graph.aFlowline) * 0.1:
                return False

            self.river_graph.aFlowline = new_flowlines

            self.build_graph()

            return True

        except Exception as e:
            logger.error(f"Incremental update failed: {e}")
            return False

    def add_vertex_if_new(self, vertex: pyvertex) -> int:
        """
        Add vertex to graph if it doesn't exist, return its ID.

        Args:
            vertex: Vertex to add

        Returns:
            int: Vertex ID
        """
        if vertex not in self.river_graph.vertex_to_id:
            vertex_id = len(self.river_graph.id_to_vertex)
            self.river_graph.vertex_to_id[vertex] = vertex_id
            self.river_graph.id_to_vertex[vertex_id] = vertex
            self.river_graph.aVertex.append(vertex)
            vertex.lVertexID = vertex_id
            return vertex_id
        else:
            return self.river_graph.vertex_to_id[vertex]

    def rebuild_graph_from_flowlines(self, flowlines: List[pyflowline]):
        """
        Completely rebuild the graph from a new set of flowlines.

        Args:
            flowlines: List of flowlines to build graph from
        """
        self.river_graph.aFlowline = flowlines

        self.build_graph()

        if self.river_graph.pVertex_outlet is not None:
            self.set_outlet_vertex_id()

        self.river_graph.state_manager.clear_pending_changes()

        logger.info(f"Rebuilt graph from {len(flowlines)} flowlines")

    def validate_graph_structure(self) -> Dict[str, any]:
        """
        Validate the internal consistency of the graph structure.

        Returns:
            Dictionary containing validation results
        """
        validation_results = {
            'is_valid': True,
            'issues': [],
            'statistics': {}
        }

        try:
            if len(self.river_graph.vertex_to_id) != len(self.river_graph.id_to_vertex):
                validation_results['issues'].append("Vertex mapping size mismatch")
                validation_results['is_valid'] = False

            for idx, flowline in enumerate(self.river_graph.aFlowline):
                if flowline is None:
                    continue

                start_vertex = flowline.pVertex_start
                end_vertex = flowline.pVertex_end

                if start_vertex not in self.river_graph.vertex_to_id:
                    validation_results['issues'].append(f"Flowline {idx} start vertex not in vertex mapping")
                    validation_results['is_valid'] = False

                if end_vertex not in self.river_graph.vertex_to_id:
                    validation_results['issues'].append(f"Flowline {idx} end vertex not in vertex mapping")
                    validation_results['is_valid'] = False

            for flowline_idx, (start_id, end_id) in self.river_graph.aFlowline_edges.items():
                if start_id not in self.river_graph.id_to_vertex:
                    validation_results['issues'].append(f"Edge {flowline_idx} references invalid start vertex {start_id}")
                    validation_results['is_valid'] = False

                if end_id not in self.river_graph.id_to_vertex:
                    validation_results['issues'].append(f"Edge {flowline_idx} references invalid end vertex {end_id}")
                    validation_results['is_valid'] = False

            validation_results['statistics'] = {
                'total_vertices': len(self.river_graph.id_to_vertex),
                'total_edges': len(self.river_graph.aFlowline_edges),
                'total_flowlines': len([fl for fl in self.river_graph.aFlowline if fl is not None]),
                'vertex_mapping_size': len(self.river_graph.vertex_to_id)
            }

        except Exception as e:
            validation_results['is_valid'] = False
            validation_results['issues'].append(f"Validation error: {e}")
            logger.error(f"Error during graph structure validation: {e}")

        return validation_results

    def get_graph_statistics(self) -> Dict[str, any]:
        """
        Get comprehensive statistics about the graph structure.

        Returns:
            Dictionary containing graph statistics
        """
        stats = {
            'vertices': {
                'total': len(self.river_graph.id_to_vertex),
                'sources': len([v_id for v_id in self.river_graph.id_to_vertex.keys()
                               if self.river_graph.in_degree[v_id] == 0]),
                'sinks': len([v_id for v_id in self.river_graph.id_to_vertex.keys()
                             if self.river_graph.out_degree[v_id] == 0]),
                'confluences': len([v_id for v_id in self.river_graph.id_to_vertex.keys()
                                  if self.river_graph.in_degree[v_id] > 1])
            },
            'edges': {
                'total': len(self.river_graph.aFlowline_edges),
                'active_flowlines': len([fl for fl in self.river_graph.aFlowline if fl is not None])
            },
            'connectivity': {
                'avg_out_degree': sum(self.river_graph.out_degree.values()) / max(len(self.river_graph.id_to_vertex), 1),
                'avg_in_degree': sum(self.river_graph.in_degree.values()) / max(len(self.river_graph.id_to_vertex), 1),
                'max_out_degree': max(self.river_graph.out_degree.values()) if self.river_graph.out_degree else 0,
                'max_in_degree': max(self.river_graph.in_degree.values()) if self.river_graph.in_degree else 0
            }
        }

        return stats

    def optimize_vertex_ordering(self):
        """
        Optimize vertex ordering for better cache performance and graph traversal.
        This method can reorder vertex IDs to improve spatial locality.
        """
        try:
            vertices_with_coords = []
            for vertex_id, vertex in self.river_graph.id_to_vertex.items():
                x = getattr(vertex, 'dLongitude_degree', 0.0)
                y = getattr(vertex, 'dLatitude_degree', 0.0)
                vertices_with_coords.append((vertex_id, vertex, x, y))

            vertices_with_coords.sort(key=lambda item: (item[2], item[3]))

            new_vertex_to_id = {}
            new_id_to_vertex = {}
            new_aVertex = []

            for new_id, (old_id, vertex, x, y) in enumerate(vertices_with_coords):
                new_vertex_to_id[vertex] = new_id
                new_id_to_vertex[new_id] = vertex
                new_aVertex.append(vertex)
                vertex.lVertexID = new_id

            new_aFlowline_edges = {}
            new_adjacency_list = defaultdict(list)
            new_in_degree = defaultdict(int)
            new_out_degree = defaultdict(int)

            for flowline_idx, (old_start_id, old_end_id) in self.river_graph.aFlowline_edges.items():
                old_start_vertex = self.river_graph.id_to_vertex[old_start_id]
                old_end_vertex = self.river_graph.id_to_vertex[old_end_id]

                new_start_id = new_vertex_to_id[old_start_vertex]
                new_end_id = new_vertex_to_id[old_end_vertex]

                new_aFlowline_edges[flowline_idx] = (new_start_id, new_end_id)
                new_adjacency_list[new_start_id].append((new_end_id, flowline_idx))
                new_out_degree[new_start_id] += 1
                new_in_degree[new_end_id] += 1

            self.river_graph.vertex_to_id = new_vertex_to_id
            self.river_graph.id_to_vertex = new_id_to_vertex
            self.river_graph.aVertex = new_aVertex
            self.river_graph.aFlowline_edges = new_aFlowline_edges
            self.river_graph.adjacency_list = new_adjacency_list
            self.river_graph.in_degree = new_in_degree
            self.river_graph.out_degree = new_out_degree

            if self.river_graph.pVertex_outlet is not None:
                self.set_outlet_vertex_id()

            logger.info("Optimized vertex ordering for improved cache performance")

        except Exception as e:
            logger.error(f"Error optimizing vertex ordering: {e}")

    def create_subgraph(self, flowline_indices: Set[int]) -> Dict[str, any]:
        """
        Create a subgraph containing only the specified flowlines.

        Args:
            flowline_indices: Set of flowline indices to include in subgraph

        Returns:
            Dictionary containing subgraph data
        """
        try:
            subgraph_flowlines = []
            for idx in flowline_indices:
                if idx < len(self.river_graph.aFlowline) and self.river_graph.aFlowline[idx] is not None:
                    subgraph_flowlines.append(self.river_graph.aFlowline[idx])

            subgraph_vertices = set()
            for flowline in subgraph_flowlines:
                subgraph_vertices.add(flowline.pVertex_start)
                subgraph_vertices.add(flowline.pVertex_end)

            subgraph_vertex_to_id = {}
            subgraph_id_to_vertex = {}
            subgraph_aVertex = []

            for vertex_id, vertex in enumerate(subgraph_vertices):
                subgraph_vertex_to_id[vertex] = vertex_id
                subgraph_id_to_vertex[vertex_id] = vertex
                subgraph_aVertex.append(vertex)

            subgraph_edges = {}
            subgraph_adjacency = defaultdict(list)
            subgraph_in_degree = defaultdict(int)
            subgraph_out_degree = defaultdict(int)

            for flowline_idx, flowline in enumerate(subgraph_flowlines):
                start_id = subgraph_vertex_to_id[flowline.pVertex_start]
                end_id = subgraph_vertex_to_id[flowline.pVertex_end]

                subgraph_edges[flowline_idx] = (start_id, end_id)
                subgraph_adjacency[start_id].append((end_id, flowline_idx))
                subgraph_out_degree[start_id] += 1
                subgraph_in_degree[end_id] += 1

            subgraph_data = {
                'flowlines': subgraph_flowlines,
                'vertices': subgraph_aVertex,
                'vertex_to_id': subgraph_vertex_to_id,
                'id_to_vertex': subgraph_id_to_vertex,
                'edges': subgraph_edges,
                'adjacency_list': dict(subgraph_adjacency),
                'in_degree': dict(subgraph_in_degree),
                'out_degree': dict(subgraph_out_degree),
                'statistics': {
                    'flowline_count': len(subgraph_flowlines),
                    'vertex_count': len(subgraph_vertices),
                    'edge_count': len(subgraph_edges)
                }
            }

            logger.info(f"Created subgraph with {len(subgraph_flowlines)} flowlines and {len(subgraph_vertices)} vertices")
            return subgraph_data

        except Exception as e:
            logger.error(f"Error creating subgraph: {e}")
            return {
                'flowlines': [],
                'vertices': [],
                'vertex_to_id': {},
                'id_to_vertex': {},
                'edges': {},
                'adjacency_list': {},
                'in_degree': {},
                'out_degree': {},
                'statistics': {'flowline_count': 0, 'vertex_count': 0, 'edge_count': 0},
                'error': str(e)
            }