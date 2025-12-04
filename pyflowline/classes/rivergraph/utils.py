"""
Utility functions for pyrivergraph.

This module provides shared utility functions used across the rivergraph package,
including geometric calculations, data structure helpers, and common algorithms.
"""

import math
import numpy as np
from typing import List, Tuple, Set, Dict, Any, Optional, Union
from collections import defaultdict, deque
import logging

logger = logging.getLogger(__name__)




def find_cycles_dfs(adjacency_dict: Dict[int, List[int]]) -> List[List[int]]:
    """
    Find all cycles in a directed graph using depth-first search.

    Args:
        adjacency_dict: Dictionary mapping node_id -> list of connected node_ids

    Returns:
        List of cycles, where each cycle is a list of node IDs
    """
    cycles = []
    visited = set()
    rec_stack = set()
    path = []

    def dfs(node: int) -> None:
        if node in rec_stack:
            # Found a cycle
            cycle_start = path.index(node)
            cycle = path[cycle_start:] + [node]
            cycles.append(cycle)
            return

        if node in visited:
            return

        visited.add(node)
        rec_stack.add(node)
        path.append(node)

        for neighbor in adjacency_dict.get(node, []):
            dfs(neighbor)

        rec_stack.remove(node)
        path.pop()

    # Start DFS from each unvisited node
    for node in adjacency_dict:
        if node not in visited:
            dfs(node)

    return cycles


def topological_sort(adjacency_dict: Dict[int, List[int]]) -> List[int]:
    """
    Perform topological sort on a directed acyclic graph.

    Args:
        adjacency_dict: Dictionary mapping node_id -> list of connected node_ids

    Returns:
        Topologically sorted list of node IDs

    Raises:
        ValueError: If the graph contains cycles
    """
    # Calculate in-degrees
    in_degree = defaultdict(int)
    all_nodes = set(adjacency_dict.keys())

    for node in adjacency_dict:
        for neighbor in adjacency_dict[node]:
            in_degree[neighbor] += 1
            all_nodes.add(neighbor)

    # Initialize queue with nodes having no incoming edges
    queue = deque([node for node in all_nodes if in_degree[node] == 0])
    result = []

    while queue:
        node = queue.popleft()
        result.append(node)

        # Remove this node from graph and update in-degrees
        for neighbor in adjacency_dict.get(node, []):
            in_degree[neighbor] -= 1
            if in_degree[neighbor] == 0:
                queue.append(neighbor)

    # Check if all nodes were processed (no cycles)
    if len(result) != len(all_nodes):
        raise ValueError("Graph contains cycles - topological sort not possible")

    return result


def find_strongly_connected_components(adjacency_dict: Dict[int, List[int]]) -> List[List[int]]:
    """
    Find strongly connected components using Tarjan's algorithm.

    Args:
        adjacency_dict: Dictionary mapping node_id -> list of connected node_ids

    Returns:
        List of strongly connected components, each as a list of node IDs
    """
    index_counter = [0]
    stack = []
    lowlinks = {}
    index = {}
    on_stack = {}
    components = []

    def strongconnect(node: int) -> None:
        # Set the depth index for this node to the smallest unused index
        index[node] = index_counter[0]
        lowlinks[node] = index_counter[0]
        index_counter[0] += 1
        stack.append(node)
        on_stack[node] = True

        # Consider successors of node
        for neighbor in adjacency_dict.get(node, []):
            if neighbor not in index:
                # Successor has not yet been visited; recurse on it
                strongconnect(neighbor)
                lowlinks[node] = min(lowlinks[node], lowlinks[neighbor])
            elif on_stack.get(neighbor, False):
                # Successor is in stack and hence in the current SCC
                lowlinks[node] = min(lowlinks[node], index[neighbor])

        # If node is a root node, pop the stack and create an SCC
        if lowlinks[node] == index[node]:
            component = []
            while True:
                w = stack.pop()
                on_stack[w] = False
                component.append(w)
                if w == node:
                    break
            components.append(component)

    # Run the algorithm for each unvisited node
    for node in adjacency_dict:
        if node not in index:
            strongconnect(node)

    return components


def merge_adjacent_segments(segments: List[List[Tuple[float, float]]],
                           tolerance: float = 1e-6) -> List[List[Tuple[float, float]]]:
    """
    Merge adjacent line segments that share endpoints.

    Args:
        segments: List of line segments, each as a list of points
        tolerance: Distance tolerance for considering points as identical

    Returns:
        List of merged segments
    """
    if not segments:
        return []

    # Build adjacency based on shared endpoints
    adjacency = defaultdict(list)

    for i, segment in enumerate(segments):
        if len(segment) < 2:
            continue

        start_point = segment[0]
        end_point = segment[-1]

        # Find segments that share endpoints
        for j, other_segment in enumerate(segments):
            if i == j or len(other_segment) < 2:
                continue

            other_start = other_segment[0]
            other_end = other_segment[-1]

            # Check if segments can be connected
            if calculate_distance_2d(end_point, other_start) < tolerance:
                adjacency[i].append(j)
            elif calculate_distance_2d(end_point, other_end) < tolerance:
                adjacency[i].append(-j - 1)  # Negative index indicates reverse
            elif calculate_distance_2d(start_point, other_start) < tolerance:
                adjacency[-i - 1].append(j)  # Reverse current segment
            elif calculate_distance_2d(start_point, other_end) < tolerance:
                adjacency[-i - 1].append(-j - 1)

    # Merge connected segments
    visited = set()
    merged_segments = []

    for i in range(len(segments)):
        if i in visited:
            continue

        # Start a new merged segment
        current_segment = segments[i][:]
        visited.add(i)

        # Follow the chain of connected segments
        current_index = i
        while current_index in adjacency:
            next_indices = adjacency[current_index]
            if not next_indices:
                break

            next_index = next_indices[0]  # Take first connection
            abs_next_index = abs(next_index) - 1 if next_index < 0 else next_index

            if abs_next_index in visited:
                break

            visited.add(abs_next_index)
            next_segment = segments[abs_next_index]

            if next_index < 0:
                # Reverse the next segment
                next_segment = next_segment[::-1]

            # Merge segments (remove duplicate endpoint)
            current_segment.extend(next_segment[1:])
            current_index = abs_next_index

        merged_segments.append(current_segment)

    return merged_segments

