"""
Graph-based algorithm for removing braided rivers and resolving flowline loops.

This module provides advanced algorithms using directed graph theory to handle
complex river networks with braided channels, loops, and parallel streams.
"""

import logging
from typing import List, Dict, Set, Tuple, Optional, DefaultDict
from collections import defaultdict, deque
import numpy as np
from pyflowline.classes.flowline import pyflowline
from pyflowline.classes.vertex import pyvertex

# Set up logger
logger = logging.getLogger(__name__)


class RiverNetworkGraph:
    """
    Directed graph representation of a river network for braided channel analysis.

    This class models the flowline network as a directed graph where:
    - Nodes represent vertices (confluence/divergence points)
    - Edges represent flowlines connecting vertices
    - Multiple edges between same nodes represent braided channels
    """

    def __init__(self, flowlines: List[pyflowline]):
        """
        Initialize the river network graph from flowlines.

        Args:
            flowlines: List of flowline objects representing the river network
        """
        self.flowlines = flowlines
        self.vertex_to_id: Dict[pyvertex, int] = {}
        self.id_to_vertex: Dict[int, pyvertex] = {}
        self.adjacency_list: DefaultDict[int, List[Tuple[int, int]]] = defaultdict(list)
        self.flowline_edges: Dict[int, Tuple[int, int]] = {}
        self.in_degree: DefaultDict[int, int] = defaultdict(int)
        self.out_degree: DefaultDict[int, int] = defaultdict(int)

        self._build_graph()

    def _build_graph(self):
        """Build the directed graph from flowlines."""
        vertex_id_counter = 0

        # Create vertex-to-ID mapping
        for i, flowline in enumerate(self.flowlines):
            start_vertex = flowline.pVertex_start
            end_vertex = flowline.pVertex_end

            # Assign IDs to vertices
            if start_vertex not in self.vertex_to_id:
                self.vertex_to_id[start_vertex] = vertex_id_counter
                self.id_to_vertex[vertex_id_counter] = start_vertex
                vertex_id_counter += 1

            if end_vertex not in self.vertex_to_id:
                self.vertex_to_id[end_vertex] = vertex_id_counter
                self.id_to_vertex[vertex_id_counter] = end_vertex
                vertex_id_counter += 1

            start_id = self.vertex_to_id[start_vertex]
            end_id = self.vertex_to_id[end_vertex]

            # Add edge to adjacency list
            self.adjacency_list[start_id].append((end_id, i))
            self.flowline_edges[i] = (start_id, end_id)

            # Update degree counts
            self.out_degree[start_id] += 1
            self.in_degree[end_id] += 1

    def find_braided_channels(self) -> Dict[Tuple[int, int], List[int]]:
        """
        Find braided channels (multiple flowlines between same vertex pair).

        Returns:
            Dictionary mapping (start_vertex_id, end_vertex_id) to list of flowline indices
        """
        channel_groups: DefaultDict[Tuple[int, int], List[int]] = defaultdict(list)

        for flowline_idx, (start_id, end_id) in self.flowline_edges.items():
            channel_groups[(start_id, end_id)].append(flowline_idx)

        # Return only groups with multiple channels (braided)
        return {key: indices for key, indices in channel_groups.items() if len(indices) > 1}

    def find_parallel_paths(self) -> List[List[int]]:
        """
        Find parallel paths between vertices (alternative routes).

        Returns:
            List of parallel path groups, each containing flowline indices
        """
        # Use DFS to find all paths between vertex pairs
        parallel_groups = []
        visited_pairs = set()

        for start_id in self.adjacency_list:
            for end_id, _ in self.adjacency_list[start_id]:
                if (start_id, end_id) in visited_pairs:
                    continue

                # Find all paths from start_id to end_id
                paths = self._find_all_paths(start_id, end_id, max_depth=5)
                if len(paths) > 1:
                    # Convert paths to flowline indices
                    flowline_groups = []
                    for path in paths:
                        path_flowlines = self._path_to_flowlines(path)
                        if path_flowlines:
                            flowline_groups.extend(path_flowlines)

                    if len(set(flowline_groups)) > 1:
                        parallel_groups.append(list(set(flowline_groups)))

                visited_pairs.add((start_id, end_id))

        return parallel_groups

    def _find_all_paths(self, start: int, end: int, max_depth: int = 5) -> List[List[int]]:
        """Find all paths between two vertices with limited depth."""
        paths = []

        def dfs(current: int, target: int, path: List[int], depth: int):
            if depth > max_depth:
                return
            if current == target:
                paths.append(path.copy())
                return

            for neighbor, _ in self.adjacency_list[current]:
                if neighbor not in path:  # Avoid cycles
                    path.append(neighbor)
                    dfs(neighbor, target, path, depth + 1)
                    path.pop()

        dfs(start, end, [start], 0)
        return paths

    def _path_to_flowlines(self, path: List[int]) -> List[int]:
        """Convert a vertex path to flowline indices."""
        flowline_indices = []
        for i in range(len(path) - 1):
            start_id, end_id = path[i], path[i + 1]
            for neighbor_id, flowline_idx in self.adjacency_list[start_id]:
                if neighbor_id == end_id:
                    flowline_indices.append(flowline_idx)
                    break
        return flowline_indices

    def detect_cycles(self) -> List[List[int]]:
        """
        Detect cycles in the network using DFS with recursion stack.

        Returns:
            List of cycles, where each cycle is a list of vertex IDs
        """
        visited = set()
        rec_stack = set()
        cycles = []

        def dfs_cycle_detection(node_id: int, path: List[int]) -> bool:
            try:
                visited.add(node_id)
                rec_stack.add(node_id)

                for neighbor_id, _ in self.adjacency_list[node_id]:
                    try:
                        if neighbor_id in rec_stack:
                            # Found a cycle - back edge to a node in recursion stack
                            try:
                                # Try to find where the cycle starts in the current path
                                cycle_start_idx = path.index(neighbor_id)
                                # Extract the cycle from the path
                                cycle = path[cycle_start_idx:] + [neighbor_id]
                                cycles.append(cycle)
                                logger.debug(f"Detected cycle: {cycle}")
                            except ValueError:
                                # neighbor_id not in current path - this means we have a back edge
                                # to a node that's in the recursion stack but not in current path
                                # This can happen with complex graph structures
                                cycle = [node_id, neighbor_id]
                                cycles.append(cycle)
                                logger.debug(f"Detected simple back-edge cycle: {cycle}")
                            except Exception as e:
                                logger.warning(f"Error processing cycle from {node_id} to {neighbor_id}: {e}")
                                # Fallback: create minimal cycle representation
                                cycles.append([node_id, neighbor_id])
                            return True
                        elif neighbor_id not in visited:
                            try:
                                new_path = path + [neighbor_id]
                                if dfs_cycle_detection(neighbor_id, new_path):
                                    return True
                            except RecursionError:
                                logger.error(f"Recursion limit reached during cycle detection at node {neighbor_id}")
                                # Add current path as potential cycle to avoid infinite recursion
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
                # Ensure cleanup even if error occurs
                try:
                    rec_stack.discard(node_id)  # Use discard to avoid KeyError
                except Exception:
                    pass
                return False

        try:
            # Create a list of node IDs to avoid "dictionary changed size during iteration" error
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

        logger.info(f"Cycle detection completed. Found {len(cycles)} cycles.")
        return cycles

    def get_sources(self) -> List[int]:
        """Get source nodes (headwaters) with no incoming edges."""
        return [node_id for node_id in self.id_to_vertex.keys() if self.in_degree[node_id] == 0]

    def get_sinks(self) -> List[int]:
        """Get sink nodes (outlets) with no outgoing edges."""
        return [node_id for node_id in self.id_to_vertex.keys() if self.out_degree[node_id] == 0]


def remove_braided_rivers(flowlines_input: List[pyflowline], outlet_vertex: Optional[pyvertex] = None) -> List[pyflowline]:
    """
    Remove braided rivers and resolve loops using advanced graph algorithms.

    This function employs sophisticated graph theory techniques to handle complex
    river networks with braided channels, parallel paths, and loops. It works with
    minimal flowline data, requiring only start and end vertices.

    1. **Graph Construction**: Models the river network as a directed graph
    2. **Braided Channel Detection**: Identifies multiple flowlines between same vertices
    3. **Parallel Path Analysis**: Finds alternative routes through the network
    4. **Cycle Resolution**: Detects and breaks cycles to create a proper tree structure
    5. **Optimal Selection**: Chooses best flowlines based on available criteria
    6. **Topological Ordering**: Ensures proper upstream-to-downstream flow sequence

    Args:
        flowlines_input (List[pyflowline]): Input flowlines that may contain braided
                                          channels, loops, and parallel paths.
                                          Only pVertex_start and pVertex_end are required.
                                          Optional attributes (iStream_order, dDrainage_area,
                                          dLength) will be initialized with defaults if missing.
        outlet_vertex (Optional[pyvertex]): Known outlet location to help guide the algorithm.
                                          If provided, flowlines connected to this vertex will
                                          be prioritized during selection.

    Returns:
        List[pyflowline]: Cleaned flowlines with braided channels resolved, loops
                         removed, and proper flow ordering

    Raises:
        TypeError: If input is not a list or contains non-flowline objects
        ValueError: If input list is empty or flowlines lack required vertex attributes

    Example:
        >>> # Minimal flowlines with only vertices
        >>> flowlines = [flowline1, flowline2, flowline3]  # Only need pVertex_start/end
        >>> clean_network = remove_braided_rivers(flowlines)
        >>> # Returns simplified network with optimal flowlines selected

    Algorithm Details:
        - **Required Attributes**: Only pVertex_start and pVertex_end are mandatory
        - **Optional Attributes**: iStream_order, dDrainage_area, dLength are auto-initialized
        - **Braided Resolution**: Selects flowline with highest available criteria
        - **Cycle Breaking**: Removes lowest priority flowlines from detected cycles
        - **Path Optimization**: For parallel paths, selects the most significant route
        - **Topological Sort**: Uses Kahn's algorithm for proper flow ordering
        - **Outlet Support**: Can work with outlet location information when available
    """
    # Input validation
    if not isinstance(flowlines_input, list):
        raise TypeError(f"Expected list of flowlines, got {type(flowlines_input)}")

    if not flowlines_input:
        raise ValueError("Input flowline list cannot be empty")

    # Validate flowline objects and required attributes
    # Only pVertex_start and pVertex_end are strictly required
    # Other attributes (iStream_order, dDrainage_area, dLength) have defaults in pyflowline class
    required_attrs = ['pVertex_start', 'pVertex_end']

    for i, flowline in enumerate(flowlines_input):
        if not isinstance(flowline, pyflowline):
            raise TypeError(f"Item at index {i} is not a pyflowline object: {type(flowline)}")

        # Check required attributes
        for attr in required_attrs:
            if not hasattr(flowline, attr):
                raise ValueError(f"Flowline at index {i} missing required attribute: {attr}")

    # If outlet vertex is provided, mark flowlines connected to it
    outlet_connected_flowlines = set()
    if outlet_vertex is not None:
        logger.info(f"Using outlet vertex information for enhanced selection")
        for i, flowline in enumerate(flowlines_input):
            if (flowline.pVertex_start == outlet_vertex or
                flowline.pVertex_end == outlet_vertex):
                outlet_connected_flowlines.add(i)
                # Boost priority for outlet-connected flowlines
                if hasattr(flowline, 'iStream_order'):
                    flowline.iStream_order = max(flowline.iStream_order, 2)
        logger.info(f"Found {len(outlet_connected_flowlines)} flowlines connected to outlet")

    logger.info(f"Starting braided river removal for {len(flowlines_input)} flowlines")

    # Step 1: Build directed graph representation
    network_graph = RiverNetworkGraph(flowlines_input)
    logger.info(f"Built network graph with {len(network_graph.vertex_to_id)} vertices")

    # Step 2: Identify and resolve braided channels
    braided_channels = network_graph.find_braided_channels()
    logger.info(f"Found {len(braided_channels)} braided channel groups")

    selected_flowlines = []
    processed_indices = set()

    # Resolve braided channels by selecting optimal flowlines
    for (start_id, end_id), flowline_indices in braided_channels.items():
        if any(idx in processed_indices for idx in flowline_indices):
            continue

        candidates = [flowlines_input[i] for i in flowline_indices]
        best_flowline = select_optimal_flowline(candidates)

        selected_flowlines.append(best_flowline)
        processed_indices.update(flowline_indices)

        logger.debug(f"Selected flowline (order={best_flowline.iStream_order}, "
                    f"area={best_flowline.dDrainage_area:.2e}) from {len(candidates)} braided options")

    # Step 3: Add non-braided flowlines
    for i, flowline in enumerate(flowlines_input):
        if i not in processed_indices:
            selected_flowlines.append(flowline)

    logger.info(f"After braided channel resolution: {len(selected_flowlines)} flowlines")

    # Step 4: Detect and resolve parallel paths
    if len(selected_flowlines) > 1:
        selected_flowlines = resolve_parallel_paths(selected_flowlines)
        logger.info(f"After parallel path resolution: {len(selected_flowlines)} flowlines")

    # Step 5: Detect and break cycles
    if len(selected_flowlines) > 1:
        selected_flowlines = break_network_cycles(selected_flowlines)
        logger.info(f"After cycle resolution: {len(selected_flowlines)} flowlines")

    # Step 6: Topological ordering for proper flow sequence
    if selected_flowlines:
        selected_flowlines = topological_sort_flowlines(selected_flowlines)

    # Step 7: Assign sequential indices
    for i, flowline in enumerate(selected_flowlines):
        flowline.lFlowlineIndex = i

    logger.info(f"Braided river removal complete: {len(selected_flowlines)} flowlines "
               f"selected from {len(flowlines_input)} input flowlines")

    return selected_flowlines


def select_optimal_flowline(candidates: List[pyflowline]) -> pyflowline:
    """
    Select the optimal flowline from candidates using hydrological criteria.

    Selection priority:
    1. Highest stream order (Strahler order)
    2. Largest drainage area
    3. Longest flowline length
    4. Lowest sinuosity (more direct path)
    5. First in list (deterministic fallback)

    Args:
        candidates: List of candidate flowlines to choose from

    Returns:
        The optimal flowline based on selection criteria
    """
    if len(candidates) == 1:
        return candidates[0]

    def priority_score(flowline: pyflowline) -> Tuple:
        """Calculate priority score for flowline selection."""
        # Only use optional attributes if they have valid values (> 0)
        stream_order = getattr(flowline, 'iStream_order', -1)
        stream_order = stream_order if stream_order > 0 else 1  # Use 1 as neutral default

        drainage_area = getattr(flowline, 'dDrainage_area', 0.0)
        drainage_area = drainage_area if drainage_area > 0 else 0.0  # Use 0 as neutral default

        length = getattr(flowline, 'dLength', 0.0)
        length = length if length > 0 else 0.0  # Use 0 as neutral default

        # Calculate sinuosity if not already computed, with fallback
        try:
            sinuosity = flowline.dSinuosity if hasattr(flowline, 'dSinuosity') and flowline.dSinuosity > 0 else flowline.get_sinuosity()
        except (AttributeError, ValueError):
            # Fallback to 1.0 (straight line) if sinuosity cannot be calculated
            sinuosity = 1.0

        return (
            stream_order,                     # Higher is better
            drainage_area,                    # Larger is better
            length,                           # Longer is better
            -sinuosity                        # Lower sinuosity is better (more direct)
        )

    return max(candidates, key=priority_score)


def resolve_parallel_paths(flowlines: List[pyflowline]) -> List[pyflowline]:
    """
    Resolve parallel paths by selecting the most significant route.

    Identifies alternative routes between distant vertices and selects
    the path with the highest cumulative hydrological significance.

    Args:
        flowlines: List of flowlines that may contain parallel paths

    Returns:
        Flowlines with parallel paths resolved
    """
    if len(flowlines) <= 1:
        return flowlines

    network_graph = RiverNetworkGraph(flowlines)
    parallel_groups = network_graph.find_parallel_paths()

    if not parallel_groups:
        return flowlines

    logger.info(f"Found {len(parallel_groups)} parallel path groups")

    flowlines_to_remove = set()

    for group in parallel_groups:
        if len(group) <= 1:
            continue

        # Calculate path significance for each flowline in the group
        path_scores = []
        for flowline_idx in group:
            flowline = flowlines[flowline_idx]
            # Only use attributes if they have valid values (> 0)
            stream_order = getattr(flowline, 'iStream_order', -1)
            stream_order = stream_order if stream_order > 0 else 0

            drainage_area = getattr(flowline, 'dDrainage_area', 0.0)
            drainage_area = drainage_area if drainage_area > 0 else 0.0

            length = getattr(flowline, 'dLength', 0.0)
            length = length if length > 0 else 0.0

            # Calculate score based on available valid attributes
            if stream_order > 0 and drainage_area > 0:
                score = stream_order * drainage_area
            elif stream_order > 0:
                score = stream_order * 1000  # Boost stream order when drainage area unavailable
            elif drainage_area > 0:
                score = drainage_area
            else:
                score = length  # Fallback to length if no hydrological data available

            path_scores.append((score, flowline_idx, flowline))

        # Keep the highest scoring flowline, remove others
        path_scores.sort(reverse=True)
        best_score, best_idx, best_flowline = path_scores[0]

        for score, idx, flowline in path_scores[1:]:
            flowlines_to_remove.add(flowline)

        logger.debug(f"Selected path with score {best_score:.2e} from {len(group)} parallel options")

    # Return flowlines with parallel paths removed
    result = [f for f in flowlines if f not in flowlines_to_remove]
    logger.info(f"Removed {len(flowlines_to_remove)} flowlines to resolve parallel paths")

    return result


def break_network_cycles(flowlines: List[pyflowline]) -> List[pyflowline]:
    """
    Detect and break cycles in the network by removing lowest priority flowlines.

    Uses depth-first search to detect cycles and removes the flowline with
    the lowest hydrological significance from each cycle.

    Args:
        flowlines: List of flowlines that may contain cycles

    Returns:
        Flowlines with cycles broken (acyclic network)
    """
    if len(flowlines) <= 1:
        return flowlines

    network_graph = RiverNetworkGraph(flowlines)
    cycles = network_graph.detect_cycles()

    if not cycles:
        return flowlines

    logger.info(f"Detected {len(cycles)} cycles in network")

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
            for flowline_idx, (s_id, e_id) in network_graph.flowline_edges.items():
                if s_id == start_vertex_id and e_id == end_vertex_id:
                    cycle_flowlines.append(flowlines[flowline_idx])

        if cycle_flowlines:
            # Remove the lowest priority flowline from the cycle
            def cycle_priority(flowline: pyflowline) -> Tuple:
                # Only use attributes if they have valid values (> 0)
                stream_order = getattr(flowline, 'iStream_order', -1)
                stream_order = stream_order if stream_order > 0 else 1

                drainage_area = getattr(flowline, 'dDrainage_area', 0.0)
                drainage_area = drainage_area if drainage_area > 0 else 0.0

                length = getattr(flowline, 'dLength', 0.0)
                length = length if length > 0 else 0.0

                return (stream_order, drainage_area, length)

            worst_flowline = min(cycle_flowlines, key=cycle_priority)
            flowlines_to_remove.add(worst_flowline)

            logger.debug(f"Removing flowline (order={worst_flowline.iStream_order}) to break cycle")

    # Return flowlines with cycle-causing ones removed
    result = [f for f in flowlines if f not in flowlines_to_remove]
    logger.info(f"Removed {len(flowlines_to_remove)} flowlines to break cycles")

    return result


def topological_sort_flowlines(flowlines: List[pyflowline]) -> List[pyflowline]:
    """
    Perform topological sort to order flowlines from upstream to downstream.

    Uses Kahn's algorithm to ensure proper flow ordering in the network.
    Headwater flowlines appear first, outlet flowlines appear last.

    Args:
        flowlines: List of flowlines to sort topologically

    Returns:
        Topologically sorted flowlines (upstream to downstream order)
    """
    if len(flowlines) <= 1:
        return flowlines

    network_graph = RiverNetworkGraph(flowlines)

    # Kahn's algorithm for topological sorting
    in_degree = network_graph.in_degree.copy()
    queue = deque(network_graph.get_sources())
    sorted_nodes = []

    while queue:
        current_node = queue.popleft()
        sorted_nodes.append(current_node)

        # Process all outgoing edges
        for neighbor_id, flowline_idx in network_graph.adjacency_list[current_node]:
            in_degree[neighbor_id] -= 1
            if in_degree[neighbor_id] == 0:
                queue.append(neighbor_id)

    # Map sorted nodes back to flowlines
    node_to_flowlines = defaultdict(list)
    for flowline in flowlines:
        start_id = network_graph.vertex_to_id[flowline.pVertex_start]
        node_to_flowlines[start_id].append(flowline)

    # Build result maintaining topological order
    result = []
    for node_id in sorted_nodes:
        # Sort flowlines at each node by priority for deterministic ordering
        def sort_key(f):
            stream_order = getattr(f, 'iStream_order', -1)
            stream_order = stream_order if stream_order > 0 else 1

            drainage_area = getattr(f, 'dDrainage_area', 0.0)
            drainage_area = drainage_area if drainage_area > 0 else 0.0

            return (stream_order, drainage_area)

        node_flowlines = sorted(node_to_flowlines[node_id], key=sort_key, reverse=True)
        result.extend(node_flowlines)

    # Add any remaining flowlines (isolated components)
    added_flowlines = set(result)
    remaining = [f for f in flowlines if f not in added_flowlines]
    if remaining:
        logger.warning(f"Found {len(remaining)} isolated flowlines, adding to end")
        result.extend(remaining)

    logger.info(f"Topological sort complete: {len(result)} flowlines ordered")
    return result