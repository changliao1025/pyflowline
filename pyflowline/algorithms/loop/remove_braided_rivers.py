import logging
import os
from typing import List, Dict, Set, Tuple, Optional, DefaultDict
from collections import defaultdict, deque
import numpy as np
from pyflowline.classes.flowline import pyflowline
from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.rivergraph import pyrivergraph
from pyflowline.formats.export_flowline import export_flowline_to_geojson

# Set up logger
logger = logging.getLogger(__name__)

def _write_debug_flowlines(flowlines: List[pyflowline], debug_folder: str, filename: str, step_description: str = ""):
    """
    Write flowlines to debug folder as GeoJSON for debugging purposes.

    Args:
        flowlines: List of flowlines to export
        debug_folder: Path to debug output folder
        filename: Name of the output file (without extension)
        step_description: Description of the processing step for logging
    """
    if not flowlines:
        logger.warning(f"No flowlines to write for debug step: {step_description}")
        return

    try:
        # Ensure debug folder exists
        os.makedirs(debug_folder, exist_ok=True)

        # Create full file path
        output_path = os.path.join(debug_folder, f"{filename}.geojson")

        # Export flowlines to GeoJSON
        export_flowline_to_geojson(flowlines, output_path)

        logger.info(f"Debug output: {step_description} - {len(flowlines)} flowlines written to {output_path}")

    except Exception as e:
        logger.error(f"Failed to write debug output for {step_description}: {e}")


def remove_braided_rivers(flowlines_input: List[pyflowline], outlet_vertex: Optional[pyvertex] = None, sFolder_out: Optional[str] = None, iFlag_merge_segments: bool = True) -> List[pyflowline]:
    """
    Remove braided rivers and resolve loops using advanced graph algorithms.

    This function employs sophisticated graph theory techniques to handle complex
    river networks with braided channels, parallel paths, and loops. It works with
    minimal flowline data, requiring only start and end vertices.

    1. **Graph Construction**: Models the river network as a directed graph
    2. **Linear Segment Merging**: Connects artificially segmented flowlines (optional)
    3. **Braided Channel Detection**: Identifies multiple flowlines between same vertices
    4. **Parallel Path Analysis**: Finds alternative routes through the network
    5. **Cycle Resolution**: Detects and breaks cycles to create a proper tree structure
    6. **Optimal Selection**: Chooses best flowlines based on available criteria
    7. **Topological Ordering**: Ensures proper upstream-to-downstream flow sequence

    Args:
        flowlines_input (List[pyflowline]): Input flowlines that may contain braided
                                          channels, loops, and parallel paths.
                                          Only pVertex_start and pVertex_end are required.
                                          Optional attributes (iStream_order, dDrainage_area,
                                          dLength) will be initialized with defaults if missing.
        outlet_vertex (Optional[pyvertex]): Known outlet location to help guide the algorithm.
                                          If provided, flowlines connected to this vertex will
                                          be prioritized during selection.
        sFolder_out (Optional[str]): Path to folder for debug output files. If provided,
                                    intermediate processing steps will be saved as GeoJSON files
                                    for debugging and analysis purposes.
        iFlag_merge_segments (bool): Whether to merge linear segments before processing.
                                   Default True. Set to False to skip segment merging.

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

        >>> # With debug output and segment merging
        >>> clean_network = remove_braided_rivers(flowlines, sFolder_out="./debug_output")
        >>> # Creates debug files showing each processing step

        >>> # Skip segment merging if already done
        >>> clean_network = remove_braided_rivers(flowlines, iFlag_merge_segments=False)

    Algorithm Details:
        - **Required Attributes**: Only pVertex_start and pVertex_end are mandatory
        - **Optional Attributes**: iStream_order, dDrainage_area, dLength are auto-initialized
        - **Linear Segment Merging**: Connects degree-2 vertices to form continuous reaches
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

    # Debug output: Initial input flowlines
    if sFolder_out:
        _write_debug_flowlines(flowlines_input, sFolder_out, "01_input_flowlines",
                              "Initial input flowlines before processing")

    # Step 1: Build directed graph representation
    network_graph = pyrivergraph(flowlines_input)
    logger.info(f"Built network graph with {len(network_graph.vertex_to_id)} vertices")

    # Step 1.5: Merge linear segments if requested
    current_flowlines = flowlines_input
    if iFlag_merge_segments:
        logger.info("Merging linear segments before braided channel processing")
        current_flowlines = network_graph.merge_linear_segments(flowlines_input)
        logger.info(f"After linear segment merging: {len(current_flowlines)} flowlines")

        # Rebuild graph with merged flowlines for subsequent processing
        if len(current_flowlines) != len(flowlines_input):
            network_graph = pyrivergraph(current_flowlines)
            logger.info(f"Rebuilt network graph with {len(network_graph.vertex_to_id)} vertices after merging")

    # Step 2: Identify and resolve braided channels
    braided_channels = network_graph.find_braided_channels()
    logger.info(f"Found {len(braided_channels)} braided channel groups")

    selected_flowlines = []
    processed_indices = set()

    # Resolve braided channels by selecting optimal flowlines
    for (start_id, end_id), flowline_indices in braided_channels.items():
        if any(idx in processed_indices for idx in flowline_indices):
            continue

        candidates = [current_flowlines[i] for i in flowline_indices]
        best_flowline = select_optimal_flowline(candidates)

        selected_flowlines.append(best_flowline)
        processed_indices.update(flowline_indices)

        logger.debug(f"Selected flowline (order={best_flowline.iStream_order}, "
                    f"area={best_flowline.dDrainage_area:.2e}) from {len(candidates)} braided options")

    # Step 3: Add non-braided flowlines
    for i, flowline in enumerate(current_flowlines):
        if i not in processed_indices:
            selected_flowlines.append(flowline)

    logger.info(f"After braided channel resolution: {len(selected_flowlines)} flowlines")

    # Debug output: After braided channel resolution
    if sFolder_out:
        _write_debug_flowlines(selected_flowlines, sFolder_out, "02_after_braided_resolution",
                              "Flowlines after braided channel resolution")

    # Step 4: Detect and resolve parallel paths
    if len(selected_flowlines) > 1:
        selected_flowlines = resolve_parallel_paths(selected_flowlines)
        logger.info(f"After parallel path resolution: {len(selected_flowlines)} flowlines")

        # Debug output: After parallel path resolution
        if sFolder_out:
            _write_debug_flowlines(selected_flowlines, sFolder_out, "03_after_parallel_resolution",
                                  "Flowlines after parallel path resolution")

    # Step 5: Detect and break cycles
    if len(selected_flowlines) > 1:
        selected_flowlines = break_network_cycles(selected_flowlines)
        logger.info(f"After cycle resolution: {len(selected_flowlines)} flowlines")

        # Debug output: After cycle resolution
        if sFolder_out:
            _write_debug_flowlines(selected_flowlines, sFolder_out, "04_after_cycle_resolution",
                                  "Flowlines after cycle resolution")

    # Step 6: Topological ordering for proper flow sequence
    if selected_flowlines:
        selected_flowlines = topological_sort_flowlines(selected_flowlines)

        # Debug output: After topological sorting
        if sFolder_out:
            _write_debug_flowlines(selected_flowlines, sFolder_out, "05_after_topological_sort",
                                  "Flowlines after topological sorting")

    # Step 7: Assign sequential indices
    for i, flowline in enumerate(selected_flowlines):
        flowline.lFlowlineIndex = i

    # Debug output: Final result
    if sFolder_out:
        _write_debug_flowlines(selected_flowlines, sFolder_out, "06_final_result",
                              "Final processed flowlines with indices assigned")

    logger.info(f"Braided river removal complete: {len(selected_flowlines)} flowlines "
                f"selected from {len(flowlines_input)} original input flowlines "
                f"({len(current_flowlines)} after linear segment merging)")

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
        logger.debug("Skipping parallel path resolution: insufficient flowlines")
        return flowlines

    try:
        network_graph = pyrivergraph(flowlines)
        parallel_groups = network_graph.find_parallel_paths()

        if not parallel_groups:
            logger.debug("No parallel paths found")
            return flowlines

        logger.info(f"Found {len(parallel_groups)} parallel path groups")
    except Exception as e:
        logger.error(f"Error building network graph or finding parallel paths: {e}")
        return flowlines

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
                if flowline_idx >= len(flowlines):
                    logger.warning(f"Invalid flowline index {flowline_idx}, skipping")
                    continue

                flowline = flowlines[flowline_idx]
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
        result = [flowlines[i] for i in range(len(flowlines)) if i not in flowlines_to_remove]
        logger.info(f"Removed {len(flowlines_to_remove)} flowlines to resolve parallel paths")

        # Validate result
        if len(result) == 0 and len(flowlines) > 0:
            logger.warning("All flowlines were removed during parallel path resolution, returning original")
            return flowlines

        return result
    except Exception as e:
        logger.error(f"Error filtering flowlines during parallel path resolution: {e}")
        return flowlines


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

    network_graph = pyrivergraph(flowlines)
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
            for flowline_idx, (s_id, e_id) in network_graph.aFlowline_edges.items():
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
    Perform topological sort to order flowlines from downstream to upstream.

    Uses reverse topological ordering to ensure proper flow ordering in the network.
    Outlet flowlines appear first, headwater flowlines appear last.

    Args:
        flowlines: List of flowlines to sort topologically

    Returns:
        Topologically sorted flowlines (downstream to upstream order)
    """
    if len(flowlines) <= 1:
        return flowlines

    network_graph = pyrivergraph(flowlines)

    # Reverse topological sorting - start from sinks (outlets) instead of sources
    out_degree = network_graph.out_degree.copy()
    queue = deque(network_graph.get_sinks())  # Start from outlets (sinks)
    sorted_nodes = []

    while queue:
        current_node = queue.popleft()
        sorted_nodes.append(current_node)

        # Process all incoming edges (reverse direction)
        # Create list to avoid dictionary modification issues during iteration
        start_ids = list(network_graph.adjacency_list.keys())
        for start_id in start_ids:
            neighbors = list(network_graph.adjacency_list[start_id])
            for neighbor_id, flowline_idx in neighbors:
                if neighbor_id == current_node:
                    out_degree[start_id] -= 1
                    if out_degree[start_id] == 0:
                        queue.append(start_id)

    # Map sorted nodes back to flowlines - use end vertex instead of start vertex
    node_to_flowlines = defaultdict(list)
    for flowline in flowlines:
        end_id = network_graph.vertex_to_id[flowline.pVertex_end]
        node_to_flowlines[end_id].append(flowline)

    # Build result maintaining reverse topological order (downstream first)
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

    logger.info(f"Reverse topological sort complete: {len(result)} flowlines ordered (downstream first)")
    return result