"""
Fast duplicate removal algorithms for pyflowline objects.

This module provides efficient methods to remove duplicate flowlines, vertices,
and other pyflowline objects using multiple strategies for optimal performance.
"""

import numpy as np
from typing import List, Set, Tuple, Union, Any, Dict
from collections import defaultdict
import hashlib

def remove_duplicate_flowlines_fast(aFlowline_in: List) -> List:
    """
    Efficiently remove duplicate flowlines based purely on geometry (vertices).

    This function ignores all IDs and focuses only on the geometric properties
    of flowlines by comparing all vertices in the flowline path.

    Args:
        aFlowline_in: List of pyflowline objects

    Returns:
        List of unique pyflowline objects (based on geometry only)

    Performance: O(n) average case for simple geometries

    Example:
        >>> flowlines = [flowline1, flowline2, flowline1]  # flowline1 is duplicate
        >>> unique_flowlines = remove_duplicate_flowlines_fast(flowlines)
        >>> len(unique_flowlines)  # Should be 2
    """
    if not aFlowline_in:
        return []

    unique_flowlines = []
    seen_geometries = set()

    for flowline in aFlowline_in:
        # Create geometric hash from ALL vertices in the flowline
        geom_hash = _create_complete_geometric_hash(flowline)

        if geom_hash not in seen_geometries:
            seen_geometries.add(geom_hash)
            unique_flowlines.append(flowline)

    return unique_flowlines

def remove_duplicate_flowlines_by_geometry(aFlowline_in: List, tolerance: float = 1e-10) -> List:
    """
    Remove duplicate flowlines based on geometric similarity with tolerance.

    Args:
        aFlowline_in: List of pyflowline objects
        tolerance: Coordinate tolerance for considering flowlines as duplicates

    Returns:
        List of unique pyflowline objects
    """
    if not aFlowline_in:
        return []

    unique_flowlines = []

    for flowline in aFlowline_in:
        is_duplicate = False

        for existing in unique_flowlines:
            if _are_geometrically_similar(flowline, existing, tolerance):
                is_duplicate = True
                break

        if not is_duplicate:
            unique_flowlines.append(flowline)

    return unique_flowlines

def remove_duplicate_vertices_fast(aVertex_in: List) -> List:
    """
    Efficiently remove duplicate vertices from a list.

    Args:
        aVertex_in: List of pyvertex objects

    Returns:
        List of unique pyvertex objects
    """
    if not aVertex_in:
        return []

    # Use the built-in __hash__ and __eq__ methods of pyvertex
    # which are based on coordinates with proper precision handling
    return list(set(aVertex_in))

def remove_duplicate_edges_fast(aEdge_in: List) -> List:
    """
    Efficiently remove duplicate edges from a list.

    Args:
        aEdge_in: List of pyedge objects

    Returns:
        List of unique pyedge objects
    """
    if not aEdge_in:
        return []

    # Use the built-in __hash__ and __eq__ methods of pyedge
    return list(set(aEdge_in))

def remove_duplicates_generic(objects: List,
                            id_attrs: List[str] = None,
                            coord_attrs: List[str] = None,
                            tolerance: float = 1e-10) -> List:
    """
    Generic function to remove duplicates from any pyflowline object list.

    Args:
        objects: List of objects to deduplicate
        id_attrs: List of ID attribute names to check (e.g., ['lFlowlineID', 'lFlowlineIndex'])
        coord_attrs: List of coordinate attribute names for geometric comparison
        tolerance: Tolerance for coordinate comparison

    Returns:
        List of unique objects
    """
    if not objects:
        return []

    # Try to use built-in hash/eq methods first (most efficient)
    try:
        return list(set(objects))
    except TypeError:
        # Fallback to manual deduplication if objects don't support hashing
        pass

    unique_objects = []
    seen_hashes = set()

    for obj in objects:
        obj_hash = _create_object_hash(obj, id_attrs, coord_attrs, tolerance)

        if obj_hash not in seen_hashes:
            seen_hashes.add(obj_hash)
            unique_objects.append(obj)

    return unique_objects

def _create_geometric_hash(flowline) -> str:
    """
    Create a hash based on flowline geometry (start and end vertices only).

    Args:
        flowline: pyflowline object

    Returns:
        str: Hash string representing the geometry
    """
    if not hasattr(flowline, 'pVertex_start') or not hasattr(flowline, 'pVertex_end'):
        return str(hash(flowline))

    start_coords = (
        round(flowline.pVertex_start.dLongitude_degree, 10),
        round(flowline.pVertex_start.dLatitude_degree, 10)
    )
    end_coords = (
        round(flowline.pVertex_end.dLongitude_degree, 10),
        round(flowline.pVertex_end.dLatitude_degree, 10)
    )

    # Create hash from coordinate tuple
    coord_str = f"{start_coords[0]},{start_coords[1]},{end_coords[0]},{end_coords[1]}"
    return hashlib.md5(coord_str.encode()).hexdigest()

def _create_complete_geometric_hash(flowline) -> str:
    """
    Create a hash based on ALL vertices in the flowline geometry.

    This function creates a comprehensive hash by including all vertices
    in the flowline path, not just start and end points. Direction matters -
    flowlines with opposite directions are considered different.

    Args:
        flowline: pyflowline object

    Returns:
        str: Hash string representing the complete geometry with direction
    """
    if not hasattr(flowline, 'aVertex') or not flowline.aVertex:
        # Fallback to start/end vertices if aVertex is not available
        return _create_geometric_hash(flowline)

    # Extract all vertex coordinates in order (direction matters)
    coord_list = []
    for vertex in flowline.aVertex:
        if hasattr(vertex, 'dLongitude_degree') and hasattr(vertex, 'dLatitude_degree'):
            coord_list.extend([
                round(vertex.dLongitude_degree, 10),
                round(vertex.dLatitude_degree, 10)
            ])

    # Create hash from coordinates in their original order
    # Direction matters - reversed flowlines will have different hashes
    coord_str = ','.join(map(str, coord_list))

    return hashlib.md5(coord_str.encode()).hexdigest()

def _are_geometrically_similar(flowline1, flowline2, tolerance: float) -> bool:
    """
    Check if two flowlines are geometrically similar within tolerance.

    This function compares ALL vertices in both flowlines in the same direction.
    Direction matters - flowlines with opposite directions are considered different.

    Args:
        flowline1: First pyflowline object
        flowline2: Second pyflowline object
        tolerance: Coordinate tolerance

    Returns:
        bool: True if flowlines are geometrically similar in the same direction
    """
    # Try to use all vertices first
    if (hasattr(flowline1, 'aVertex') and hasattr(flowline2, 'aVertex') and
        flowline1.aVertex and flowline2.aVertex):

        vertices1 = flowline1.aVertex
        vertices2 = flowline2.aVertex

        # Check if they have the same number of vertices
        if len(vertices1) != len(vertices2):
            return False

        # Check forward direction only (direction matters)
        for v1, v2 in zip(vertices1, vertices2):
            if (abs(v1.dLongitude_degree - v2.dLongitude_degree) >= tolerance or
                abs(v1.dLatitude_degree - v2.dLatitude_degree) >= tolerance):
                return False

        return True

    # Fallback to start/end vertex comparison if aVertex is not available
    if not all(hasattr(fl, 'pVertex_start') and hasattr(fl, 'pVertex_end')
               for fl in [flowline1, flowline2]):
        return False

    # Check if start and end vertices are within tolerance (same direction only)
    start_similar = (
        abs(flowline1.pVertex_start.dLongitude_degree - flowline2.pVertex_start.dLongitude_degree) < tolerance and
        abs(flowline1.pVertex_start.dLatitude_degree - flowline2.pVertex_start.dLatitude_degree) < tolerance
    )

    end_similar = (
        abs(flowline1.pVertex_end.dLongitude_degree - flowline2.pVertex_end.dLongitude_degree) < tolerance and
        abs(flowline1.pVertex_end.dLatitude_degree - flowline2.pVertex_end.dLatitude_degree) < tolerance
    )

    return start_similar and end_similar

def _create_object_hash(obj, id_attrs: List[str] = None, coord_attrs: List[str] = None,
                       tolerance: float = 1e-10) -> str:
    """
    Create a hash for any object based on specified attributes.

    Args:
        obj: Object to hash
        id_attrs: List of ID attribute names
        coord_attrs: List of coordinate attribute names
        tolerance: Tolerance for coordinate rounding

    Returns:
        str: Hash string for the object
    """
    hash_components = []

    # Add ID attributes to hash
    if id_attrs:
        for attr in id_attrs:
            if hasattr(obj, attr):
                val = getattr(obj, attr)
                if val is not None and val >= 0:  # Valid ID/index
                    hash_components.append(f"{attr}:{val}")

    # Add coordinate attributes to hash
    if coord_attrs:
        for attr in coord_attrs:
            if hasattr(obj, attr):
                val = getattr(obj, attr)
                if val is not None:
                    if isinstance(val, (list, np.ndarray)):
                        # Handle coordinate arrays
                        rounded_coords = [round(float(x), int(-np.log10(tolerance))) for x in val]
                        hash_components.append(f"{attr}:{rounded_coords}")
                    else:
                        # Handle single coordinate values
                        rounded_val = round(float(val), int(-np.log10(tolerance)))
                        hash_components.append(f"{attr}:{rounded_val}")

    # Fallback: use string representation if no specific attributes
    if not hash_components:
        hash_components.append(str(obj))

    # Create hash from components
    hash_string = '|'.join(hash_components)
    return hashlib.md5(hash_string.encode()).hexdigest()

def remove_duplicates_by_spatial_clustering(objects: List,
                                          coord_attrs: List[str],
                                          cluster_tolerance: float = 1e-6) -> List:
    """
    Remove objects that are spatially clustered together.

    This is useful for removing near-duplicate coordinates that might not be
    exactly identical due to floating point precision issues.

    Args:
        objects: List of objects with spatial coordinates
        coord_attrs: List of coordinate attribute names (e.g., ['dLongitude_degree', 'dLatitude_degree'])
        cluster_tolerance: Spatial tolerance for clustering

    Returns:
        List of objects with spatial duplicates removed
    """
    if not objects or not coord_attrs:
        return objects

    unique_objects = []

    for obj in objects:
        is_duplicate = False

        # Extract coordinates for this object
        obj_coords = []
        for attr in coord_attrs:
            if hasattr(obj, attr):
                obj_coords.append(getattr(obj, attr))
            else:
                obj_coords.append(0.0)  # Default if attribute missing

        # Check against existing unique objects
        for existing in unique_objects:
            existing_coords = []
            for attr in coord_attrs:
                if hasattr(existing, attr):
                    existing_coords.append(getattr(existing, attr))
                else:
                    existing_coords.append(0.0)

            # Calculate spatial distance
            distance = np.sqrt(sum((a - b) ** 2 for a, b in zip(obj_coords, existing_coords)))

            if distance < cluster_tolerance:
                is_duplicate = True
                break

        if not is_duplicate:
            unique_objects.append(obj)

    return unique_objects

# Convenience function that matches the original API
def remove_duplicate_flowline(aFlowline_in: List) -> List:
    """
    Remove duplicate flowlines - improved version of the original function.

    This function maintains compatibility with the original API while providing
    much better performance and accuracy.

    Args:
        aFlowline_in: List of pyflowline objects

    Returns:
        List of unique pyflowline objects
    """
    return remove_duplicate_flowlines_fast(aFlowline_in)