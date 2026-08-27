"""
PySquare class for representing square mesh cells in flowline networks.

This module defines the pysquare class which extends pymeshcell to provide
square-specific mesh cell functionality for flowline network modeling.
"""

import json
from json import JSONEncoder
from typing import List, Optional, Any, Tuple
import numpy as np

from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area
from pyearthriver.classes.vertex import pyvertex
from pyearthriver.classes.edge import pyedge
from pyearthriver.classes.flowline import pyflowline
from pyflowline.classes.meshcell import pymeshcell


class SquareClassEncoder(JSONEncoder):
    """
    Custom JSON encoder for pysquare objects.

    Handles numpy data types, pyvertex objects, and other complex types,
    converting them to native Python types for JSON serialization.
    """

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, (np.float32, np.float64)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, list):
            return obj
        if isinstance(obj, pyvertex):
            return json.loads(obj.tojson())
        if pyedge and isinstance(obj, pyedge):
            return getattr(obj, "lEdgeID", str(obj))
        if pyflowline and isinstance(obj, pyflowline):
            return getattr(obj, "lFlowlineID", str(obj))
        if isinstance(obj, pysquare):
            return obj.lCellID
        return JSONEncoder.default(self, obj)


class pysquare(pymeshcell):
    """
    Square mesh cell class for flowline network representation.

    Extends the pymeshcell class to provide square-specific functionality for
    mesh cells in flowline topology and stream network analysis. A square cell
    represents a rectangular polygonal area with 4 vertices and 4 edges.

    This class is designed for regular grid-based mesh representations where
    cells have uniform square geometry, commonly used in structured grid
    hydrological models.

    Attributes:
        All attributes from pymeshcell plus square-specific properties.
        The square cell is defined by exactly 4 edges and 4 vertices forming
        a rectangular boundary.

    Args:
        dLon (float): The longitude of the square center in degrees
        dLat (float): The latitude of the square center in degrees
        aEdge (List): A list of exactly 4 edges that define the square boundary
        aVertex (List): A list of exactly 4 vertices that define the square boundary

    Raises:
        ValueError: If the cell doesn't have exactly 4 edges and 4 vertices
        TypeError: If input parameters are not of the expected types

    Example:
        >>> center_lon, center_lat = -77.0, 38.0
        >>> edges = [edge1, edge2, edge3, edge4]  # List of 4 pyedge objects
        >>> vertices = [vertex1, vertex2, vertex3, vertex4]  # List of 4 vertices
        >>> square = pysquare(center_lon, center_lat, edges, vertices)
        >>> square.set_cell_id(100)
        >>> print(square)
        pysquare(ID=100, Center=(-77.0, 38.0), Area=1234.56m²)
    """

    def __init__(self, dLon: float, dLat: float, aEdge: List, aVertex: List) -> None:
        """
        Initialize a square mesh cell object.

        Args:
            dLon (float): The longitude of the square center in degrees
            dLat (float): The latitude of the square center in degrees
            aEdge (List): A list of exactly 4 edges that define the square boundary
            aVertex (List): A list of exactly 4 vertices that define the square boundary

        Raises:
            ValueError: If the cell doesn't have exactly 4 edges and vertices
            TypeError: If input parameters are not of the expected types
        """
        # Input validation
        if not isinstance(dLon, (int, float, np.number)):
            raise TypeError(f"dLon must be a number, got {type(dLon)}")
        if not isinstance(dLat, (int, float, np.number)):
            raise TypeError(f"dLat must be a number, got {type(dLat)}")
        if not isinstance(aEdge, list):
            raise TypeError(f"aEdge must be a list, got {type(aEdge)}")
        if not isinstance(aVertex, list):
            raise TypeError(f"aVertex must be a list, got {type(aVertex)}")

        # Validate coordinate ranges
        if not (-180 <= dLon <= 180):
            raise ValueError(f"Longitude must be between -180 and 180, got {dLon}")
        if not (-90 <= dLat <= 90):
            raise ValueError(f"Latitude must be between -90 and 90, got {dLat}")

        # Validate square geometry requirements
        nEdge = len(aEdge)
        nVertex = len(aVertex)

        if nEdge != 4:
            raise ValueError(f"Square cell must have exactly 4 edges, got {nEdge}")
        if nVertex != 4:
            raise ValueError(f"Square cell must have exactly 4 vertices, got {nVertex}")

        # Validate edges and vertices are not None
        for i, edge in enumerate(aEdge):
            if edge is None:
                raise ValueError(f"Edge at index {i} cannot be None")

        for i, vertex in enumerate(aVertex):
            if vertex is None:
                raise ValueError(f"Vertex at index {i} cannot be None")

        # Initialize parent class
        super().__init__(dLon, dLat, aEdge, aVertex)

        # Square-specific initialization
        self.nEdge = 4
        self.nVertex = 4  # Note: Different from nPoint which includes closing vertex

        # Calculate initial properties
        self.calculate_cell_area()
        self.calculate_edge_length()

    def __repr__(self) -> str:
        """
        Return a detailed string representation of the square cell.

        Returns:
            str: Detailed representation including cell ID, center coordinates, and area
        """
        return (
            f"pysquare(ID={self.lCellID}, "
            f"Center=({self.dLongitude_center_degree:.6f}, {self.dLatitude_center_degree:.6f}), "
            f"Area={self.dArea:.2f}m², Side={self.dLength:.2f}m, "
            f"Flowlines={self.nFlowline}, Neighbors={self.nNeighbor})"
        )

    def __str__(self) -> str:
        """
        Return a concise string representation of the square cell.

        Returns:
            str: Concise representation with ID, center coordinates, and area
        """
        return (
            f"pysquare(ID={self.lCellID}, "
            f"Center=({self.dLongitude_center_degree:.6f}, {self.dLatitude_center_degree:.6f}), "
            f"Area={self.dArea:.2f}m²)"
        )

    def calculate_cell_bound(self) -> Tuple[float, float, float, float]:
        """
        Calculate the bounding box of the square cell.

        Computes the minimum and maximum longitude and latitude values
        from all vertices to create a bounding rectangle for spatial indexing.

        Returns:
            Tuple[float, float, float, float]: Bounding box as (min_lon, min_lat, max_lon, max_lat)
        """
        if not self.aVertex or len(self.aVertex) == 0:
            raise ValueError("Cannot calculate bounds: no vertices defined")

        dLat_min = 90.0
        dLat_max = -90.0
        dLon_min = 180.0
        dLon_max = -180.0

        for vertex in self.aVertex:
            if hasattr(vertex, "dLongitude_degree") and hasattr(
                vertex, "dLatitude_degree"
            ):
                dLon_max = max(dLon_max, vertex.dLongitude_degree)
                dLon_min = min(dLon_min, vertex.dLongitude_degree)
                dLat_max = max(dLat_max, vertex.dLatitude_degree)
                dLat_min = min(dLat_min, vertex.dLatitude_degree)
            else:
                raise ValueError(
                    "Vertex must have dLongitude_degree and dLatitude_degree attributes"
                )

        self.pBound = (dLon_min, dLat_min, dLon_max, dLat_max)
        return self.pBound

    def has_this_edge(self, pEdge_in) -> int:
        """
        Check whether the square contains a specific edge.

        Args:
            pEdge_in: The edge to be checked

        Returns:
            int: 1 if found, 0 if not found

        Raises:
            ValueError: If pEdge_in is None
        """
        if pEdge_in is None:
            raise ValueError("Edge to check cannot be None")

        iFlag_found = 0
        for pEdge in self.aEdge:
            if hasattr(pEdge, "is_overlap") and pEdge.is_overlap(pEdge_in):
                iFlag_found = 1
                break

        return iFlag_found

    def which_edge_cross_this_vertex(self, pVertex_in) -> Tuple[int, Optional[Any]]:
        """
        Find which edge overlaps with a vertex.

        Args:
            pVertex_in: The vertex to be checked

        Returns:
            Tuple[int, Optional[Any]]: (1 if found with edge object, 0 if not found with None)

        Raises:
            ValueError: If pVertex_in is None
        """
        if pVertex_in is None:
            raise ValueError("Vertex to check cannot be None")

        iFlag_found = 0
        pEdge_out = None

        for pEdge in self.aEdge:
            if hasattr(pEdge, "check_vertex_on_edge"):
                iFlag, dummy, diff = pEdge.check_vertex_on_edge(pVertex_in)
                if iFlag == 1:
                    iFlag_found = 1
                    pEdge_out = pEdge
                    break

        return iFlag_found, pEdge_out

    def calculate_cell_area(self) -> float:
        """
        Calculate the area of the square cell.

        Uses the polygon area calculation based on the vertex coordinates.
        For a perfect square, this should give area = side_length².

        Returns:
            float: The area in square meters

        Raises:
            ValueError: If vertices are not properly defined
        """
        if not self.aVertex or len(self.aVertex) == 0:
            raise ValueError("Cannot calculate area: no vertices defined")

        lons = []
        lats = []

        for vertex in self.aVertex:
            if hasattr(vertex, "dLongitude_degree") and hasattr(
                vertex, "dLatitude_degree"
            ):
                lons.append(vertex.dLongitude_degree)
                lats.append(vertex.dLatitude_degree)
            else:
                raise ValueError(
                    "Vertex must have dLongitude_degree and dLatitude_degree attributes"
                )

        if len(lons) != 4 or len(lats) != 4:
            raise ValueError(f"Square must have exactly 4 vertices, got {len(lons)}")

        self.dArea = calculate_polygon_area(lons, lats)
        return self.dArea

    def calculate_edge_length(self) -> float:
        """
        Calculate the effective side length of the square cell.

        For a square cell, the effective length is the square root of the area,
        which represents the side length of an equivalent square.

        Returns:
            float: The effective side length in meters

        Raises:
            ValueError: If area is negative or not calculated
        """
        if self.dArea < 0:
            raise ValueError("Cannot calculate edge length: area is negative")
        if self.dArea == 0:
            # Calculate area first if not done
            self.calculate_cell_area()

        dLength_edge = np.sqrt(self.dArea)
        self.dLength = dLength_edge
        return dLength_edge

    def share_edge(self, other: "pysquare") -> int:
        """
        Check whether this square cell shares an edge with another square cell.

        Two square cells share an edge if any of their boundary edges overlap.
        This is useful for determining adjacency in mesh topology.

        Args:
            other (pysquare): The other square cell to check

        Returns:
            int: 1 if cells share an edge, 0 if not

        Raises:
            TypeError: If other is not a pysquare instance
            ValueError: If other is None
        """
        if other is None:
            raise ValueError("Other cell cannot be None")
        if not isinstance(other, pysquare):
            raise TypeError(
                f"Other cell must be a pysquare instance, got {type(other)}"
            )

        iFlag_share = 0
        for pEdge in self.aEdge:
            for pEdge2 in other.aEdge:
                if hasattr(pEdge, "is_overlap") and hasattr(pEdge2, "is_overlap"):
                    if pEdge.is_overlap(pEdge2) == 1:
                        iFlag_share = 1
                        break
            if iFlag_share == 1:
                break

        return iFlag_share

    def get_square_dimensions(self) -> Tuple[float, float]:
        """
        Get the dimensions of the square cell.

        For a perfect square, both dimensions should be equal to the side length.
        For rectangular cells, returns width and height.

        Returns:
            Tuple[float, float]: (width, height) of the square in meters
        """
        if self.pBound is None:
            self.calculate_cell_bound()

        dLon_min, dLat_min, dLon_max, dLat_max = self.pBound

        # Convert coordinate differences to approximate meters
        # This is a rough approximation - for precise calculations,
        # geodetic computations should be used
        dWidth_degrees = dLon_max - dLon_min
        dHeight_degrees = dLat_max - dLat_min

        # Approximate conversion (111 km per degree)
        dWidth_meters = (
            dWidth_degrees * 111000 * np.cos(np.radians(self.dLatitude_center_degree))
        )
        dHeight_meters = dHeight_degrees * 111000

        return (dWidth_meters, dHeight_meters)

    def is_regular_square(self, tolerance: float = 0.01) -> bool:
        """
        Check if the cell is a regular square (all sides equal, all angles 90°).

        Args:
            tolerance (float): Tolerance for side length comparison (default: 0.01)

        Returns:
            bool: True if the cell is a regular square within tolerance
        """
        if len(self.aEdge) != 4:
            return False

        # Check if all edges have similar length
        edge_lengths = []
        for edge in self.aEdge:
            if hasattr(edge, "calculate_length"):
                edge_lengths.append(edge.calculate_length())
            elif hasattr(edge, "dLength"):
                edge_lengths.append(edge.dLength)
            else:
                # Cannot determine edge lengths
                return False

        if len(edge_lengths) != 4:
            return False

        # Check if all sides are approximately equal
        mean_length = np.mean(edge_lengths)
        for length in edge_lengths:
            if abs(length - mean_length) / mean_length > tolerance:
                return False

        return True

    def get_corner_coordinates(self) -> List[Tuple[float, float]]:
        """
        Get the coordinates of all four corners of the square.

        Returns:
            List[Tuple[float, float]]: List of (longitude, latitude) tuples for each corner
        """
        corners = []
        for vertex in self.aVertex:
            if hasattr(vertex, "dLongitude_degree") and hasattr(
                vertex, "dLatitude_degree"
            ):
                corners.append((vertex.dLongitude_degree, vertex.dLatitude_degree))
            else:
                raise ValueError(
                    "Vertex must have dLongitude_degree and dLatitude_degree attributes"
                )

        return corners

    def is_valid(self) -> bool:
        """
        Check if the square cell has valid attributes.

        A square cell is considered valid if:
        - It has valid center coordinates
        - It has exactly 4 boundary vertices and edges
        - It has a positive area
        - It has a valid cell ID
        - All vertices and edges are properly defined

        Returns:
            bool: True if cell has valid square attributes
        """
        # Check parent class validity
        if not super().is_valid():
            return False

        # Check square-specific requirements
        has_four_edges = len(self.aEdge) == 4
        has_four_vertices = len(self.aVertex) == 4
        has_positive_side_length = self.dLength > 0

        # Check that all vertices have required attributes
        vertices_valid = True
        for vertex in self.aVertex:
            if not (
                hasattr(vertex, "dLongitude_degree")
                and hasattr(vertex, "dLatitude_degree")
            ):
                vertices_valid = False
                break

        return (
            has_four_edges
            and has_four_vertices
            and has_positive_side_length
            and vertices_valid
        )

    def copy(self) -> "pysquare":
        """
        Create a deep copy of the square cell.

        Returns:
            pysquare: A new square cell object with the same attributes
        """
        # Create new square with same basic parameters
        new_square = pysquare(
            self.dLongitude_center_degree,
            self.dLatitude_center_degree,
            self.aEdge.copy(),
            self.aVertex.copy(),
        )

        # Copy all attributes from parent
        for attr_name, attr_value in self.__dict__.items():
            if hasattr(new_square, attr_name):
                if isinstance(attr_value, list) and attr_value is not None:
                    setattr(new_square, attr_name, attr_value.copy())
                else:
                    setattr(new_square, attr_name, attr_value)

        return new_square

    def tojson(self) -> str:
        """
        Convert square cell object to a JSON string.

        Serializes all cell attributes including neighbor information and
        geometric properties. Uses the custom SquareClassEncoder to handle
        numpy data types and complex objects.

        Returns:
            str: JSON string representation of the square cell

        Example:
            >>> square = pysquare(-77.0, 38.0, edges, vertices)
            >>> json_str = square.tojson()
        """
        aSkip = [
            "aEdge",
            "aFlowline",
            "dLongitude_radian",
            "dLatitude_radian",
            "wkt",
            "pPoint_start",
            "pPoint_end",
            "pBound",
        ]

        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)

        sJson = json.dumps(
            obj, sort_keys=True, indent=4, ensure_ascii=True, cls=SquareClassEncoder
        )
        return sJson
