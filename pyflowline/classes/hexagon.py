"""
PyHexagon class for representing hexagonal mesh cells in flowline networks.

This module defines the pyhexagon class which extends pymeshcell to provide
hexagon-specific mesh cell functionality for flowline network modeling.
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


class HexagonClassEncoder(JSONEncoder):
    """
    Custom JSON encoder for pyhexagon objects.

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
        if pyvertex and isinstance(obj, pyvertex):
            return json.loads(obj.tojson())
        if pyedge and isinstance(obj, pyedge):
            return getattr(obj, "lEdgeID", str(obj))
        if pyflowline and isinstance(obj, pyflowline):
            return getattr(obj, "lFlowlineID", str(obj))
        if isinstance(obj, pyhexagon):
            return obj.lCellID
        return JSONEncoder.default(self, obj)


class pyhexagon(pymeshcell):
    """
    Hexagonal mesh cell class for flowline network representation.

    Extends the pymeshcell class to provide hexagon-specific functionality for
    mesh cells in flowline topology and stream network analysis. A hexagon cell
    represents a hexagonal polygonal area with 6 vertices and 6 edges.

    This class is designed for hexagonal grid-based mesh representations which
    are particularly useful for unstructured grid hydrological models due to
    their geometric properties and connectivity patterns.

    Attributes:
        All attributes from pymeshcell plus hexagon-specific properties.
        The hexagon cell is defined by exactly 6 edges and 6 vertices forming
        a regular or irregular hexagonal boundary.
        iFlag_watershed_boundary_burned (int): Flag for watershed boundary processing

    Args:
        dLon (float): The longitude of the hexagon center in degrees
        dLat (float): The latitude of the hexagon center in degrees
        aEdge (List): A list of exactly 6 edges that define the hexagon boundary
        aVertex (List): A list of exactly 6 vertices that define the hexagon boundary

    Raises:
        ValueError: If the cell doesn't have exactly 6 edges and 6 vertices
        TypeError: If input parameters are not of the expected types

    Example:
        >>> center_lon, center_lat = -77.0, 38.0
        >>> edges = [edge1, edge2, edge3, edge4, edge5, edge6]  # List of 6 edges
        >>> vertices = [v1, v2, v3, v4, v5, v6]  # List of 6 vertices
        >>> hexagon = pyhexagon(center_lon, center_lat, edges, vertices)
        >>> hexagon.set_cell_id(100)
        >>> print(hexagon)
        pyhexagon(ID=100, Center=(-77.0, 38.0), Area=1234.56m²)
    """

    def __init__(self, dLon: float, dLat: float, aEdge: List, aVertex: List) -> None:
        """
        Initialize a hexagonal mesh cell object.

        Args:
            dLon (float): The longitude of the hexagon center in degrees
            dLat (float): The latitude of the hexagon center in degrees
            aEdge (List): A list of exactly 6 edges that define the hexagon boundary
            aVertex (List): A list of exactly 6 vertices that define the hexagon boundary

        Raises:
            ValueError: If the cell doesn't have exactly 6 edges and vertices
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

        # Validate hexagon geometry requirements
        nEdge = len(aEdge)
        nVertex = len(aVertex)

        if nEdge != 6:
            raise ValueError(f"Hexagon cell must have exactly 6 edges, got {nEdge}")
        if nVertex != 6:
            raise ValueError(
                f"Hexagon cell must have exactly 6 vertices, got {nVertex}"
            )

        # Validate edges and vertices are not None
        for i, edge in enumerate(aEdge):
            if edge is None:
                raise ValueError(f"Edge at index {i} cannot be None")

        for i, vertex in enumerate(aVertex):
            if vertex is None:
                raise ValueError(f"Vertex at index {i} cannot be None")

        # Initialize parent class
        super().__init__(dLon, dLat, aEdge, aVertex)

        # Hexagon-specific initialization
        self.nEdge = 6
        self.nVertex = 6  # Note: Different from nPoint which includes closing vertex

        # Hexagon-specific attributes
        self.iFlag_watershed_boundary_burned: int = 0

        # Create center vertex if pyvertex is available
        if pyvertex:
            pVertex_params = {
                "dLongitude_degree": self.dLongitude_center_degree,
                "dLatitude_degree": self.dLatitude_center_degree,
            }
            self.pVertex_center = pyvertex(pVertex_params)
        else:
            self.pVertex_center = None

        # Calculate initial properties
        self.calculate_cell_area()
        self.calculate_edge_length()

    def __repr__(self) -> str:
        """
        Return a detailed string representation of the hexagon cell.

        Returns:
            str: Detailed representation including cell ID, center coordinates, and area
        """
        return (
            f"pyhexagon(ID={self.lCellID}, "
            f"Center=({self.dLongitude_center_degree:.6f}, {self.dLatitude_center_degree:.6f}), "
            f"Area={self.dArea:.2f}m², EdgeLength={self.dLength:.2f}m, "
            f"Flowlines={self.nFlowline}, Neighbors={self.nNeighbor})"
        )

    def __str__(self) -> str:
        """
        Return a concise string representation of the hexagon cell.

        Returns:
            str: Concise representation with ID, center coordinates, and area
        """
        return (
            f"pyhexagon(ID={self.lCellID}, "
            f"Center=({self.dLongitude_center_degree:.6f}, {self.dLatitude_center_degree:.6f}), "
            f"Area={self.dArea:.2f}m²)"
        )

    def calculate_cell_bound(self) -> Tuple[float, float, float, float]:
        """
        Calculate the bounding box of the hexagon cell.

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
        Check whether the hexagon contains a specific edge.

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
        Calculate the area of the hexagon cell.

        Uses the polygon area calculation based on the vertex coordinates.
        For a regular hexagon, this should give area = (3√3/2) * side_length².

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

        if len(lons) != 6 or len(lats) != 6:
            raise ValueError(f"Hexagon must have exactly 6 vertices, got {len(lons)}")

        self.dArea = calculate_polygon_area(lons, lats)
        return self.dArea

    def calculate_edge_length(self) -> float:
        """
        Calculate the effective edge length of the hexagon cell.

        For a regular hexagon, the edge length can be derived from the area using:
        edge_length = sqrt(2 * area / (3 * sqrt(3)))

        This formula comes from the regular hexagon area formula:
        area = (3 * sqrt(3) / 2) * edge_length²

        Returns:
            float: The effective edge length in meters

        Raises:
            ValueError: If area is negative or not calculated
        """
        if self.dArea < 0:
            raise ValueError("Cannot calculate edge length: area is negative")
        if self.dArea == 0:
            # Calculate area first if not done
            self.calculate_cell_area()

        # Formula derived from regular hexagon area: A = (3√3/2) * s²
        # Solving for s: s = √(2A / (3√3))
        dLength_edge = np.sqrt(2.0 * self.dArea / (3.0 * np.sqrt(3.0)))
        self.dLength = dLength_edge
        return dLength_edge

    def share_edge(self, other: "pyhexagon") -> int:
        """
        Check whether this hexagon cell shares an edge with another hexagon cell.

        Two hexagon cells share an edge if any of their boundary edges overlap.
        This is useful for determining adjacency in mesh topology.

        Args:
            other (pyhexagon): The other hexagon cell to check

        Returns:
            int: 1 if cells share an edge, 0 if not

        Raises:
            TypeError: If other is not a pyhexagon instance
            ValueError: If other is None
        """
        if other is None:
            raise ValueError("Other cell cannot be None")
        if not isinstance(other, pyhexagon):
            raise TypeError(
                f"Other cell must be a pyhexagon instance, got {type(other)}"
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

    def get_hexagon_properties(self) -> dict:
        """
        Get geometric properties specific to the hexagon cell.

        Returns:
            dict: Dictionary containing hexagon-specific properties including:
                - edge_length: Effective edge length
                - area: Cell area
                - perimeter: Approximate perimeter (6 * edge_length for regular hexagon)
                - circumradius: Distance from center to vertex (for regular hexagon)
                - inradius: Distance from center to edge midpoint (for regular hexagon)
        """
        properties = {
            "edge_length": self.dLength,
            "area": self.dArea,
            "perimeter": 6.0 * self.dLength,  # Approximate for regular hexagon
            "circumradius": self.dLength,  # For regular hexagon
            "inradius": self.dLength * np.sqrt(3.0) / 2.0,  # For regular hexagon
        }
        return properties

    def is_regular_hexagon(self, tolerance: float = 0.01) -> bool:
        """
        Check if the cell is a regular hexagon (all sides equal, all angles 120°).

        Args:
            tolerance (float): Tolerance for side length comparison (default: 0.01)

        Returns:
            bool: True if the cell is a regular hexagon within tolerance
        """
        if len(self.aEdge) != 6:
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

        if len(edge_lengths) != 6:
            return False

        # Check if all sides are approximately equal
        mean_length = np.mean(edge_lengths)
        for length in edge_lengths:
            if abs(length - mean_length) / mean_length > tolerance:
                return False

        return True

    def get_corner_coordinates(self) -> List[Tuple[float, float]]:
        """
        Get the coordinates of all six corners of the hexagon.

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

    def set_watershed_boundary_flag(self, iFlag: int) -> None:
        """
        Set the watershed boundary flag for the hexagon cell.

        Args:
            iFlag (int): Flag value (0 or 1)

        Raises:
            TypeError: If iFlag is not an integer
            ValueError: If iFlag is not 0 or 1
        """
        if not isinstance(iFlag, (int, np.integer)):
            raise TypeError(f"Flag must be an integer, got {type(iFlag)}")
        if iFlag not in [0, 1]:
            raise ValueError(f"Flag must be 0 or 1, got {iFlag}")

        self.iFlag_watershed_boundary_burned = int(iFlag)

    def is_watershed_boundary(self) -> bool:
        """
        Check if the hexagon cell is on a watershed boundary.

        Returns:
            bool: True if cell is on watershed boundary
        """
        return self.iFlag_watershed_boundary_burned == 1

    def is_valid(self) -> bool:
        """
        Check if the hexagon cell has valid attributes.

        A hexagon cell is considered valid if:
        - It has valid center coordinates
        - It has exactly 6 boundary vertices and edges
        - It has a positive area
        - It has a valid cell ID
        - All vertices and edges are properly defined

        Returns:
            bool: True if cell has valid hexagon attributes
        """
        # Check parent class validity
        if not super().is_valid():
            return False

        # Check hexagon-specific requirements
        has_six_edges = len(self.aEdge) == 6
        has_six_vertices = len(self.aVertex) == 6
        has_positive_edge_length = self.dLength > 0

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
            has_six_edges
            and has_six_vertices
            and has_positive_edge_length
            and vertices_valid
        )

    def copy(self) -> "pyhexagon":
        """
        Create a deep copy of the hexagon cell.

        Returns:
            pyhexagon: A new hexagon cell object with the same attributes
        """
        # Create new hexagon with same basic parameters
        new_hexagon = pyhexagon(
            self.dLongitude_center_degree,
            self.dLatitude_center_degree,
            self.aEdge.copy(),
            self.aVertex.copy(),
        )

        # Copy all attributes from parent and this class
        for attr_name, attr_value in self.__dict__.items():
            if hasattr(new_hexagon, attr_name):
                if isinstance(attr_value, list) and attr_value is not None:
                    setattr(new_hexagon, attr_name, attr_value.copy())
                else:
                    setattr(new_hexagon, attr_name, attr_value)

        return new_hexagon

    def tojson(self) -> str:
        """
        Convert hexagon cell object to a JSON string.

        Serializes all cell attributes including neighbor information and
        geometric properties. Uses the custom HexagonClassEncoder to handle
        numpy data types and complex objects.

        Returns:
            str: JSON string representation of the hexagon cell

        Example:
            >>> hexagon = pyhexagon(-77.0, 38.0, edges, vertices)
            >>> json_str = hexagon.tojson()
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
            obj, sort_keys=True, indent=4, ensure_ascii=True, cls=HexagonClassEncoder
        )
        return sJson
