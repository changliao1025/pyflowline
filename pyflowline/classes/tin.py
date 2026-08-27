"""
PyTIN class for representing Triangulated Irregular Network mesh cells in flowline networks.

This module defines the pytin class which extends pymeshcell to provide
TIN-specific mesh cell functionality for flowline network modeling.
TIN (Triangulated Irregular Network) uses triangular cells for unstructured
mesh representation.
"""

import json
from json import JSONEncoder
from typing import List, Optional, Any, Tuple, Dict
import numpy as np

from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area
from pyearth.gis.geometry.calculate_spherical_triangle_area import (
    calculate_spherical_triangle_area,
)

from pyearthriver.classes.vertex import pyvertex
from pyearthriver.classes.edge import pyedge
from pyearthriver.classes.flowline import pyflowline
from pyflowline.classes.meshcell import pymeshcell


class TINClassEncoder(JSONEncoder):
    """
    Custom JSON encoder for pytin objects.

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
        if isinstance(obj, pytin):
            return obj.lCellID
        return JSONEncoder.default(self, obj)


class pytin(pymeshcell):
    """
    Triangulated Irregular Network (TIN) mesh cell class for flowline network representation.

    Extends the pymeshcell class to provide TIN-specific functionality for
    mesh cells in flowline topology and stream network analysis. A TIN cell
    represents a triangular area with exactly 3 vertices and 3 edges.

    TIN (Triangulated Irregular Network) uses triangular cells that can adapt
    to complex terrain and irregular boundaries. This makes TIN particularly
    suitable for high-resolution terrain modeling and hydrological applications
    where terrain complexity needs to be preserved.

    Attributes:
        All attributes from pymeshcell plus TIN-specific properties.
        The TIN cell is always triangular with exactly 3 edges and 3 vertices.
        use_spherical_calculation (bool): Whether to use spherical triangle area calculation
        triangle_quality (float): Quality metric for triangle shape (0-1, 1 is equilateral)

    Args:
        dLon (float): The longitude of the TIN triangle center in degrees
        dLat (float): The latitude of the TIN triangle center in degrees
        aEdge (List): A list of exactly 3 edges that define the triangle boundary
        aVertex (List): A list of exactly 3 vertices that define the triangle boundary

    Raises:
        ValueError: If the cell doesn't have exactly 3 edges and 3 vertices
        TypeError: If input parameters are not of the expected types

    Example:
        >>> center_lon, center_lat = -77.0, 38.0
        >>> edges = [edge1, edge2, edge3]  # List of 3 edges
        >>> vertices = [vertex1, vertex2, vertex3]  # List of 3 vertices
        >>> tin_cell = pytin(center_lon, center_lat, edges, vertices)
        >>> tin_cell.set_cell_id(100)
        >>> print(tin_cell)
        pytin(ID=100, Center=(-77.0, 38.0), Area=1234.56m², Quality=0.85)
    """

    def __init__(self, dLon: float, dLat: float, aEdge: List, aVertex: List) -> None:
        """
        Initialize a TIN triangular mesh cell object.

        Args:
            dLon (float): The longitude of the TIN triangle center in degrees
            dLat (float): The latitude of the TIN triangle center in degrees
            aEdge (List): A list of exactly 3 edges that define the triangle boundary
            aVertex (List): A list of exactly 3 vertices that define the triangle boundary

        Raises:
            ValueError: If the cell doesn't have exactly 3 edges and vertices
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

        # Validate TIN geometry requirements (exactly 3 edges and vertices)
        nEdge = len(aEdge)
        nVertex = len(aVertex)

        if nEdge != 3:
            raise ValueError(f"TIN cell must have exactly 3 edges, got {nEdge}")
        if nVertex != 3:
            raise ValueError(f"TIN cell must have exactly 3 vertices, got {nVertex}")

        # Validate edges and vertices are not None
        for i, edge in enumerate(aEdge):
            if edge is None:
                raise ValueError(f"Edge at index {i} cannot be None")

        for i, vertex in enumerate(aVertex):
            if vertex is None:
                raise ValueError(f"Vertex at index {i} cannot be None")

        # Initialize parent class
        super().__init__(dLon, dLat, aEdge, aVertex)

        # TIN-specific initialization
        self.nEdge = 3
        self.nVertex = 3

        # TIN-specific attributes
        self.use_spherical_calculation: bool = True
        self.triangle_quality: float = -1.0  # Will be calculated

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
        if self.dArea > 0:
            self.calculate_edge_length()
            self.calculate_triangle_quality()

    def __repr__(self) -> str:
        """
        Return a detailed string representation of the TIN triangle cell.

        Returns:
            str: Detailed representation including cell ID, center coordinates, and area
        """
        return (
            f"pytin(ID={self.lCellID}, "
            f"Center=({self.dLongitude_center_degree:.6f}, {self.dLatitude_center_degree:.6f}), "
            f"Area={self.dArea:.2f}m², Quality={self.triangle_quality:.3f}, "
            f"EdgeLength={self.dLength:.2f}m, Flowlines={self.nFlowline}, "
            f"Neighbors={self.nNeighbor})"
        )

    def __str__(self) -> str:
        """
        Return a concise string representation of the TIN triangle cell.

        Returns:
            str: Concise representation with ID, center coordinates, and area
        """
        return (
            f"pytin(ID={self.lCellID}, "
            f"Center=({self.dLongitude_center_degree:.6f}, {self.dLatitude_center_degree:.6f}), "
            f"Area={self.dArea:.2f}m², Quality={self.triangle_quality:.3f})"
        )

    def set_spherical_calculation(self, use_spherical: bool) -> None:
        """
        Set whether to use spherical triangle area calculation.

        Args:
            use_spherical (bool): True to use spherical calculation, False for planar

        Raises:
            TypeError: If use_spherical is not a boolean
        """
        if not isinstance(use_spherical, bool):
            raise TypeError(
                f"use_spherical must be a boolean, got {type(use_spherical)}"
            )

        self.use_spherical_calculation = use_spherical
        # Recalculate area with new method
        self.calculate_cell_area()

    def calculate_cell_bound(self) -> Tuple[float, float, float, float]:
        """
        Calculate the bounding box of the TIN triangle cell.

        Computes the minimum and maximum longitude and latitude values
        from all three vertices to create a bounding rectangle for spatial indexing.

        Returns:
            Tuple[float, float, float, float]: Bounding box as (min_lon, min_lat, max_lon, max_lat)
        """
        if not self.aVertex or len(self.aVertex) != 3:
            raise ValueError(
                "Cannot calculate bounds: triangle must have exactly 3 vertices"
            )

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
        Check whether the TIN triangle contains a specific edge.

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
        Calculate the area of the TIN triangle cell.

        Uses either spherical triangle area calculation (for large triangles
        or global applications) or planar polygon area calculation based on
        the use_spherical_calculation setting.

        Returns:
            float: The area in square meters

        Raises:
            ValueError: If vertices are not properly defined
        """
        if not self.aVertex or len(self.aVertex) != 3:
            raise ValueError(
                "Cannot calculate area: triangle must have exactly 3 vertices"
            )

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

        if self.use_spherical_calculation:
            self.dArea = calculate_spherical_triangle_area(lons, lats)
        else:
            self.dArea = calculate_polygon_area(lons, lats)

        return self.dArea

    def calculate_edge_length(self) -> float:
        """
        Calculate the characteristic edge length of the TIN triangle.

        For a triangle, the characteristic length is calculated as the
        square root of the area, which gives a representative size measure.

        Returns:
            float: The characteristic edge length in meters

        Raises:
            ValueError: If area is negative or not calculated
        """
        if self.dArea < 0:
            raise ValueError("Cannot calculate edge length: area is negative")
        if self.dArea == 0:
            # Calculate area first if not done
            self.calculate_cell_area()

        self.dLength = np.sqrt(self.dArea)
        return self.dLength

    def calculate_triangle_quality(self) -> float:
        """
        Calculate the quality metric of the triangle.

        Triangle quality is measured as the ratio of the area to the
        area of an equilateral triangle with the same perimeter.
        Quality ranges from 0 to 1, where 1 represents a perfect equilateral triangle.

        Returns:
            float: Triangle quality metric (0-1, higher is better)

        Raises:
            ValueError: If edges cannot be calculated
        """
        try:
            # Calculate edge lengths if available
            edge_lengths = []
            for edge in self.aEdge:
                if hasattr(edge, "calculate_length"):
                    edge_lengths.append(edge.calculate_length())
                elif hasattr(edge, "dLength"):
                    edge_lengths.append(edge.dLength)

            if len(edge_lengths) == 3:
                # Calculate perimeter
                perimeter = sum(edge_lengths)

                # Area of equilateral triangle with same perimeter
                side_length = perimeter / 3.0
                equilateral_area = (np.sqrt(3.0) / 4.0) * side_length**2

                # Quality metric
                if equilateral_area > 0:
                    self.triangle_quality = self.dArea / equilateral_area
                else:
                    self.triangle_quality = 0.0
            else:
                # Fallback: use aspect ratio approximation
                corners = self.get_corner_coordinates()
                if len(corners) == 3:
                    # Calculate distances between vertices
                    d1 = self._calculate_distance(corners[0], corners[1])
                    d2 = self._calculate_distance(corners[1], corners[2])
                    d3 = self._calculate_distance(corners[2], corners[0])

                    # Aspect ratio approximation
                    max_side = max(d1, d2, d3)
                    min_side = min(d1, d2, d3)

                    if max_side > 0:
                        self.triangle_quality = min_side / max_side
                    else:
                        self.triangle_quality = 0.0
                else:
                    self.triangle_quality = 0.0

        except Exception:
            self.triangle_quality = 0.0

        return self.triangle_quality

    def _calculate_distance(
        self, point1: Tuple[float, float], point2: Tuple[float, float]
    ) -> float:
        """
        Calculate approximate distance between two points in degrees.

        Args:
            point1: (longitude, latitude) tuple
            point2: (longitude, latitude) tuple

        Returns:
            float: Approximate distance
        """
        lon1, lat1 = point1
        lon2, lat2 = point2
        return np.sqrt((lon2 - lon1) ** 2 + (lat2 - lat1) ** 2)

    def share_edge(self, other: "pytin") -> int:
        """
        Check whether this TIN triangle shares an edge with another TIN triangle.

        Two TIN triangles share an edge if any of their boundary edges overlap.
        This is essential for determining adjacency in triangular meshes.

        Args:
            other (pytin): The other TIN triangle to check

        Returns:
            int: 1 if triangles share an edge, 0 if not

        Raises:
            TypeError: If other is not a pytin instance
            ValueError: If other is None
        """
        if other is None:
            raise ValueError("Other cell cannot be None")
        if not isinstance(other, pytin):
            raise TypeError(f"Other cell must be a pytin instance, got {type(other)}")

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

    def get_tin_properties(self) -> Dict[str, Any]:
        """
        Get properties specific to the TIN triangle cell.

        Returns:
            dict: Dictionary containing TIN-specific properties including:
                - area: Triangle area
                - quality: Triangle quality metric (0-1)
                - use_spherical: Whether using spherical calculations
                - perimeter: Triangle perimeter (if calculable)
                - aspect_ratio: Ratio of longest to shortest side
                - centroid: Triangle centroid coordinates
        """
        properties = {
            "area": self.dArea,
            "quality": self.triangle_quality,
            "use_spherical": self.use_spherical_calculation,
            "characteristic_length": self.dLength,
        }

        # Calculate additional properties if possible
        try:
            corners = self.get_corner_coordinates()
            if len(corners) == 3:
                # Calculate centroid
                centroid_lon = sum(corner[0] for corner in corners) / 3.0
                centroid_lat = sum(corner[1] for corner in corners) / 3.0
                properties["centroid"] = (centroid_lon, centroid_lat)

                # Calculate edge lengths and perimeter
                d1 = self._calculate_distance(corners[0], corners[1])
                d2 = self._calculate_distance(corners[1], corners[2])
                d3 = self._calculate_distance(corners[2], corners[0])

                properties["perimeter"] = d1 + d2 + d3
                properties["edge_lengths"] = [d1, d2, d3]

                # Aspect ratio
                max_side = max(d1, d2, d3)
                min_side = min(d1, d2, d3)
                properties["aspect_ratio"] = (
                    max_side / min_side if min_side > 0 else float("inf")
                )

        except Exception:
            # If calculation fails, use defaults
            properties["centroid"] = (
                self.dLongitude_center_degree,
                self.dLatitude_center_degree,
            )
            properties["perimeter"] = 0.0
            properties["aspect_ratio"] = 1.0

        return properties

    def is_degenerate(self, tolerance: float = 1e-10) -> bool:
        """
        Check if the triangle is degenerate (has zero or very small area).

        Args:
            tolerance (float): Minimum area threshold for non-degenerate triangle

        Returns:
            bool: True if triangle is degenerate
        """
        return self.dArea < tolerance

    def is_well_shaped(self, quality_threshold: float = 0.3) -> bool:
        """
        Check if the triangle has good shape quality.

        Args:
            quality_threshold (float): Minimum quality threshold (0-1)

        Returns:
            bool: True if triangle quality exceeds threshold
        """
        if self.triangle_quality < 0:
            self.calculate_triangle_quality()

        return self.triangle_quality >= quality_threshold

    def get_corner_coordinates(self) -> List[Tuple[float, float]]:
        """
        Get the coordinates of all three corners of the triangle.

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

    def get_triangle_angles(self) -> List[float]:
        """
        Calculate the three interior angles of the triangle.

        Returns:
            List[float]: List of three angles in degrees

        Note:
            This is an approximation using spherical geometry for small triangles
        """
        corners = self.get_corner_coordinates()
        if len(corners) != 3:
            return [0.0, 0.0, 0.0]

        angles = []
        try:
            # Calculate angles using law of cosines approximation
            for i in range(3):
                p1 = corners[i]
                p2 = corners[(i + 1) % 3]
                p3 = corners[(i + 2) % 3]

                # Calculate side lengths
                a = self._calculate_distance(p2, p3)  # Opposite side
                b = self._calculate_distance(p1, p3)  # Adjacent side 1
                c = self._calculate_distance(p1, p2)  # Adjacent side 2

                # Law of cosines: cos(angle) = (b² + c² - a²) / (2bc)
                if b > 0 and c > 0:
                    cos_angle = (b * b + c * c - a * a) / (2 * b * c)
                    # Clamp to valid range for acos
                    cos_angle = max(-1.0, min(1.0, cos_angle))
                    angle_rad = np.arccos(cos_angle)
                    angle_deg = np.degrees(angle_rad)
                    angles.append(angle_deg)
                else:
                    angles.append(0.0)

        except Exception:
            angles = [60.0, 60.0, 60.0]  # Default to equilateral triangle

        return angles

    def is_valid(self) -> bool:
        """
        Check if the TIN triangle has valid attributes.

        A TIN triangle is considered valid if:
        - It has valid center coordinates
        - It has exactly 3 boundary vertices and edges
        - It has a positive area
        - It has a valid cell ID
        - All vertices and edges are properly defined
        - Triangle is not degenerate

        Returns:
            bool: True if cell has valid TIN attributes
        """
        # Check parent class validity
        if not super().is_valid():
            return False

        # Check TIN-specific requirements
        has_three_edges = len(self.aEdge) == 3
        has_three_vertices = len(self.aVertex) == 3
        has_positive_area = self.dArea > 0
        is_not_degenerate = not self.is_degenerate()

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
            has_three_edges
            and has_three_vertices
            and has_positive_area
            and is_not_degenerate
            and vertices_valid
        )

    def copy(self) -> "pytin":
        """
        Create a deep copy of the TIN triangle.

        Returns:
            pytin: A new TIN triangle object with the same attributes
        """
        # Create new TIN triangle with same basic parameters
        new_tin = pytin(
            self.dLongitude_center_degree,
            self.dLatitude_center_degree,
            self.aEdge.copy(),
            self.aVertex.copy(),
        )

        # Copy all attributes from parent and this class
        for attr_name, attr_value in self.__dict__.items():
            if hasattr(new_tin, attr_name):
                if isinstance(attr_value, list) and attr_value is not None:
                    setattr(new_tin, attr_name, attr_value.copy())
                else:
                    setattr(new_tin, attr_name, attr_value)

        return new_tin

    def tojson(self) -> str:
        """
        Convert TIN triangle object to a JSON string.

        Serializes all cell attributes including neighbor information and
        geometric properties. Uses the custom TINClassEncoder to handle
        numpy data types and complex objects.

        Returns:
            str: JSON string representation of the TIN triangle

        Example:
            >>> tin_cell = pytin(-77.0, 38.0, edges, vertices)
            >>> json_str = tin_cell.tojson()
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
            obj, sort_keys=True, indent=4, ensure_ascii=True, cls=TINClassEncoder
        )
        return sJson
