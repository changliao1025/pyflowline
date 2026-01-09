"""
PyDGGRID class for representing Discrete Global Grid mesh cells in flowline networks.

This module defines the pydggrid class which extends pymeshcell to provide
DGGRID-specific mesh cell functionality for flowline network modeling.
DGGRID (Discrete Global Grid) systems use hierarchical grids with various
cell shapes including hexagons, triangles, and diamonds.
"""

import json
from json import JSONEncoder
from typing import List, Optional, Any, Tuple, Dict
import numpy as np

from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area

from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.edge import pyedge
from pyflowline.classes.flowline import pyflowline
from pyflowline.classes.meshcell import pymeshcell


class DggridClassEncoder(JSONEncoder):
    """
    Custom JSON encoder for pydggrid objects.

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
        if isinstance(obj, pydggrid):
            return obj.lCellID
        return JSONEncoder.default(self, obj)


class pydggrid(pymeshcell):
    """
    Discrete Global Grid (DGGRID) mesh cell class for flowline network representation.

    Extends the pymeshcell class to provide DGGRID-specific functionality for
    mesh cells in flowline topology and stream network analysis. A DGGRID cell
    represents a polygonal area in a hierarchical discrete global grid system.

    DGGRID systems use multi-resolution grids that can have various cell shapes
    including hexagons, triangles, diamonds, and other polygons depending on
    the grid configuration and hierarchical level. This flexibility makes DGGRID
    suitable for global modeling at multiple scales.

    Attributes:
        All attributes from pymeshcell plus DGGRID-specific properties.
        The DGGRID cell can have variable numbers of edges and vertices
        depending on the grid system configuration.
        dggrid_level (int): Hierarchical level in the DGGRID system
        dggrid_id (str): Unique identifier within the DGGRID system
        dggrid_resolution (float): Resolution class of the DGGRID cell

    Args:
        dLon (float): The longitude of the DGGRID cell center in degrees
        dLat (float): The latitude of the DGGRID cell center in degrees
        aEdge (List): A list of edges that define the DGGRID cell boundary
        aVertex (List): A list of vertices that define the DGGRID cell boundary

    Raises:
        ValueError: If the cell doesn't have at least 3 edges and vertices
        TypeError: If input parameters are not of the expected types

    Example:
        >>> center_lon, center_lat = -77.0, 38.0
        >>> edges = [edge1, edge2, edge3, edge4, edge5, edge6]  # List of edges
        >>> vertices = [v1, v2, v3, v4, v5, v6]  # List of vertices
        >>> dggrid_cell = pydggrid(center_lon, center_lat, edges, vertices)
        >>> dggrid_cell.set_cell_id(100)
        >>> dggrid_cell.set_dggrid_properties(level=5, grid_id="N1234", resolution=1000.0)
        >>> print(dggrid_cell)
        pydggrid(ID=100, Center=(-77.0, 38.0), Level=5, Edges=6, Area=1234.56m²)
    """

    def __init__(self, dLon: float, dLat: float, aEdge: List, aVertex: List) -> None:
        """
        Initialize a DGGRID mesh cell object.

        Args:
            dLon (float): The longitude of the DGGRID cell center in degrees
            dLat (float): The latitude of the DGGRID cell center in degrees
            aEdge (List): A list of edges that define the DGGRID cell boundary
            aVertex (List): A list of vertices that define the DGGRID cell boundary

        Raises:
            ValueError: If the cell doesn't have at least 3 edges and vertices
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

        # Validate DGGRID geometry requirements (minimum 3 edges for polygon)
        nEdge = len(aEdge)
        nVertex = len(aVertex)

        if nEdge < 3:
            raise ValueError(f"DGGRID cell must have at least 3 edges, got {nEdge}")
        if nVertex < 3:
            raise ValueError(
                f"DGGRID cell must have at least 3 vertices, got {nVertex}"
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

        # DGGRID-specific initialization
        self.nEdge = nEdge
        self.nVertex = nVertex

        # DGGRID-specific attributes
        self.dggrid_level: int = -1
        self.dggrid_id: str = ""
        self.dggrid_resolution: float = -1.0
        self.dggrid_shape_type: str = "unknown"

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

        # Determine shape type based on number of edges
        self._determine_shape_type()

    def __repr__(self) -> str:
        """
        Return a detailed string representation of the DGGRID cell.

        Returns:
            str: Detailed representation including cell ID, center coordinates, level, and area
        """
        return (
            f"pydggrid(ID={self.lCellID}, "
            f"Center=({self.dLongitude_center_degree:.6f}, {self.dLatitude_center_degree:.6f}), "
            f"Level={self.dggrid_level}, Edges={self.nEdge}, Shape={self.dggrid_shape_type}, "
            f"Area={self.dArea:.2f}m², Resolution={self.dLength:.2f}m, "
            f"Flowlines={self.nFlowline}, Neighbors={self.nNeighbor})"
        )

    def __str__(self) -> str:
        """
        Return a concise string representation of the DGGRID cell.

        Returns:
            str: Concise representation with ID, center coordinates, level, and area
        """
        return (
            f"pydggrid(ID={self.lCellID}, "
            f"Center=({self.dLongitude_center_degree:.6f}, {self.dLatitude_center_degree:.6f}), "
            f"Level={self.dggrid_level}, Shape={self.dggrid_shape_type}, Area={self.dArea:.2f}m²)"
        )

    def _determine_shape_type(self) -> None:
        """
        Determine the shape type based on the number of edges.

        Sets the dggrid_shape_type attribute based on edge count.
        """
        shape_map = {
            3: "triangle",
            4: "quadrilateral",
            5: "pentagon",
            6: "hexagon",
            7: "heptagon",
            8: "octagon",
        }
        self.dggrid_shape_type = shape_map.get(
            self.nEdge, f"{self.nEdge}-sided polygon"
        )

    def set_dggrid_properties(
        self, level: int = -1, grid_id: str = "", resolution: float = -1.0
    ) -> None:
        """
        Set DGGRID-specific properties for the cell.

        Args:
            level (int): Hierarchical level in the DGGRID system
            grid_id (str): Unique identifier within the DGGRID system
            resolution (float): Resolution class of the DGGRID cell

        Raises:
            TypeError: If parameters are not of the expected types
            ValueError: If level or resolution are negative when provided
        """
        if not isinstance(level, (int, np.integer)):
            raise TypeError(f"Level must be an integer, got {type(level)}")
        if not isinstance(grid_id, str):
            raise TypeError(f"Grid ID must be a string, got {type(grid_id)}")
        if not isinstance(resolution, (int, float, np.number)):
            raise TypeError(f"Resolution must be a number, got {type(resolution)}")

        if level < -1:
            raise ValueError(f"Level must be -1 (unset) or positive, got {level}")
        if resolution < -1.0:
            raise ValueError(
                f"Resolution must be -1.0 (unset) or positive, got {resolution}"
            )

        self.dggrid_level = int(level)
        self.dggrid_id = str(grid_id)
        self.dggrid_resolution = float(resolution)

    def calculate_cell_bound(self) -> Tuple[float, float, float, float]:
        """
        Calculate the bounding box of the DGGRID cell.

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
        Check whether the DGGRID cell contains a specific edge.

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
        Calculate the area of the DGGRID cell.

        Uses the polygon area calculation based on the vertex coordinates.
        This method works for any polygon shape with 3 or more vertices.

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

        if len(lons) < 3:
            raise ValueError(
                f"DGGRID cell must have at least 3 vertices, got {len(lons)}"
            )

        self.dArea = calculate_polygon_area(lons, lats)
        return self.dArea

    def calculate_edge_length(self) -> float:
        """
        Calculate the effective edge length of the DGGRID cell.

        The calculation method depends on the shape type:
        - For hexagonal cells: Uses hexagon-specific formula
        - For other shapes: Uses square root of area as characteristic length

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

        if self.dggrid_shape_type == "hexagon":
            # Formula for regular hexagon: A = (3√3/2) * s²
            # Solving for s: s = √(2A / (3√3))
            dLength_edge = np.sqrt(2.0 * self.dArea / (3.0 * np.sqrt(3.0)))
        else:
            # For other shapes, use square root of area as characteristic length
            dLength_edge = np.sqrt(self.dArea)

        self.dLength = dLength_edge
        return dLength_edge

    def share_edge(self, other: "pydggrid") -> int:
        """
        Check whether this DGGRID cell shares an edge with another DGGRID cell.

        Two DGGRID cells share an edge if any of their boundary edges overlap.
        This is essential for determining adjacency in hierarchical grids.

        Args:
            other (pydggrid): The other DGGRID cell to check

        Returns:
            int: 1 if cells share an edge, 0 if not

        Raises:
            TypeError: If other is not a pydggrid instance
            ValueError: If other is None
        """
        if other is None:
            raise ValueError("Other cell cannot be None")
        if not isinstance(other, pydggrid):
            raise TypeError(
                f"Other cell must be a pydggrid instance, got {type(other)}"
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

    def get_dggrid_properties(self) -> Dict[str, Any]:
        """
        Get properties specific to the DGGRID cell.

        Returns:
            dict: Dictionary containing DGGRID-specific properties including:
                - level: Hierarchical level
                - grid_id: DGGRID identifier
                - resolution: Cell resolution
                - shape_type: Polygon shape description
                - edge_count: Number of edges/vertices
                - area: Cell area
                - characteristic_length: Effective cell size
        """
        properties = {
            "level": self.dggrid_level,
            "grid_id": self.dggrid_id,
            "resolution": self.dggrid_resolution,
            "shape_type": self.dggrid_shape_type,
            "edge_count": self.nEdge,
            "vertex_count": self.nVertex,
            "area": self.dArea,
            "characteristic_length": self.dLength,
        }
        return properties

    def get_hierarchical_neighbors(self, level_offset: int = 0) -> List[str]:
        """
        Get potential neighbor cell IDs at a different hierarchical level.

        Args:
            level_offset (int): Offset from current level (0 = same level, +1 = finer, -1 = coarser)

        Returns:
            List[str]: List of potential neighbor grid IDs (implementation dependent)

        Note:
            This is a placeholder method. Actual implementation would depend on
            the specific DGGRID system being used (ISEA, H3, etc.)
        """
        # This would need to be implemented based on the specific DGGRID system
        # For now, return empty list as placeholder
        target_level = self.dggrid_level + level_offset
        if target_level < 0:
            return []

        # Placeholder implementation
        return []

    def is_at_pole(self, tolerance: float = 0.01) -> bool:
        """
        Check if the DGGRID cell is near a pole.

        Args:
            tolerance (float): Tolerance in degrees for pole proximity

        Returns:
            bool: True if cell is within tolerance of north or south pole
        """
        return (
            abs(self.dLatitude_center_degree - 90.0) < tolerance
            or abs(self.dLatitude_center_degree + 90.0) < tolerance
        )

    def get_corner_coordinates(self) -> List[Tuple[float, float]]:
        """
        Get the coordinates of all corners of the DGGRID cell.

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
        Check if the DGGRID cell has valid attributes.

        A DGGRID cell is considered valid if:
        - It has valid center coordinates
        - It has at least 3 boundary vertices and edges
        - It has a positive area
        - It has a valid cell ID
        - All vertices and edges are properly defined
        - DGGRID-specific properties are consistent

        Returns:
            bool: True if cell has valid DGGRID attributes
        """
        # Check parent class validity
        if not super().is_valid():
            return False

        # Check DGGRID-specific requirements
        has_min_edges = len(self.aEdge) >= 3
        has_min_vertices = len(self.aVertex) >= 3
        has_positive_length = self.dLength > 0
        has_valid_level = self.dggrid_level >= -1  # -1 means unset

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
            has_min_edges
            and has_min_vertices
            and has_positive_length
            and has_valid_level
            and vertices_valid
        )

    def copy(self) -> "pydggrid":
        """
        Create a deep copy of the DGGRID cell.

        Returns:
            pydggrid: A new DGGRID cell object with the same attributes
        """
        # Create new DGGRID cell with same basic parameters
        new_dggrid = pydggrid(
            self.dLongitude_center_degree,
            self.dLatitude_center_degree,
            self.aEdge.copy(),
            self.aVertex.copy(),
        )

        # Copy all attributes from parent and this class
        for attr_name, attr_value in self.__dict__.items():
            if hasattr(new_dggrid, attr_name):
                if isinstance(attr_value, list) and attr_value is not None:
                    setattr(new_dggrid, attr_name, attr_value.copy())
                else:
                    setattr(new_dggrid, attr_name, attr_value)

        return new_dggrid

    def tojson(self) -> str:
        """
        Convert DGGRID cell object to a JSON string.

        Serializes all cell attributes including neighbor information and
        geometric properties. Uses the custom DggridClassEncoder to handle
        numpy data types and complex objects.

        Returns:
            str: JSON string representation of the DGGRID cell

        Example:
            >>> dggrid_cell = pydggrid(-77.0, 38.0, edges, vertices)
            >>> json_str = dggrid_cell.tojson()
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
            obj, sort_keys=True, indent=4, ensure_ascii=True, cls=DggridClassEncoder
        )
        return sJson
