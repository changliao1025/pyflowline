"""
PyLatlon class for representing latitude-longitude grid cells in flowline networks.

This module extends the pymeshcell class from pyflowline to add latitude-longitude
specific cell attributes and functionality. A lat-lon cell represents a rectangular
grid cell defined by longitude and latitude boundaries.
"""

import json
from json import JSONEncoder
from typing import List, Dict, Any, Optional, Union, Tuple
import numpy as np
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area
from pyearthriver.classes.vertex import pyvertex
from pyearthriver.classes.edge import pyedge
from pyearthriver.classes.flowline import pyflowline
from pyflowline.classes.meshcell import pymeshcell


class LatlonClassEncoder(JSONEncoder):
    """
    Custom JSON encoder for pylatlon objects.

    Handles numpy data types, pyvertex, and pyedge objects, converting them to
    native Python types for JSON serialization.
    """

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.float32):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, list):
            return obj
        if isinstance(obj, pyvertex):
            return json.loads(obj.tojson())
        if isinstance(obj, pyedge):
            return obj.lEdgeID
        if isinstance(obj, pyflowline):
            return obj.lFlowlineID
        if isinstance(obj, pylatlon):
            return obj.lCellID

        return JSONEncoder.default(self, obj)


class pylatlon(pymeshcell):
    """
    Latitude-longitude grid cell class for flowline network representation.

    Extends the pymeshcell class to represent rectangular grid cells defined by
    longitude and latitude boundaries. A lat-lon cell is typically a quadrilateral
    with four edges and four vertices forming a rectangular grid element.

    Attributes:
        lCellID (int): Unique identifier for the lat-lon cell (default: -1)
        nFlowline (int): Number of flowlines intersecting this cell (default: 0)
        nVertex (int): Number of vertices (should be 4 for lat-lon cells)
        nEdge (int): Number of edges (should be 4 for lat-lon cells)
        dLength (float): Effective length of the cell in meters (default: 0.0)
        dArea (float): Area of the cell in square meters (default: 0.0)
        dX_center_meter (float): X coordinate of center in meters (default: 0.0)
        dY_center_meter (float): Y coordinate of center in meters (default: 0.0)
        dz_center (float): Z coordinate of center in meters (default: 0.0)
        dLongitude_center_degree (float): Longitude of cell center in degrees
        dLatitude_center_degree (float): Latitude of cell center in degrees
        dElevation_mean (float): Mean elevation of the cell (default: -9999.0)
        dElevation_profile0 (float): Elevation profile parameter (default: 0.0)
        dLength_flowline (float): Total length of flowlines in cell (default: 0.0)
        iFlag_intersected (int): Flag indicating intersection status (default: -1)
        iFlag_coast (int): Flag indicating if cell is coastal (default: 0)
        lCellID_downstream_burned (int): ID of downstream cell (default: -1)
        iStream_order_burned (int): Stream order after burning (default: -1)
        iStream_segment_burned (int): Stream segment ID after burning (default: -1)
        aEdge (List[pyedge]): List of edges defining the cell boundary
        aEdgeID (Optional[List[int]]): List of edge IDs
        aVertex (List[pyvertex]): List of vertices defining the cell boundary
        aVertexID (Optional[List[int]]): List of vertex IDs
        pVertex_center (pyvertex): Center vertex of the cell
        aFlowline (Optional[List]): List of flowlines intersecting the cell
        nNeighbor (int): Total number of neighboring cells (default: -1)
        nNeighbor_land (int): Number of land neighboring cells (default: -1)
        nNeighbor_ocean (int): Number of ocean neighboring cells (default: -1)
        nNeighbor_land_virtual (int): Number of virtual land neighbors (default: -1)
        aNeighbor_land_virtual (Optional[List]): List of virtual land neighbors
        aNeighbor (Optional[List[int]]): List of all neighbor cell IDs
        aNeighbor_land (Optional[List[int]]): List of land neighbor cell IDs
        aNeighbor_ocean (Optional[List[int]]): List of ocean neighbor cell IDs
        aNeighbor_distance (Optional[List[float]]): List of distances to neighbors
        pBound (Optional[Tuple]): Boundary tuple for spatial indexing

    Args:
        dLon (float): The longitude of the cell center in degrees
        dLat (float): The latitude of the cell center in degrees
        aEdge (List[pyedge]): A list of 4 edges that define the lat-lon cell
        aVertex (List[pyvertex]): A list of 4 vertices that define the lat-lon cell

    Raises:
        TypeError: If input parameters are not of the expected types
        ValueError: If the cell doesn't have exactly 4 edges and vertices
        ValueError: If coordinates are outside valid ranges

    Example:
        >>> # Create vertices for a lat-lon cell
        >>> v1 = pyvertex({'dLongitude_degree': -77.0, 'dLatitude_degree': 38.0})
        >>> v2 = pyvertex({'dLongitude_degree': -76.9, 'dLatitude_degree': 38.0})
        >>> v3 = pyvertex({'dLongitude_degree': -76.9, 'dLatitude_degree': 38.1})
        >>> v4 = pyvertex({'dLongitude_degree': -77.0, 'dLatitude_degree': 38.1})
        >>> # Create edges
        >>> e1 = pyedge(v1, v2)
        >>> e2 = pyedge(v2, v3)
        >>> e3 = pyedge(v3, v4)
        >>> e4 = pyedge(v4, v1)
        >>> # Create lat-lon cell
        >>> cell = pylatlon(-76.95, 38.05, [e1, e2, e3, e4], [v1, v2, v3, v4, v1])
        >>> cell.set_cell_id(100)
        >>> print(cell)
        pylatlon(ID=100, Center=(-76.95, 38.05), Area=123.45m²)
    """

    def __init__(
        self, dLon: float, dLat: float, aEdge: List[pyedge], aVertex: List[pyvertex]
    ) -> None:
        """
        Initialize a lat-lon cell object.

        Args:
            dLon (float): The longitude of the cell center in degrees (-180 to 180)
            dLat (float): The latitude of the cell center in degrees (-90 to 90)
            aEdge (List[pyedge]): A list of exactly 4 edges that define the cell
            aVertex (List[pyvertex]): A list of vertices that define the cell (first and last should be the same)

        Raises:
            TypeError: If input parameters are not of the expected types
            ValueError: If the cell doesn't have exactly 4 edges or valid coordinates
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
            raise ValueError(
                f"Longitude must be between -180 and 180 degrees, got {dLon}"
            )
        if not (-90 <= dLat <= 90):
            raise ValueError(f"Latitude must be between -90 and 90 degrees, got {dLat}")

        # Validate lat-lon cell requirements (must be rectangular)
        if len(aEdge) != 4:
            raise ValueError(
                f"Lat-lon cell must have exactly 4 edges, got {len(aEdge)}"
            )
        if len(aVertex) < 4:
            raise ValueError(
                f"Lat-lon cell must have at least 4 vertices, got {len(aVertex)}"
            )

        # Validate edges are not None
        for i, edge in enumerate(aEdge):
            if edge is None:
                raise ValueError(f"Edge at index {i} cannot be None")
            if not hasattr(edge, "pVertex_start") or not hasattr(edge, "pVertex_end"):
                raise ValueError(f"Edge at index {i} must be a valid pyedge object")

        # Validate vertices are not None
        for i, vertex in enumerate(aVertex):
            if vertex is None:
                raise ValueError(f"Vertex at index {i} cannot be None")

        # Initialize parent class (pymeshcell)
        super().__init__(dLon, dLat, aEdge, aVertex)

        # Override parent initialization with lat-lon specific values
        self.nEdge: int = 4
        self.nVertex: int = 4
        self.dLongitude_center_degree: float = float(dLon)
        self.dLatitude_center_degree: float = float(dLat)

        # Set geometry
        self.aEdge: List[pyedge] = aEdge
        self.aVertex: List[pyvertex] = aVertex

        # Initialize lat-lon specific attributes
        self.aEdgeID: Optional[List[int]] = None
        self.aVertexID: Optional[List[int]] = None

        # Create center vertex
        pVertex_params = {
            "dLongitude_degree": self.dLongitude_center_degree,
            "dLatitude_degree": self.dLatitude_center_degree,
        }
        self.pVertex_center: pyvertex = pyvertex(pVertex_params)

        # Initialize with default values for stream processing
        self.dElevation_mean: float = -9999.0

        # Calculate cell boundary for spatial indexing
        self.calculate_cell_bound()

    def __repr__(self) -> str:
        """
        Return a detailed string representation of the lat-lon cell.

        Returns:
            str: Detailed representation including cell ID, center coordinates, and area
        """
        return (
            f"pylatlon(ID={self.lCellID}, "
            f"Center=({self.dLongitude_center_degree:.6f}, {self.dLatitude_center_degree:.6f}), "
            f"Area={self.dArea:.2f}m², Flowlines={self.nFlowline}, "
            f"StreamOrder={self.iStream_order_burned})"
        )

    def __str__(self) -> str:
        """
        Return a concise string representation of the lat-lon cell.

        Returns:
            str: Concise representation with ID, center coordinates, and area
        """
        return (
            f"pylatlon(ID={self.lCellID}, "
            f"Center=({self.dLongitude_center_degree:.6f}, {self.dLatitude_center_degree:.6f}), "
            f"Area={self.dArea:.2f}m²)"
        )

    def __eq__(self, other: Any) -> bool:
        """
        Check if two lat-lon cells are equal based on their center coordinates and geometry.

        Args:
            other: Another object to compare with

        Returns:
            bool: True if cells have the same center coordinates and geometry
        """
        if not isinstance(other, pylatlon):
            return NotImplemented

        # Check center coordinates
        coord_match = (
            abs(self.dLongitude_center_degree - other.dLongitude_center_degree) < 1e-9
            and abs(self.dLatitude_center_degree - other.dLatitude_center_degree) < 1e-9
        )

        # Use parent class geometry comparison
        return coord_match and super().__eq__(other)

    def __hash__(self) -> int:
        """
        Return hash value for the lat-lon cell.

        Returns:
            int: Hash value based on center coordinates
        """
        return hash(
            (
                round(self.dLongitude_center_degree, 9),
                round(self.dLatitude_center_degree, 9),
            )
        )

    def calculate_cell_bound(self) -> Tuple[float, float, float, float]:
        """
        Calculate the bounding box of the lat-lon cell.

        Returns:
            Tuple[float, float, float, float]: Bounding box as (lon_min, lat_min, lon_max, lat_max)

        Raises:
            ValueError: If cell has no vertices or invalid geometry
        """
        if not self.aVertex or self.nVertex == 0:
            raise ValueError("Cannot calculate bounds: cell has no vertices")

        dLat_min = 90.0
        dLat_max = -90.0
        dLon_min = 180.0
        dLon_max = -180.0

        for i in range(self.nVertex):
            vertex = self.aVertex[i]
            if not hasattr(vertex, "dLongitude_degree") or not hasattr(
                vertex, "dLatitude_degree"
            ):
                raise ValueError(f"Vertex at index {i} missing coordinate attributes")

            dLon_max = max(dLon_max, vertex.dLongitude_degree)
            dLon_min = min(dLon_min, vertex.dLongitude_degree)
            dLat_max = max(dLat_max, vertex.dLatitude_degree)
            dLat_min = min(dLat_min, vertex.dLatitude_degree)

        self.pBound = (dLon_min, dLat_min, dLon_max, dLat_max)
        return self.pBound

    def has_this_edge(self, pEdge_in: pyedge) -> bool:
        """
        Check whether the lat-lon cell contains a specific edge.

        Args:
            pEdge_in (pyedge): The edge to be checked

        Returns:
            bool: True if the edge is found in the cell, False otherwise

        Raises:
            TypeError: If pEdge_in is not a pyedge object
            ValueError: If the cell has no edges
        """
        if not isinstance(pEdge_in, pyedge):
            raise TypeError(f"Expected pyedge object, got {type(pEdge_in)}")

        if not self.aEdge:
            raise ValueError("Cell has no edges to check against")

        for pEdge in self.aEdge:
            if hasattr(pEdge, "is_overlap") and pEdge.is_overlap(pEdge_in):
                return True

        return False

    def which_edge_crosses_vertex(
        self, pVertex_in: pyvertex
    ) -> Tuple[bool, Optional[pyedge]]:
        """
        Find which edge crosses or contains a specific vertex.

        Args:
            pVertex_in (pyvertex): The vertex to be checked

        Returns:
            Tuple[bool, Optional[pyedge]]: (True, edge) if found, (False, None) if not found

        Raises:
            TypeError: If pVertex_in is not a pyvertex object
            ValueError: If the cell has no edges
        """
        if not isinstance(pVertex_in, pyvertex):
            raise TypeError(f"Expected pyvertex object, got {type(pVertex_in)}")

        if not self.aEdge:
            raise ValueError("Cell has no edges to check against")

        for pEdge in self.aEdge:
            if hasattr(pEdge, "check_vertex_on_edge"):
                try:
                    iFlag, _, _ = pEdge.check_vertex_on_edge(pVertex_in)
                    if iFlag == 1:
                        return True, pEdge
                except Exception:
                    # Skip edges that don't support this operation
                    continue

        return False, None

    def calculate_cell_area(self) -> float:
        """
        Calculate the area of the lat-lon cell using spherical geometry.

        Returns:
            float: The area in square meters

        Raises:
            ValueError: If the cell has insufficient vertices for area calculation
        """
        if not self.aVertex or self.nVertex < 3:
            raise ValueError(
                "Cannot calculate area: cell must have at least 3 vertices"
            )

        try:
            lons = [vertex.dLongitude_degree for vertex in self.aVertex[: self.nVertex]]
            lats = [vertex.dLatitude_degree for vertex in self.aVertex[: self.nVertex]]

            self.dArea = calculate_polygon_area(lons, lats)
            return self.dArea

        except Exception as e:
            raise ValueError(f"Failed to calculate cell area: {e}")

    def calculate_effective_length(self) -> float:
        """
        Calculate the effective length of the lat-lon cell.

        The effective length is the square root of the cell area, providing
        a characteristic length scale for the cell.

        Returns:
            float: The effective length in meters

        Raises:
            ValueError: If area has not been calculated or is negative
        """
        if self.dArea <= 0:
            self.calculate_cell_area()

        if self.dArea <= 0:
            raise ValueError(
                "Cannot calculate effective length: cell area is zero or negative"
            )

        self.dLength = np.sqrt(self.dArea)
        return self.dLength

    def shares_edge_with(self, other: "pylatlon") -> bool:
        """
        Check whether this lat-lon cell shares an edge with another lat-lon cell.

        Args:
            other (pylatlon): The other lat-lon cell to check against

        Returns:
            bool: True if cells share an edge, False otherwise

        Raises:
            TypeError: If other is not a pylatlon object
            ValueError: If either cell has no edges
        """
        if not isinstance(other, pylatlon):
            raise TypeError(f"Expected pylatlon object, got {type(other)}")

        if not self.aEdge or not other.aEdge:
            raise ValueError("Both cells must have edges to check for sharing")

        for pEdge1 in self.aEdge:
            for pEdge2 in other.aEdge:
                if hasattr(pEdge1, "is_overlap") and pEdge1.is_overlap(pEdge2):
                    return True

        return False

    def is_valid(self) -> bool:
        """
        Check if the lat-lon cell has valid attributes and geometry.

        Returns:
            bool: True if cell is valid, False otherwise
        """
        # Check basic validity from parent class
        if not super().is_valid():
            return False

        # Check lat-lon specific requirements
        if self.nEdge != 4 or self.nVertex != 4:
            return False

        # Check coordinate validity
        if not (-180 <= self.dLongitude_center_degree <= 180):
            return False
        if not (-90 <= self.dLatitude_center_degree <= 90):
            return False

        # Check that we have required geometry
        if not self.aEdge or len(self.aEdge) != 4:
            return False
        if not self.aVertex or len(self.aVertex) < 4:
            return False

        return True

    def get_center_coordinates(self) -> Tuple[float, float]:
        """
        Get the center coordinates of the lat-lon cell.

        Returns:
            Tuple[float, float]: (longitude, latitude) of the cell center in degrees
        """
        return (self.dLongitude_center_degree, self.dLatitude_center_degree)

    def get_bounds(self) -> Tuple[float, float, float, float]:
        """
        Get the bounding box of the lat-lon cell.

        Returns:
            Tuple[float, float, float, float]: Bounding box as (lon_min, lat_min, lon_max, lat_max)
        """
        if self.pBound is None:
            self.calculate_cell_bound()
        return self.pBound

    def set_cell_id(self, lCellID: int) -> None:
        """
        Set the lat-lon cell ID.

        Args:
            lCellID (int): New cell ID

        Raises:
            TypeError: If lCellID is not an integer
        """
        if not isinstance(lCellID, (int, np.integer)):
            raise TypeError(f"Cell ID must be an integer, got {type(lCellID)}")
        self.lCellID = int(lCellID)

    def copy(self) -> "pylatlon":
        """
        Create a deep copy of the lat-lon cell.

        Returns:
            pylatlon: A new lat-lon cell object with the same attributes
        """
        # Create copies of edges and vertices
        aEdge_copy = [
            edge.copy() if hasattr(edge, "copy") else edge for edge in self.aEdge
        ]
        aVertex_copy = [
            vertex.copy() if hasattr(vertex, "copy") else vertex
            for vertex in self.aVertex
        ]

        # Create new cell
        new_cell = pylatlon(
            self.dLongitude_center_degree,
            self.dLatitude_center_degree,
            aEdge_copy,
            aVertex_copy,
        )

        # Copy all attributes from parent
        if hasattr(super(), "copy"):
            parent_copy = super().copy()
            for attr in parent_copy.__dict__:
                if hasattr(new_cell, attr):
                    setattr(new_cell, attr, getattr(parent_copy, attr))

        return new_cell

    def tojson(self) -> str:
        """
        Convert lat-lon cell object to a JSON string.

        Serializes most cell attributes while excluding complex objects
        that are better handled separately.

        Returns:
            str: JSON string representation of the lat-lon cell

        Example:
            >>> cell = pylatlon(-77.0, 38.0, edges, vertices)
            >>> json_str = cell.tojson()
        """
        aSkip = ["aEdge", "aFlowline", "dLongitude_radian", "dLatitude_radian", "wkt"]

        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)

        sJson = json.dumps(
            obj, sort_keys=True, indent=4, ensure_ascii=True, cls=LatlonClassEncoder
        )
        return sJson
