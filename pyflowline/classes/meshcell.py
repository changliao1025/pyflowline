"""
PyMeshCell class for representing mesh cells in flowline networks.

This module defines mesh cell types and implements the pymeshcell class
that extends pypolygon to add flowline-specific mesh cell attributes
and functionality.
"""

import enum
import json
from json import JSONEncoder
from typing import Dict, Any, Optional, List
import numpy as np
from pyearth.toolbox.mesh.polygon import pypolygon
from pyearth.toolbox.mesh.point import pypoint
from pyearth.toolbox.mesh.line import pyline


from .vertex import pyvertex


class CellClassEncoder(JSONEncoder):
    """
    Custom JSON encoder for pymeshcell objects.

    Handles numpy data types and pyvertex objects, converting them to
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
        return JSONEncoder.default(self, obj)


class celltype(enum.Enum):
    """
    Enumeration of supported mesh cell types.

    Defines the different types of mesh cells that can be used in
    flowline network modeling and analysis.

    Attributes:
        hexagon (int): Hexagonal mesh cells (value: 1)
        square (int): Square/rectangular mesh cells (value: 2)
        latlon (int): Latitude-longitude grid cells (value: 3)
        mpas (int): MPAS (Model for Prediction Across Scales) cells (value: 4)
        dggrid (int): Discrete Global Grid cells (value: 5)
        tin (int): Triangulated Irregular Network cells (value: 6)
    """

    hexagon = 1
    square = 2
    latlon = 3
    mpas = 4
    dggrid = 5
    tin = 6

    def __str__(self) -> str:
        """Return string representation of cell type."""
        return self.name

    def __repr__(self) -> str:
        """Return detailed string representation of cell type."""
        return f"celltype.{self.name}"


class pymeshcell(pypolygon):
    """
    Mesh cell class for flowline network representation.

    Extends the pypolygon class to include mesh cell-specific attributes used
    in flowline topology and stream network analysis. A mesh cell represents
    a polygonal area that can contain flowlines and participate in flow
    routing calculations.

    Attributes:
        lCellID (int): Unique identifier for the mesh cell (default: -1)
        nFlowline (int): Number of flowlines intersecting this cell (default: 0)
        nPoint (int): Number of vertices defining the cell boundary
        nEdge (int): Number of edges defining the cell boundary
        dLength (float): Total length of cell perimeter in meters (default: 0.0)
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
        nNeighbor (int): Total number of neighboring cells (default: -1)
        nNeighbor_land (int): Number of land neighboring cells (default: -1)
        nNeighbor_ocean (int): Number of ocean neighboring cells (default: -1)
        nNeighbor_land_virtual (int): Number of virtual land neighbors (default: -1)
        aNeighbor_land_virtual (Optional[List]): List of virtual land neighbors
        aNeighbor (Optional[List]): List of all neighbor cell IDs
        aNeighbor_land (Optional[List]): List of land neighbor cell IDs
        aNeighbor_ocean (Optional[List]): List of ocean neighbor cell IDs
        aNeighbor_distance (Optional[List]): List of distances to neighbors
        pBound (Optional[Any]): Boundary object for spatial indexing
        aEdge (List): List of edges defining the cell boundary
        aVertex (List): List of vertices defining the cell boundary
        pPoint_center (pypoint): Center point of the cell

    Args:
        dLon (float): The longitude of the cell center in degrees
        dLat (float): The latitude of the cell center in degrees
        aEdge (List): A list of edges that define the cell boundary
        aVertex (List): A list of points that define the cell boundary

    Raises:
        ValueError: If cell has fewer than 3 edges (invalid polygon)
        TypeError: If input parameters are not of the expected types

    Example:
        >>> center_lon, center_lat = -77.0, 38.0
        >>> edges = [edge1, edge2, edge3, edge4]  # List of pyedge objects
        >>> points = [point1, point2, point3, point4, point1]  # Closed polygon
        >>> cell = pymeshcell(center_lon, center_lat, edges, points)
        >>> cell.set_cell_id(100)
        >>> print(cell)
        pymeshcell(ID=100, Center=(-77.0, 38.0), Area=1234.56m²)
    """

    def __init__(self, dLon: float, dLat: float, aEdge: List, aVertex: List) -> None:
        """
        Initialize a mesh cell object.

        Args:
            dLon (float): The longitude of the cell center in degrees
            dLat (float): The latitude of the cell center in degrees
            aEdge (List): A list of edges that define the cell boundary
            aVertex (List): A list of points that define the cell boundary

        Raises:
            ValueError: If cell has fewer than 3 edges or invalid coordinates
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

        # Validate minimum polygon requirements
        nEdge = len(aEdge)
        if nEdge < 3:
            raise ValueError(
                f"Cell must have at least 3 edges to form a valid polygon, got {nEdge}"
            )

        if len(aVertex) < 4:  # Need at least 3 unique points + 1 closing point
            raise ValueError(
                f"Cell must have at least 4 points (including closing point), got {len(aVertex)}"
            )

        # Convert aEdge to aLine with proper type handling
        aLine = []
        for i, edge in enumerate(aEdge):
            if edge is None:
                raise ValueError(f"Edge at index {i} cannot be None")
            # Convert to pyline using vertex attributes
            if hasattr(edge, "pVertex_start") and hasattr(edge, "pVertex_end"):
                point_start = pypoint(edge.pVertex_start.__dict__)
                point_end = pypoint(edge.pVertex_end.__dict__)
                line = pyline(point_start, point_end)
                aLine.append(line)
            else:
                raise ValueError(
                    f"Edge at index {i} must have pVertex_start and pVertex_end attributes"
                )

        # Convert aVertex to aPoint with proper type handling
        aPoint = []
        for i, vertex in enumerate(aVertex):
            if vertex is None:
                raise ValueError(f"Vertex at index {i} cannot be None")
            # Convert to pypoint using vertex attributes
            if hasattr(vertex, "__dict__"):
                point = pypoint(vertex.__dict__)
                aPoint.append(point)
            else:
                raise ValueError(f"Vertex at index {i} must be convertible to pypoint")

        # Initialize parent class with converted types
        # we need to add the closing point, this is different from the mesh vertex list
        aPoint.append(aPoint[0])
        super().__init__(dLon, dLat, aLine, aPoint)

        # Basic cell identification
        self.lCellID: int = -1

        # Flowline and geometry counts
        self.nFlowline: int = 0
        self.nVertex: int = len(aVertex)
        self.nEdge: int = len(aEdge)

        # Geometric properties
        self.dLength: float = 0.0  # Perimeter length
        self.dArea: float = 0.0
        self.dX_center_meter: float = 0.0
        self.dY_center_meter: float = 0.0
        self.dz_center: float = 0.0
        self.dLongitude_center_degree: float = float(dLon)
        self.dLatitude_center_degree: float = float(dLat)

        # Elevation properties
        self.dElevation_mean: float = -9999.0
        self.dElevation_profile0: float = 0.0

        # Flowline properties
        self.dLength_flowline: float = 0.0

        # Status flags
        self.iFlag_intersected: int = -1
        self.iFlag_coast: int = 0

        # Downstream connectivity (after burning)
        self.lCellID_downstream_burned: int = -1
        self.iStream_order_burned: int = -1
        self.iStream_segment_burned: int = -1

        # Neighbor information
        self.nNeighbor: int = -1
        self.nNeighbor_land: int = -1
        self.nNeighbor_ocean: int = -1
        self.nNeighbor_land_virtual: int = -1

        # Neighbor lists (initialized as None, will be populated as needed)
        self.aNeighbor_land_virtual: Optional[List] = None
        self.aNeighbor: Optional[List] = None  # Global IDs of all neighbors
        self.aNeighbor_land: Optional[List] = None  # Global IDs of land neighbors
        self.aNeighbor_ocean: Optional[List] = None  # Global IDs of ocean neighbors
        self.aNeighbor_distance: Optional[List] = None

        # Spatial indexing
        self.pBound: Optional[Any] = None

        # Store geometry
        self.aEdge: List = aEdge
        self.aVertex: List = aVertex  # First and last points should be the same

        # Create center point object
        center_params = {
            "dLongitude_degree": self.dLongitude_center_degree,
            "dLatitude_degree": self.dLatitude_center_degree,
        }
        self.pVertex_center: pyvertex = pyvertex(center_params)

        # Calculate cell boundary for spatial indexing
        self.calculate_cell_bound()

    def __repr__(self) -> str:
        """
        Return a detailed string representation of the mesh cell.

        Returns:
            str: Detailed representation including cell ID, center coordinates, and area
        """
        return (
            f"pymeshcell(ID={self.lCellID}, "
            f"Center=({self.dLongitude_center_degree:.6f}, {self.dLatitude_center_degree:.6f}), "
            f"Area={self.dArea:.2f}m², Flowlines={self.nFlowline}, "
            f"Neighbors={self.nNeighbor})"
        )

    def __str__(self) -> str:
        """
        Return a concise string representation of the mesh cell.

        Returns:
            str: Concise representation with ID, center coordinates, and area
        """
        return (
            f"pymeshcell(ID={self.lCellID}, "
            f"Center=({self.dLongitude_center_degree:.6f}, {self.dLatitude_center_degree:.6f}), "
            f"Area={self.dArea:.2f}m²)"
        )

    def __hash__(self) -> int:
        """
        Return hash value for the mesh cell.

        Uses parent class hash which is based on the polygon geometry.
        This allows cells to be used in sets and as dictionary keys.

        Returns:
            int: Hash value based on geometry
        """
        return super().__hash__()

    def __eq__(self, other: Any) -> bool:
        """
        Check if two mesh cells are equal based on their geometry.

        Cells are considered equal if their boundary geometry matches
        within a precision threshold, regardless of their IDs or other attributes.

        Args:
            other: Another object to compare with

        Returns:
            bool: True if cells have the same geometry
        """
        if not isinstance(other, pymeshcell):
            return NotImplemented
        return super().__eq__(other)

    def tojson(self) -> str:
        """
        Convert mesh cell object to a JSON string.

        Serializes all cell attributes including neighbor information and
        geometric properties. Uses the custom CellClassEncoder to handle
        numpy data types and complex objects.

        Returns:
            str: JSON string representation of the mesh cell

        Example:
            >>> cell = pymeshcell(-77.0, 38.0, edges, points)
            >>> json_str = cell.tojson()
        """
        aSkip = [
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
            obj, sort_keys=True, indent=4, ensure_ascii=True, cls=CellClassEncoder
        )
        return sJson

    def set_cell_id(self, lCellID: int) -> None:
        """
        Set the mesh cell ID.

        Args:
            lCellID (int): New cell ID

        Raises:
            TypeError: If lCellID is not an integer
        """
        if not isinstance(lCellID, (int, np.integer)):
            raise TypeError(f"Cell ID must be an integer, got {type(lCellID)}")
        self.lCellID = int(lCellID)

    def set_downstream_cell(self, lCellID_downstream: int) -> None:
        """
        Set the downstream cell ID after flow routing.

        Args:
            lCellID_downstream (int): ID of the downstream cell

        Raises:
            TypeError: If lCellID_downstream is not an integer
        """
        if not isinstance(lCellID_downstream, (int, np.integer)):
            raise TypeError(
                f"Downstream cell ID must be an integer, got {type(lCellID_downstream)}"
            )
        self.lCellID_downstream_burned = int(lCellID_downstream)

    def set_stream_properties(self, iStream_order: int, iStream_segment: int) -> None:
        """
        Set stream order and segment after flow processing.

        Args:
            iStream_order (int): Stream order value
            iStream_segment (int): Stream segment ID

        Raises:
            TypeError: If parameters are not integers
        """
        if not isinstance(iStream_order, (int, np.integer)):
            raise TypeError(
                f"Stream order must be an integer, got {type(iStream_order)}"
            )
        if not isinstance(iStream_segment, (int, np.integer)):
            raise TypeError(
                f"Stream segment must be an integer, got {type(iStream_segment)}"
            )

        self.iStream_order_burned = int(iStream_order)
        self.iStream_segment_burned = int(iStream_segment)

    def add_neighbor(
        self,
        lNeighborID: int,
        neighbor_type: str = "land",
        distance: Optional[float] = None,
    ) -> None:
        """
        Add a neighboring cell to the appropriate neighbor list.

        Args:
            lNeighborID (int): ID of the neighboring cell
            neighbor_type (str): Type of neighbor ('land', 'ocean', 'virtual')
            distance (Optional[float]): Distance to the neighbor cell

        Raises:
            ValueError: If neighbor_type is not recognized
            TypeError: If lNeighborID is not an integer
        """
        if not isinstance(lNeighborID, (int, np.integer)):
            raise TypeError(f"Neighbor ID must be an integer, got {type(lNeighborID)}")

        valid_types = ["land", "ocean", "virtual"]
        if neighbor_type not in valid_types:
            raise ValueError(
                f"neighbor_type must be one of {valid_types}, got {neighbor_type}"
            )

        # Initialize neighbor lists if they don't exist
        if self.aNeighbor is None:
            self.aNeighbor = []
        if self.aNeighbor_distance is None:
            self.aNeighbor_distance = []

        # Add to main neighbor list
        if lNeighborID not in self.aNeighbor:
            self.aNeighbor.append(lNeighborID)
            self.aNeighbor_distance.append(distance if distance is not None else 0.0)
            self.nNeighbor = len(self.aNeighbor)

        # Add to specific type lists
        if neighbor_type == "land":
            if self.aNeighbor_land is None:
                self.aNeighbor_land = []
            if lNeighborID not in self.aNeighbor_land:
                self.aNeighbor_land.append(lNeighborID)
                self.nNeighbor_land = len(self.aNeighbor_land)

        elif neighbor_type == "ocean":
            if self.aNeighbor_ocean is None:
                self.aNeighbor_ocean = []
            if lNeighborID not in self.aNeighbor_ocean:
                self.aNeighbor_ocean.append(lNeighborID)
                self.nNeighbor_ocean = len(self.aNeighbor_ocean)

        elif neighbor_type == "virtual":
            if self.aNeighbor_land_virtual is None:
                self.aNeighbor_land_virtual = []
            if lNeighborID not in self.aNeighbor_land_virtual:
                self.aNeighbor_land_virtual.append(lNeighborID)
                self.nNeighbor_land_virtual = len(self.aNeighbor_land_virtual)

    def is_valid(self) -> bool:
        """
        Check if the mesh cell has valid attributes.

        A cell is considered valid if:
        - It has valid center coordinates
        - It has at least 3 boundary points
        - It has a positive area (if calculated)
        - It has a valid cell ID

        Returns:
            bool: True if cell has valid attributes
        """
        has_valid_coords = (
            -180 <= self.dLongitude_center_degree <= 180
            and -90 <= self.dLatitude_center_degree <= 90
        )
        has_valid_geometry = len(self.aVertex) >= 4 and len(self.aEdge) >= 3
        has_valid_id = self.lCellID > 0
        has_valid_area = self.dArea >= 0  # Area should be non-negative

        return (
            has_valid_coords and has_valid_geometry and has_valid_id and has_valid_area
        )

    def copy(self) -> "pymeshcell":
        """
        Create a deep copy of the mesh cell.

        Returns:
            pymeshcell: A new mesh cell object with the same attributes
        """
        # Create new cell with same basic parameters
        new_cell = pymeshcell(
            self.dLongitude_center_degree,
            self.dLatitude_center_degree,
            self.aEdge.copy(),
            self.aVertex.copy(),
        )

        # Copy all attributes
        new_cell.lCellID = self.lCellID
        new_cell.nFlowline = self.nFlowline
        new_cell.dLength = self.dLength
        new_cell.dArea = self.dArea
        new_cell.dX_center_meter = self.dX_center_meter
        new_cell.dY_center_meter = self.dY_center_meter
        new_cell.dz_center = self.dz_center
        new_cell.dElevation_mean = self.dElevation_mean
        new_cell.dElevation_profile0 = self.dElevation_profile0
        new_cell.dLength_flowline = self.dLength_flowline
        new_cell.iFlag_intersected = self.iFlag_intersected
        new_cell.iFlag_coast = self.iFlag_coast
        new_cell.lCellID_downstream_burned = self.lCellID_downstream_burned
        new_cell.iStream_order_burned = self.iStream_order_burned
        new_cell.iStream_segment_burned = self.iStream_segment_burned

        # Copy neighbor information
        new_cell.nNeighbor = self.nNeighbor
        new_cell.nNeighbor_land = self.nNeighbor_land
        new_cell.nNeighbor_ocean = self.nNeighbor_ocean
        new_cell.nNeighbor_land_virtual = self.nNeighbor_land_virtual

        # Deep copy neighbor lists
        new_cell.aNeighbor = (
            self.aNeighbor.copy() if self.aNeighbor is not None else None
        )
        new_cell.aNeighbor_land = (
            self.aNeighbor_land.copy() if self.aNeighbor_land is not None else None
        )
        new_cell.aNeighbor_ocean = (
            self.aNeighbor_ocean.copy() if self.aNeighbor_ocean is not None else None
        )
        new_cell.aNeighbor_land_virtual = (
            self.aNeighbor_land_virtual.copy()
            if self.aNeighbor_land_virtual is not None
            else None
        )
        new_cell.aNeighbor_distance = (
            self.aNeighbor_distance.copy()
            if self.aNeighbor_distance is not None
            else None
        )

        return new_cell

    def is_coastal(self) -> bool:
        """
        Check if the mesh cell is coastal (has ocean neighbors).

        Returns:
            bool: True if cell has ocean neighbors or is flagged as coastal
        """
        return self.iFlag_coast == 1 or (
            self.nNeighbor_ocean is not None and self.nNeighbor_ocean > 0
        )

    def has_flowlines(self) -> bool:
        """
        Check if the mesh cell contains flowlines.

        Returns:
            bool: True if cell has intersecting flowlines
        """
        return self.nFlowline > 0 and self.iFlag_intersected > 0

    def get_neighbor_count(self, neighbor_type: str = "all") -> int:
        """
        Get the count of neighbors of a specific type.

        Args:
            neighbor_type (str): Type of neighbors to count ('all', 'land', 'ocean', 'virtual')

        Returns:
            int: Number of neighbors of the specified type

        Raises:
            ValueError: If neighbor_type is not recognized
        """
        valid_types = ["all", "land", "ocean", "virtual"]
        if neighbor_type not in valid_types:
            raise ValueError(
                f"neighbor_type must be one of {valid_types}, got {neighbor_type}"
            )

        if neighbor_type == "all":
            return self.nNeighbor if self.nNeighbor >= 0 else 0
        elif neighbor_type == "land":
            return self.nNeighbor_land if self.nNeighbor_land >= 0 else 0
        elif neighbor_type == "ocean":
            return self.nNeighbor_ocean if self.nNeighbor_ocean >= 0 else 0
        elif neighbor_type == "virtual":
            return (
                self.nNeighbor_land_virtual if self.nNeighbor_land_virtual >= 0 else 0
            )

    def get_center_coordinates(self) -> tuple:
        """
        Get the center coordinates of the mesh cell.

        Returns:
            tuple: (longitude, latitude) of the cell center
        """
        return (self.dLongitude_center_degree, self.dLatitude_center_degree)

    def update_flowline_stats(self, nFlowlines: int, dLength_total: float) -> None:
        """
        Update flowline statistics for the mesh cell.

        Args:
            nFlowlines (int): Number of flowlines intersecting the cell
            dLength_total (float): Total length of flowlines in the cell

        Raises:
            TypeError: If parameters are not of the expected types
            ValueError: If values are negative
        """
        if not isinstance(nFlowlines, (int, np.integer)):
            raise TypeError(f"nFlowlines must be an integer, got {type(nFlowlines)}")
        if not isinstance(dLength_total, (int, float, np.number)):
            raise TypeError(
                f"dLength_total must be a number, got {type(dLength_total)}"
            )

        if nFlowlines < 0:
            raise ValueError(f"nFlowlines must be non-negative, got {nFlowlines}")
        if dLength_total < 0:
            raise ValueError(f"dLength_total must be non-negative, got {dLength_total}")

        self.nFlowline = int(nFlowlines)
        self.dLength_flowline = float(dLength_total)
        self.iFlag_intersected = 1 if nFlowlines > 0 else 0
