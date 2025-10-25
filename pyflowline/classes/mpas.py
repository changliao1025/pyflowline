
"""
PyMPAS class for representing MPAS mesh cells in flowline networks.

This module defines the pympas class which extends pymeshcell to provide
MPAS-specific mesh cell functionality for flowline network modeling.
MPAS (Model for Prediction Across Scales) uses unstructured meshes with
variable polygon shapes.
"""

import json
from json import JSONEncoder
from typing import List, Optional, Any, Tuple
import numpy as np

from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area

from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.edge import pyedge
from pyflowline.classes.flowline import pyflowline
from pyflowline.classes.meshcell import pymeshcell

try:
    from .flowline import pyflowline
except ImportError:
    pyflowline = None


class MpasClassEncoder(JSONEncoder):
    """
    Custom JSON encoder for pympas objects.

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
            return getattr(obj, 'lEdgeID', str(obj))
        if pyflowline and isinstance(obj, pyflowline):
            return getattr(obj, 'lFlowlineID', str(obj))
        if isinstance(obj, pympas):
            return obj.lCellID
        return JSONEncoder.default(self, obj)


class pympas(pymeshcell):
    """
    MPAS mesh cell class for flowline network representation.

    Extends the pymeshcell class to provide MPAS-specific functionality for
    mesh cells in flowline topology and stream network analysis. An MPAS cell
    represents a polygonal area with variable number of vertices and edges (3-10).

    MPAS (Model for Prediction Across Scales) uses unstructured meshes that can
    adapt to varying geometric requirements. This flexibility makes MPAS cells
    suitable for complex geographical features and multi-resolution modeling.

    Attributes:
        All attributes from pymeshcell plus MPAS-specific properties.
        The MPAS cell can have 3-10 edges and vertices, allowing for triangular
        to decagonal shapes depending on the mesh requirements.
        iFlag_watershed_boundary_burned (int): Flag for watershed boundary processing

    Args:
        dLon (float): The longitude of the MPAS cell center in degrees
        dLat (float): The latitude of the MPAS cell center in degrees
        aEdge (List): A list of 3-10 edges that define the MPAS cell boundary
        aVertex (List): A list of 3-10 vertices that define the MPAS cell boundary

    Raises:
        ValueError: If the cell doesn't have 3-10 edges and vertices
        TypeError: If input parameters are not of the expected types

    Example:
        >>> center_lon, center_lat = -77.0, 38.0
        >>> edges = [edge1, edge2, edge3, edge4, edge5]  # List of 5 edges
        >>> vertices = [v1, v2, v3, v4, v5]  # List of 5 vertices
        >>> mpas_cell = pympas(center_lon, center_lat, edges, vertices)
        >>> mpas_cell.set_cell_id(100)
        >>> print(mpas_cell)
        pympas(ID=100, Center=(-77.0, 38.0), Edges=5, Area=1234.56m²)
    """

    def __init__(self, dLon: float, dLat: float, aEdge: List, aVertex: List) -> None:
        """
        Initialize an MPAS mesh cell object.

        Args:
            dLon (float): The longitude of the MPAS cell center in degrees
            dLat (float): The latitude of the MPAS cell center in degrees
            aEdge (List): A list of 3-10 edges that define the MPAS cell boundary
            aVertex (List): A list of 3-10 vertices that define the MPAS cell boundary

        Raises:
            ValueError: If the cell doesn't have 3-10 edges and vertices
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

        # Validate MPAS geometry requirements (3-10 edges for flexibility)
        nEdge = len(aEdge)
        nVertex = len(aVertex)

        if not (3 <= nEdge <= 10):
            raise ValueError(f"MPAS cell must have 3-10 edges, got {nEdge}")
        if not (3 <= nVertex <= 10):
            raise ValueError(f"MPAS cell must have 3-10 vertices, got {nVertex}")
        if nEdge != nVertex:
            raise ValueError(f"Number of edges ({nEdge}) must equal number of vertices ({nVertex})")

        # Validate edges and vertices are not None
        for i, edge in enumerate(aEdge):
            if edge is None:
                raise ValueError(f"Edge at index {i} cannot be None")

        for i, vertex in enumerate(aVertex):
            if vertex is None:
                raise ValueError(f"Vertex at index {i} cannot be None")

        # Initialize parent class
        super().__init__(dLon, dLat, aEdge, aVertex)

        # MPAS-specific initialization
        self.nEdge = nEdge
        self.nVertex = nVertex  # Note: Different from nPoint which includes closing vertex

        # MPAS-specific attributes
        self.iFlag_watershed_boundary_burned: int = 0
        self.dElevation_profile0: float = -9999.0

        # Create center vertex if pyvertex is available
        if pyvertex:
            pVertex_params = {
                'dLongitude_degree': self.dLongitude_center_degree,
                'dLatitude_degree': self.dLatitude_center_degree
            }
            self.pVertex_center = pyvertex(pVertex_params)
        else:
            self.pVertex_center = None

        # Calculate initial properties
        self.calculate_cell_area()
        if self.dArea > 0:
            self.calculate_edge_length()

    def __repr__(self) -> str:
        """
        Return a detailed string representation of the MPAS cell.

        Returns:
            str: Detailed representation including cell ID, center coordinates, edges, and area
        """
        return (f"pympas(ID={self.lCellID}, "
                f"Center=({self.dLongitude_center_degree:.6f}, {self.dLatitude_center_degree:.6f}), "
                f"Edges={self.nEdge}, Area={self.dArea:.2f}m², Resolution={self.dLength:.2f}m, "
                f"Flowlines={self.nFlowline}, Neighbors={self.nNeighbor})")

    def __str__(self) -> str:
        """
        Return a concise string representation of the MPAS cell.

        Returns:
            str: Concise representation with ID, center coordinates, edges, and area
        """
        return (f"pympas(ID={self.lCellID}, "
                f"Center=({self.dLongitude_center_degree:.6f}, {self.dLatitude_center_degree:.6f}), "
                f"Edges={self.nEdge}, Area={self.dArea:.2f}m²)")

    def calculate_cell_bound(self) -> Tuple[float, float, float, float]:
        """
        Calculate the bounding box of the MPAS cell.

        This method handles the International Date Line crossing case,
        which is important for global MPAS meshes. When a cell crosses
        the date line, the bounding box coordinates are adjusted accordingly.

        Returns:
            Tuple[float, float, float, float]: Bounding box as (min_lon, min_lat, max_lon, max_lat)

        Note:
            For cells crossing the International Date Line, max_lon may be > 180°
        """
        if not self.aVertex or len(self.aVertex) == 0:
            raise ValueError("Cannot calculate bounds: no vertices defined")

        dLat_min = 90.0
        dLat_max = -90.0
        dLon_min = 180.0
        dLon_max = -180.0

        for vertex in self.aVertex:
            if hasattr(vertex, 'dLongitude_degree') and hasattr(vertex, 'dLatitude_degree'):
                dLon_max = max(dLon_max, vertex.dLongitude_degree)
                dLon_min = min(dLon_min, vertex.dLongitude_degree)
                dLat_max = max(dLat_max, vertex.dLatitude_degree)
                dLat_min = min(dLat_min, vertex.dLatitude_degree)
            else:
                raise ValueError("Vertex must have dLongitude_degree and dLatitude_degree attributes")

        # Handle International Date Line crossing
        if dLon_max - dLon_min > 180:
            # Swap values and adjust for date line crossing
            tmp = dLon_max
            dLon_max = dLon_min + 360
            dLon_min = tmp
            self.pBound = (dLon_min, dLat_min, dLon_max, dLat_max)
        else:
            self.pBound = (dLon_min, dLat_min, dLon_max, dLat_max)

        return self.pBound

    def has_this_edge(self, pEdge_in) -> int:
        """
        Check whether the MPAS cell contains a specific edge.

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
            if hasattr(pEdge, 'is_overlap') and pEdge.is_overlap(pEdge_in):
                iFlag_found = 1
                break

        return iFlag_found

    def which_edge_cross_this_vertex(self, pVertex_in) -> Tuple[int, Optional[Any]]:
        """
        Find which edge overlaps with a vertex.

        When a flowline intersects with a cell, this function identifies
        which edge is intersected by the given vertex.

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
            if hasattr(pEdge, 'check_vertex_on_edge'):
                iFlag, dummy, diff = pEdge.check_vertex_on_edge(pVertex_in)
                if iFlag == 1:
                    iFlag_found = 1
                    pEdge_out = pEdge
                    break

        return iFlag_found, pEdge_out

    def calculate_cell_area(self) -> float:
        """
        Calculate the area of the MPAS cell.

        Uses the polygon area calculation based on the vertex coordinates.
        This method works for any polygon shape with 3-10 vertices.

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
            if hasattr(vertex, 'dLongitude_degree') and hasattr(vertex, 'dLatitude_degree'):
                lons.append(vertex.dLongitude_degree)
                lats.append(vertex.dLatitude_degree)
            else:
                raise ValueError("Vertex must have dLongitude_degree and dLatitude_degree attributes")

        if not (3 <= len(lons) <= 10):
            raise ValueError(f"MPAS cell must have 3-10 vertices, got {len(lons)}")

        self.dArea = calculate_polygon_area(lons, lats)
        return self.dArea

    def calculate_edge_length(self) -> float:
        """
        Calculate the effective resolution/length of the MPAS cell.

        For MPAS cells, the effective length represents the cell resolution
        and is calculated as the square root of the area. This gives a
        characteristic length scale for the cell.

        Returns:
            float: The effective cell resolution in meters

        Raises:
            ValueError: If area is invalid
        """
        if self.dArea is None:
            raise ValueError("Area is not defined. Calculate area first.")
        if self.dArea < 0:
            raise ValueError(f"Area cannot be negative: {self.dArea}")
        if self.dArea == 0:
            # Calculate area first if not done
            self.calculate_cell_area()

        self.dLength = np.sqrt(self.dArea)
        return self.dLength

    def share_edge(self, other: 'pympas') -> int:
        """
        Check whether this MPAS cell shares an edge with another MPAS cell.

        Two MPAS cells share an edge if any of their boundary edges overlap.
        This is essential for determining adjacency in unstructured meshes.

        Args:
            other (pympas): The other MPAS cell to check

        Returns:
            int: 1 if cells share an edge, 0 if not

        Raises:
            TypeError: If other is not a pympas instance
            ValueError: If other is None
        """
        if other is None:
            raise ValueError("Other cell cannot be None")
        if not isinstance(other, pympas):
            raise TypeError(f"Other cell must be a pympas instance, got {type(other)}")

        iFlag_share = 0
        for pEdge in self.aEdge:
            for pEdge2 in other.aEdge:
                if hasattr(pEdge, 'is_overlap') and hasattr(pEdge2, 'is_overlap'):
                    if pEdge.is_overlap(pEdge2) == 1:
                        iFlag_share = 1
                        break
            if iFlag_share == 1:
                break

        return iFlag_share

    def get_mpas_properties(self) -> dict:
        """
        Get properties specific to the MPAS cell.

        Returns:
            dict: Dictionary containing MPAS-specific properties including:
                - edge_count: Number of edges/vertices
                - resolution: Effective cell resolution
                - area: Cell area
                - shape_type: Description of polygon shape
                - crosses_dateline: Whether cell crosses International Date Line
        """
        # Determine shape type based on edge count
        shape_names = {
            3: 'triangle',
            4: 'quadrilateral',
            5: 'pentagon',
            6: 'hexagon',
            7: 'heptagon',
            8: 'octagon',
            9: 'nonagon',
            10: 'decagon'
        }

        # Check if crosses dateline
        crosses_dateline = False
        if self.pBound:
            dLon_min, _, dLon_max, _ = self.pBound
            crosses_dateline = dLon_max > 180 or dLon_max - dLon_min > 180

        properties = {
            'edge_count': self.nEdge,
            'vertex_count': self.nVertex,
            'resolution': self.dLength,
            'area': self.dArea,
            'shape_type': shape_names.get(self.nEdge, f'{self.nEdge}-gon'),
            'crosses_dateline': crosses_dateline,
            'watershed_boundary': self.iFlag_watershed_boundary_burned == 1
        }
        return properties

    def get_polygon_type(self) -> str:
        """
        Get the polygon type description based on the number of edges.

        Returns:
            str: Description of the polygon type (triangle, quadrilateral, etc.)
        """
        type_map = {
            3: 'triangle',
            4: 'quadrilateral',
            5: 'pentagon',
            6: 'hexagon',
            7: 'heptagon',
            8: 'octagon',
            9: 'nonagon',
            10: 'decagon'
        }
        return type_map.get(self.nEdge, f'{self.nEdge}-sided polygon')

    def crosses_international_dateline(self) -> bool:
        """
        Check if the MPAS cell crosses the International Date Line.

        Returns:
            bool: True if cell crosses the date line
        """
        if self.pBound is None:
            self.calculate_cell_bound()

        dLon_min, _, dLon_max, _ = self.pBound
        return dLon_max > 180 or (dLon_max - dLon_min > 180)

    def get_corner_coordinates(self) -> List[Tuple[float, float]]:
        """
        Get the coordinates of all corners of the MPAS cell.

        Returns:
            List[Tuple[float, float]]: List of (longitude, latitude) tuples for each corner
        """
        corners = []
        for vertex in self.aVertex:
            if hasattr(vertex, 'dLongitude_degree') and hasattr(vertex, 'dLatitude_degree'):
                corners.append((vertex.dLongitude_degree, vertex.dLatitude_degree))
            else:
                raise ValueError("Vertex must have dLongitude_degree and dLatitude_degree attributes")

        return corners

    def set_watershed_boundary_flag(self, iFlag: int) -> None:
        """
        Set the watershed boundary flag for the MPAS cell.

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
        Check if the MPAS cell is on a watershed boundary.

        Returns:
            bool: True if cell is on watershed boundary
        """
        return self.iFlag_watershed_boundary_burned == 1

    def get_mesh_quality_metrics(self) -> dict:
        """
        Calculate mesh quality metrics for the MPAS cell.

        Returns:
            dict: Dictionary with quality metrics including:
                - aspect_ratio: Ratio of longest to shortest edge
                - regularity: How close to regular polygon (0-1, 1 is perfect)
                - compactness: Area to perimeter ratio
        """
        if len(self.aEdge) < 3:
            return {'aspect_ratio': 0, 'regularity': 0, 'compactness': 0}

        # Calculate edge lengths if available
        edge_lengths = []
        for edge in self.aEdge:
            if hasattr(edge, 'calculate_length'):
                edge_lengths.append(edge.calculate_length())
            elif hasattr(edge, 'dLength'):
                edge_lengths.append(edge.dLength)

        metrics = {}

        if edge_lengths:
            # Aspect ratio
            min_length = min(edge_lengths)
            max_length = max(edge_lengths)
            metrics['aspect_ratio'] = max_length / min_length if min_length > 0 else float('inf')

            # Regularity (coefficient of variation of edge lengths)
            mean_length = np.mean(edge_lengths)
            std_length = np.std(edge_lengths)
            metrics['regularity'] = 1.0 - (std_length / mean_length) if mean_length > 0 else 0.0

            # Compactness (area to perimeter ratio)
            perimeter = sum(edge_lengths)
            metrics['compactness'] = self.dArea / perimeter if perimeter > 0 else 0.0
        else:
            metrics = {'aspect_ratio': 1.0, 'regularity': 1.0, 'compactness': 1.0}

        return metrics

    def is_valid(self) -> bool:
        """
        Check if the MPAS cell has valid attributes.

        An MPAS cell is considered valid if:
        - It has valid center coordinates
        - It has 3-10 boundary vertices and edges
        - It has a positive area
        - It has a valid cell ID
        - All vertices and edges are properly defined

        Returns:
            bool: True if cell has valid MPAS attributes
        """
        # Check parent class validity
        if not super().is_valid():
            return False

        # Check MPAS-specific requirements
        has_valid_edge_count = 3 <= len(self.aEdge) <= 10
        has_valid_vertex_count = 3 <= len(self.aVertex) <= 10
        has_matching_counts = len(self.aEdge) == len(self.aVertex)
        has_positive_resolution = self.dLength > 0

        # Check that all vertices have required attributes
        vertices_valid = True
        for vertex in self.aVertex:
            if not (hasattr(vertex, 'dLongitude_degree') and
                   hasattr(vertex, 'dLatitude_degree')):
                vertices_valid = False
                break

        return (has_valid_edge_count and has_valid_vertex_count and
                has_matching_counts and has_positive_resolution and vertices_valid)

    def copy(self) -> 'pympas':
        """
        Create a deep copy of the MPAS cell.

        Returns:
            pympas: A new MPAS cell object with the same attributes
        """
        # Create new MPAS cell with same basic parameters
        new_mpas = pympas(self.dLongitude_center_degree,
                         self.dLatitude_center_degree,
                         self.aEdge.copy(),
                         self.aVertex.copy())

        # Copy all attributes from parent and this class
        for attr_name, attr_value in self.__dict__.items():
            if hasattr(new_mpas, attr_name):
                if isinstance(attr_value, list) and attr_value is not None:
                    setattr(new_mpas, attr_name, attr_value.copy())
                else:
                    setattr(new_mpas, attr_name, attr_value)

        return new_mpas

    def tojson(self) -> str:
        """
        Convert MPAS cell object to a JSON string.

        Serializes all cell attributes including neighbor information and
        geometric properties. Uses the custom MpasClassEncoder to handle
        numpy data types and complex objects.

        Returns:
            str: JSON string representation of the MPAS cell

        Example:
            >>> mpas_cell = pympas(-77.0, 38.0, edges, vertices)
            >>> json_str = mpas_cell.tojson()
        """
        aSkip = ['aEdge', 'aFlowline', 'dLongitude_radian', 'dLatitude_radian',
                 'wkt', 'pPoint_start', 'pPoint_end', 'pBound']

        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)

        sJson = json.dumps(obj,
                          sort_keys=True,
                          indent=4,
                          ensure_ascii=True,
                          cls=MpasClassEncoder)
        return sJson
