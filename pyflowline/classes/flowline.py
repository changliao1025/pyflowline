
"""
PyFlowline class for representing flowlines in river network models.

This module extends the pypolyline class from pyearth to add flowline-specific
attributes and functionality. A flowline represents a connected sequence of edges
forming a continuous path in a stream network, with associated hydrological properties.
"""

import copy
import json
from json import JSONEncoder
from typing import List, Dict, Any, Optional, Union
import numpy as np

from pyearth.toolbox.mesh.point import pypoint
from pyearth.toolbox.mesh.line import pyline
from pyearth.toolbox.mesh.polyline import pypolyline
from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.edge import pyedge

class FlowlineClassEncoder(JSONEncoder):
    """
    Custom JSON encoder for pyflowline objects.

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

        return JSONEncoder.default(self, obj)


class pyflowline(pypolyline):
    """
    Flowline class for river network representation.

    Extends the pypolyline class to represent a connected sequence of edges
    forming a continuous path in a stream network. A flowline has hydrological
    attributes such as drainage area, stream order, and connectivity information
    with upstream and downstream flowlines.

    Attributes:
        lFlowlineID (int): Unique identifier for the flowline (default: -1)
        lFlowlineID_downstream (int): ID of the downstream flowline (default: -1)
        lFlowlineIndex (int): Array index for the flowline (default: -1)
        lIndex_upstream (int): Index of upstream connection (default: -1)
        lIndex_downstream (int): Index of downstream connection (default: -1)
        iFlag_keep (int): Flag for simplification algorithm (1=keep, 0=remove, default: 1)
        iFlag_dam (int): Flag indicating if flowline has a dam (0=no, 1=yes, default: 0)
        iFlag_endorheic (int): Flag for endorheic basin (0=no, 1=yes, default: 0)
        iFlag_right (int): Right side flag (default: 0)
        iFlag_left (int): Left side flag (default: 0)
        lNHDPlusID (int): NHDPlus dataset ID if applicable (default: -1)
        dSinuosity (float): Sinuosity ratio (flowline length / straight distance, default: 0.0)
        dDrainage_area (float): Drainage area in square meters (default: 0.0)
        iStream_segment (int): Stream segment identifier (default: -1)
        iStream_order (int): Strahler stream order (default: -1)
        nEdge (int): Number of edges in the flowline
        nVertex (int): Number of vertices in the flowline
        pVertex_start (pyvertex): Starting vertex of the flowline
        pVertex_end (pyvertex): Ending vertex of the flowline
        aEdge (List[pyedge]): List of edges forming the flowline
        aVertex (List[pyvertex]): List of vertices along the flowline
        lFlowlineIndex_downstream (int): Index of downstream flowline (default: None)
        aFlowline_upstream (List): List of upstream flowlines (default: None)
        aFlowlineID_start_start (List): IDs of flowlines starting at start (default: None)
        aFlowlineID_start_end (List): IDs of flowlines ending at start (default: None)
        aFlowlineID_end_start (List): IDs of flowlines starting at end (default: None)
        aFlowlineID_end_end (List): IDs of flowlines ending at end (default: None)

    Args:
        aEdge (List[pyedge]): A list of edge objects forming the flowline

    Raises:
        TypeError: If aEdge is not a list
        ValueError: If aEdge is empty or contains invalid edges

    Example:
        >>> v1 = pyvertex({'dLongitude_degree': -77.0, 'dLatitude_degree': 38.0})
        >>> v2 = pyvertex({'dLongitude_degree': -77.1, 'dLatitude_degree': 38.1})
        >>> v3 = pyvertex({'dLongitude_degree': -77.2, 'dLatitude_degree': 38.2})
        >>> e1 = pyedge(v1, v2)
        >>> e2 = pyedge(v2, v3)
        >>> flowline = pyflowline([e1, e2])
        >>> flowline.set_flowline_id(100)
        >>> print(flowline)
        pyflowline(ID=100, Edges=2, Length=31176.90m)
    """

    def __init__(self, aEdge: List[Union[pyedge, pyline]]) -> None:
        """
        Initialize a flowline object.

        Converts edge inputs to pyline objects for the parent class and stores
        them as pyedge objects. Initializes all flowline-specific attributes.

        Args:
            aEdge (List[pyedge]): A list of edge objects forming the flowline

        Raises:
            TypeError: If aEdge is not a list
            ValueError: If aEdge is empty or contains invalid edges
        """
        # Validate input
        if not isinstance(aEdge, list):
            raise TypeError(f"aEdge must be a list, got {type(aEdge)}")

        if len(aEdge) == 0:
            raise ValueError("aEdge cannot be empty - flowline must have at least one edge")

        # Convert all elements to pyline for parent class
        aLine = []
        for i, edge in enumerate(aEdge):
            if edge is None:
                raise ValueError(f"Edge at index {i} cannot be None")

            if hasattr(edge.pVertex_start, '__dict__') and hasattr(edge.pVertex_end, '__dict__'):
                point_start = pypoint(edge.pVertex_start.__dict__)
                point_end = pypoint(edge.pVertex_end.__dict__)
                line = pyline(point_start, point_end)
                aLine.append(line)
            else:
                raise TypeError(f"Edge at index {i} must be a pyedge object with valid vertices")

        # Initialize parent class
        super().__init__(aLine)

        # Flowline identification attributes
        self.lFlowlineID: int = -1
        self.lFlowlineID_downstream: int = -1  # if braided, then we need a list
        self.lFlowlineIndex: int = -1
        self.lIndex_upstream: int = -1
        self.lIndex_downstream: int = -1

        # Flowline status flags
        self.iFlag_keep: int = 1  # used for simplification algorithm
        self.iFlag_dam: int = 0
        self.iFlag_endorheic: int = 0  # used for endorheic basin
        self.iFlag_right: int = 0
        self.iFlag_left: int = 0

        # External dataset attributes
        self.lNHDPlusID: int = -1

        # Geometric and hydrological attributes
        self.dSinuosity: float = 0.0
        self.dDrainage_area: float = 0.0
        self.iStream_segment: int = -1
        self.iStream_order: int = -1

        # Store vertices and edges
        self.pVertex_start: pyvertex = pyvertex(self.pPoint_start.__dict__)
        self.pVertex_end: pyvertex = pyvertex(self.pPoint_end.__dict__)
        # Convert lines to edges properly
        self.aEdge = aEdge
        self.aVertex: List[pyvertex] = [pyvertex(v.__dict__) if not isinstance(v, pyvertex) else v for v in self.aPoint]

        # Counts
        self.nEdge: int = len(self.aEdge)
        self.nVertex: int = len(self.aVertex)

        # Connection lists (initialized as None, populated during network topology analysis)
        self.lFlowlineIndex_downstream: Optional[int] = None
        self.aFlowline_upstream: Optional[List] = None
        self.aFlowlineID_start_start: Optional[List[int]] = None
        self.aFlowlineID_start_end: Optional[List[int]] = None
        self.aFlowlineID_end_start: Optional[List[int]] = None
        self.aFlowlineID_end_end: Optional[List[int]] = None

    def __repr__(self) -> str:
        """
        Return a detailed string representation of the flowline.

        Returns:
            str: Detailed representation including ID, indices, edge count, and length
        """
        return (f"pyflowline(ID={self.lFlowlineID}, Index={self.lFlowlineIndex}, "
                f"Edges={self.nEdge}, Vertices={self.nVertex}, "
                f"Length={self.dLength:.2f}m, StreamOrder={self.iStream_order}, "
                f"Sinuosity={self.dSinuosity:.3f})")

    def __str__(self) -> str:
        """
        Return a concise string representation of the flowline.

        Returns:
            str: Concise representation with ID, edge count, and length
        """
        return f"pyflowline(ID={self.lFlowlineID}, Edges={self.nEdge}, Length={self.dLength:.2f}m)"

    def __hash__(self) -> int:
        """
        Return hash value for the flowline.

        Uses flowline ID for hashing to allow flowlines to be used in sets
        and as dictionary keys.

        Returns:
            int: Hash value based on flowline ID
        """
        return hash((self.pVertex_start, self.pVertex_end))

    def __eq__(self, other: Any) -> bool:
        """
        Check if two flowlines are equal based on their IDs.

        Args:
            other: Another object to compare with

        Returns:
            bool: True if flowlines have the same ID
        """
        if not isinstance(other, pyflowline):
            return NotImplemented
        return self.lFlowlineID == other.lFlowlineID

    def set_flowline_id(self, lFlowlineID: int) -> None:
        """
        Set the flowline ID.

        Args:
            lFlowlineID (int): New flowline ID

        Raises:
            TypeError: If lFlowlineID is not an integer
        """
        if not isinstance(lFlowlineID, (int, np.integer)):
            raise TypeError(f"Flowline ID must be an integer, got {type(lFlowlineID)}")
        self.lFlowlineID = int(lFlowlineID)

    def set_flowline_index(self, lFlowlineIndex: int) -> None:
        """
        Set the flowline array index.

        Args:
            lFlowlineIndex (int): New flowline index

        Raises:
            TypeError: If lFlowlineIndex is not an integer
        """
        if not isinstance(lFlowlineIndex, (int, np.integer)):
            raise TypeError(f"Flowline index must be an integer, got {type(lFlowlineIndex)}")
        self.lFlowlineIndex = int(lFlowlineIndex)

    def set_stream_order(self, iStream_order: int) -> None:
        """
        Set the Strahler stream order.

        Args:
            iStream_order (int): Stream order value

        Raises:
            TypeError: If iStream_order is not an integer
            ValueError: If iStream_order is negative
        """
        if not isinstance(iStream_order, (int, np.integer)):
            raise TypeError(f"Stream order must be an integer, got {type(iStream_order)}")
        if iStream_order < 0:
            raise ValueError(f"Stream order cannot be negative, got {iStream_order}")
        self.iStream_order = int(iStream_order)

    def set_drainage_area(self, dDrainage_area: float) -> None:
        """
        Set the drainage area.

        Args:
            dDrainage_area (float): Drainage area in square meters

        Raises:
            TypeError: If dDrainage_area is not numeric
            ValueError: If dDrainage_area is negative
        """
        if not isinstance(dDrainage_area, (int, float, np.number)):
            raise TypeError(f"Drainage area must be numeric, got {type(dDrainage_area)}")
        if dDrainage_area < 0:
            raise ValueError(f"Drainage area cannot be negative, got {dDrainage_area}")
        self.dDrainage_area = float(dDrainage_area)

    def is_valid(self) -> bool:
        """
        Check if the flowline has valid attributes.

        A flowline is considered valid if:
        - It has at least one edge
        - All edges are valid
        - All vertices are valid
        - It has positive length
        - It has a valid ID or index

        Returns:
            bool: True if flowline has valid attributes
        """
        if self.nEdge == 0:
            return False

        if self.dLength <= 0:
            return False

        # Check if all edges are valid
        if not all(edge.is_valid() for edge in self.aEdge):
            return False

        # Check if all vertices are valid
        if not all(vertex.is_valid() for vertex in self.aVertex):
            return False

        # Check if flowline has a valid ID or index
        has_valid_id = self.lFlowlineID > 0 or self.lFlowlineIndex >= 0

        return has_valid_id

    def is_headwater(self) -> bool:
        """
        Check if this flowline is a headwater (no upstream flowlines).

        Returns:
            bool: True if flowline has no upstream connections
        """
        return (self.aFlowline_upstream is None or
                len(self.aFlowline_upstream) == 0)

    def is_outlet(self) -> bool:
        """
        Check if this flowline is an outlet (no downstream flowline).

        Returns:
            bool: True if flowline has no downstream connection
        """
        return (self.lFlowlineIndex_downstream is None or
                self.lFlowlineIndex_downstream < 0)

    def check_upstream(self, other: 'pyflowline') -> bool:
        """
        Check whether another flowline is directly upstream.

        A flowline is upstream if its end vertex matches this flowline's start vertex.

        Args:
            other (pyflowline): The other flowline to check

        Returns:
            bool: True if other flowline is directly upstream, False otherwise

        Raises:
            TypeError: If other is not a pyflowline
        """
        if not isinstance(other, pyflowline):
            raise TypeError(f"Expected pyflowline, got {type(other)}")

        v0 = self.pVertex_start
        v3 = other.pVertex_end

        return v0 == v3

    def check_downstream(self, other: 'pyflowline') -> bool:
        """
        Check whether another flowline is directly downstream.

        A flowline is downstream if its start vertex matches this flowline's end vertex.

        Args:
            other (pyflowline): The other flowline to check

        Returns:
            bool: True if other flowline is directly downstream, False otherwise

        Raises:
            TypeError: If other is not a pyflowline
        """
        if not isinstance(other, pyflowline):
            raise TypeError(f"Expected pyflowline, got {type(other)}")

        v1 = self.pVertex_end
        v2 = other.pVertex_start

        return v1 == v2

    def merge_upstream(self, other: 'pyflowline') -> 'pyflowline':
        """
        Merge with an upstream flowline to create a new combined flowline.

        The upstream flowline's end vertex must match this flowline's start vertex.
        The resulting flowline inherits attributes from both flowlines, with the
        dam flag taking the maximum value.

        Args:
            other (pyflowline): The upstream flowline to merge with

        Returns:
            pyflowline: A new merged flowline

        Raises:
            TypeError: If other is not a pyflowline
            ValueError: If flowlines are not properly connected
        """
        if not isinstance(other, pyflowline):
            raise TypeError(f"Expected pyflowline, got {type(other)}")

        pFlowline_out = copy.deepcopy(other)

        iFlag_dam1 = other.iFlag_dam
        pVertex_start1 = other.pVertex_start
        pVertex_end1 = other.pVertex_end
        nVertex1 = other.nVertex
        nEdge1 = other.nEdge

        pVertex_start2 = self.pVertex_start
        pVertex_end2 = self.pVertex_end
        nVertex2 = self.nVertex
        nEdge2 = self.nEdge
        iFlag_dam2 = self.iFlag_dam

        if pVertex_end1 == pVertex_start2:
            # This is the expected connection - merge the flowlines
            nVertex = nVertex1 + nVertex2 - 1
            nEdge = nVertex - 1

            # Combine edges
            aEdge = copy.deepcopy(other.aEdge)
            for i in range(nEdge2):
                aEdge.append(self.aEdge[i])

            # Combine vertices (skip duplicate connection point)
            aVertex = copy.deepcopy(other.aVertex)
            for i in range(1, nVertex2):
                aVertex.append(self.aVertex[i])

            # Update merged flowline attributes
            pFlowline_out.iFlag_dam = max(iFlag_dam1, iFlag_dam2)
            pFlowline_out.aEdge = aEdge
            pFlowline_out.aVertex = aVertex
            pFlowline_out.nEdge = nEdge
            pFlowline_out.nVertex = nVertex
            pFlowline_out.dLength = self.dLength + other.dLength
            pFlowline_out.pVertex_start = pVertex_start1
            pFlowline_out.pVertex_end = pVertex_end2
        else:
            raise ValueError("Flowlines are not properly connected for merging. "
                           f"Upstream end ({pVertex_end1}) must match downstream start ({pVertex_start2})")

        return pFlowline_out

    def copy_attributes(self, other: 'pyflowline') -> None:
        """
        Copy attributes from another flowline.

        Copies identification, status, and hydrological attributes but not
        the geometric data (edges and vertices).

        Args:
            other (pyflowline): The source flowline to copy attributes from

        Raises:
            TypeError: If other is not a pyflowline
        """
        if not isinstance(other, pyflowline):
            raise TypeError(f"Expected pyflowline, got {type(other)}")

        self.lFlowlineID = other.lFlowlineID
        self.lFlowlineID_downstream = other.lFlowlineID_downstream
        self.lFlowlineIndex = other.lFlowlineIndex
        self.iFlag_dam = other.iFlag_dam
        self.lNHDPlusID = other.lNHDPlusID
        self.iStream_segment = other.iStream_segment
        self.iStream_order = other.iStream_order
        self.dDrainage_area = other.dDrainage_area
        self.iFlag_right = other.iFlag_right
        self.iFlag_left = other.iFlag_left
        self.iFlag_keep = other.iFlag_keep
        self.lFlowlineIndex_downstream = other.lFlowlineIndex_downstream

    def calculate_flowline_sinuosity(self) -> float:
        """
        Calculate the sinuosity of the flowline.

        Sinuosity is the ratio of the flowline length to the straight-line
        distance between start and end points. A value of 1.0 indicates a
        perfectly straight flowline, while higher values indicate more sinuous paths.

        Returns:
            float: Sinuosity ratio (>= 1.0)

        Raises:
            ValueError: If straight-line distance is zero (identical endpoints)
        """
        pVertex_start = self.pVertex_start
        pVertex_end = self.pVertex_end
        dDistance = pVertex_start.calculate_distance(pVertex_end)

        if dDistance == 0:
            raise ValueError("Cannot calculate sinuosity: start and end vertices are identical")

        self.dSinuosity = self.dLength / dDistance
        return self.dSinuosity

    def get_sinuosity(self) -> float:
        """
        Get the sinuosity value, calculating it if not already computed.

        Returns:
            float: Sinuosity ratio
        """
        if self.dSinuosity == 0.0:
            self.calculate_flowline_sinuosity()
        return self.dSinuosity

    def copy(self) -> 'pyflowline':
        """
        Create a deep copy of the flowline.

        Returns:
            pyflowline: A new flowline object with the same attributes
        """
        # Copy edges
        aEdge_copy = [edge.copy() for edge in self.aEdge]

        # Create new flowline
        new_flowline = pyflowline(aEdge_copy)

        # Copy all attributes
        new_flowline.copy_attributes(self)
        new_flowline.dSinuosity = self.dSinuosity
        new_flowline.iFlag_endorheic = self.iFlag_endorheic
        new_flowline.lIndex_upstream = self.lIndex_upstream
        new_flowline.lIndex_downstream = self.lIndex_downstream

        return new_flowline

    def reverse(self) -> 'pyflowline':
        """
        Create a reversed copy of the flowline.

        Returns a new flowline with all edges reversed and vertices in

        reverse order.

        Returns:
            pyflowline: A new flowline with reversed direction

        Raises:
            ValueError: If flowline cannot be reversed
        """
        # Reverse edges
        aEdge_reversed = [edge.reverse() for edge in reversed(self.aEdge)]

        # Create new flowline
        reversed_flowline = pyflowline(aEdge_reversed)

        # Copy attributes
        reversed_flowline.copy_attributes(self)
        reversed_flowline.pVertex_start = self.pVertex_end.copy()
        reversed_flowline.pVertex_end = self.pVertex_start.copy()
        reversed_flowline.dSinuosity = self.dSinuosity

        # Swap upstream/downstream indices
        reversed_flowline.lIndex_upstream = self.lIndex_downstream
        reversed_flowline.lIndex_downstream = self.lIndex_upstream

        return reversed_flowline

    def get_length(self) -> float:
        """
        Get the total length of the flowline in meters.

        Returns:
            float: Total length in meters
        """
        return self.dLength

    def get_edge_count(self) -> int:
        """
        Get the number of edges in the flowline.

        Returns:
            int: Number of edges
        """
        return self.nEdge

    def get_vertex_count(self) -> int:
        """
        Get the number of vertices in the flowline.

        Returns:
            int: Number of vertices
        """
        return self.nVertex

    def get_upstream_count(self) -> int:
        """
        Get the number of upstream flowlines.

        Returns:
            int: Number of upstream flowlines
        """
        if self.aFlowline_upstream is None:
            return 0
        return len(self.aFlowline_upstream)

    def has_dam(self) -> bool:
        """
        Check if the flowline has a dam.

        Returns:
            bool: True if flowline has a dam
        """
        return self.iFlag_dam == 1

    def is_endorheic(self) -> bool:
        """
        Check if the flowline is in an endorheic basin.

        Returns:
            bool: True if flowline is in an endorheic basin
        """
        return self.iFlag_endorheic == 1

    def should_keep(self) -> bool:
        """
        Check if the flowline should be kept during simplification.

        Returns:
            bool: True if flowline should be kept
        """
        return self.iFlag_keep == 1

    def mark_for_removal(self) -> None:
        """
        Mark the flowline for removal during simplification.
        """
        self.iFlag_keep = 0

    def mark_for_keeping(self) -> None:
        """
        Mark the flowline to be kept during simplification.
        """
        self.iFlag_keep = 1

    def tojson(self) -> str:
        """
        Convert pyflowline object to a JSON string.

        Serializes most flowline attributes while excluding internal arrays
        and complex objects that are better serialized separately.

        Returns:
            str: JSON string representation of the flowline

        Example:
            >>> flowline = pyflowline([edge1, edge2])
            >>> json_str = flowline.tojson()
        """
        aSkip = ['aEdge', 'aVertex', 'aFlowlineID_start_start','aLine', 'aPoint', 'pPoint_start', 'pPoint_end',
                 'aFlowlineID_start_end', 'aFlowlineID_end_start',
                 'aFlowlineID_end_end', 'aFlowline_upstream']

        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)

        sJson = json.dumps(obj,
                          sort_keys=True,
                          indent=4,
                          ensure_ascii=True,
                          cls=FlowlineClassEncoder)
        return sJson
