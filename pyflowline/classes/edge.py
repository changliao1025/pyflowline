"""
PyEdge class for representing edges in flowline networks.

This module extends the pyline class from pyearth to add flowline-specific
edge attributes and functionality. An edge represents a directed connection
between two vertices in the flowline network.
"""

import json
from json import JSONEncoder
from typing import Dict, Any, Optional
import numpy as np
from pyearth.toolbox.mesh.point import pypoint
from pyearth.toolbox.mesh.line import pyline
from pyflowline.classes.vertex import pyvertex
class EdgeClassEncoder(JSONEncoder):
    """
    Custom JSON encoder for pyedge objects.

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


class pyedge(pyline):
    """
    Edge class for flowline network representation.

    Extends the pyline class to include edge-specific attributes used
    in flowline topology and stream network analysis. An edge represents
    a directed connection between two vertices, forming the basic building
    block of flowlines in the network.

    Attributes:
        lEdgeID (int): Unique identifier for the edge (default: -1)
        lEdgeIndex (int): Array index for the edge in computational operations (default: -1)
        lIndex_upstream (int): Index of the upstream edge in the network (default: -1)
        lIndex_downstream (int): Index of the downstream edge in the network (default: -1)
        pVertex_start (pyvertex): Starting vertex of the edge
        pVertex_end (pyvertex): Ending vertex of the edge

    Inherits from pyline:
        pPoint_start (pypoint): Start point of the line
        pPoint_end (pypoint): End point of the line
        dLength (float): Length of the edge in meters
        dSlope (float): Slope of the edge

    Args:
        pVertex_start_in (pyvertex or pypoint): The starting vertex/point
        pVertex_end_in (pyvertex or pypoint): The ending vertex/point

    Raises:
        TypeError: If vertex parameters are not pyvertex or pypoint objects
        ValueError: If vertices are identical (zero-length edge)

    Example:
        >>> v1 = pyvertex({'dLongitude_degree': -77.0, 'dLatitude_degree': 38.0})
        >>> v2 = pyvertex({'dLongitude_degree': -77.1, 'dLatitude_degree': 38.1})
        >>> edge = pyedge(v1, v2)
        >>> edge.set_edge_id(100)
        >>> print(edge)
        pyedge(ID=100, Length=15588.45m)
    """

    def __init__(self, pVertex_start_in, pVertex_end_in) -> None:
        """
        Initialize a pyedge object.

        Converts vertex inputs to pypoint objects if needed and initializes
        the parent pyline class. Sets default values for edge-specific attributes.

        Args:
            pVertex_start_in (pyvertex or pypoint): The starting vertex
            pVertex_end_in (pyvertex or pypoint): The ending vertex

        Raises:
            TypeError: If inputs are not pyvertex or pypoint objects
            ValueError: If start and end vertices are identical
        """
        # Validate inputs
        if pVertex_start_in is None or pVertex_end_in is None:
            raise ValueError("Start and end vertices cannot be None")

        # Convert pVertex_start_in and pVertex_end_in to pypoint if they are not
        if not isinstance(pVertex_start_in, pypoint):
            if hasattr(pVertex_start_in, '__dict__'):
                pPoint_start_in = pypoint(pVertex_start_in.__dict__)
            else:
                raise TypeError(f"Start vertex must be pypoint or pyvertex, got {type(pVertex_start_in)}")
        else:
            pPoint_start_in = pVertex_start_in

        if not isinstance(pVertex_end_in, pypoint):
            if hasattr(pVertex_end_in, '__dict__'):
                pPoint_end_in = pypoint(pVertex_end_in.__dict__)
            else:
                raise TypeError(f"End vertex must be pypoint or pyvertex, got {type(pVertex_end_in)}")
        else:
            pPoint_end_in = pVertex_end_in

        # Check if vertices are identical
        if pPoint_start_in == pPoint_end_in:
            raise ValueError("Start and end vertices cannot be identical (zero-length edge)")

        # Initialize parent class
        super().__init__(pPoint_start_in, pPoint_end_in)

        # Edge-specific attributes with defaults
        # lEdgeID: Unique identifier for the edge in the network
        self.lEdgeID: int = -1

        # lEdgeIndex: Used for array indexing in computational operations
        self.lEdgeIndex: int = -1

        # lIndex_upstream: Index of the upstream edge connected to this edge
        self.lIndex_upstream: int = -1

        # lIndex_downstream: Index of the downstream edge connected to this edge
        self.lIndex_downstream: int = -1

        # Store the original vertex objects for later use
        self.pVertex_start = pVertex_start_in if isinstance(pVertex_start_in, pyvertex) else pyvertex(pVertex_start_in.__dict__)
        self.pVertex_end = pVertex_end_in if isinstance(pVertex_end_in, pyvertex) else pyvertex(pVertex_end_in.__dict__)

    def __repr__(self) -> str:
        """
        Return a detailed string representation of the edge.

        Returns:
            str: Detailed representation including edge ID, index, and length
        """
        return (f"pyedge(ID={self.lEdgeID}, Index={self.lEdgeIndex}, "
                f"Length={self.dLength:.2f}m, "
                f"Upstream={self.lIndex_upstream}, Downstream={self.lIndex_downstream})")

    def __str__(self) -> str:
        """
        Return a concise string representation of the edge.

        Returns:
            str: Concise representation with ID and length
        """
        return f"pyedge(ID={self.lEdgeID}, Length={self.dLength:.2f}m)"

    def __hash__(self) -> int:
        """
        Return hash value for the edge.

        Uses parent class hash which is based on the line geometry.
        This allows edges to be used in sets and as dictionary keys.

        Returns:
            int: Hash value based on geometry
        """
        return super().__hash__()

    def __eq__(self, other: Any) -> bool:
        """
        Check if two edges are equal based on their geometry.

        Edges are considered equal if their start and end points match
        within a precision threshold, regardless of their IDs or indices.

        Args:
            other: Another object to compare with

        Returns:
            bool: True if edges have the same geometry
        """
        if not isinstance(other, pyedge):
            return NotImplemented
        return super().__eq__(other)

    def tojson(self) -> str:
        """
        Convert edge object to a JSON string.

        Serializes all edge attributes including vertex information.
        Uses the custom EdgeClassEncoder to handle numpy data types
        and pyvertex objects.

        Returns:
            str: JSON string representation of the edge

        Example:
            >>> edge = pyedge(v1, v2)
            >>> json_str = edge.tojson()
        """
        aSkip = ['dLongitude_radian', 'dLatitude_radian', 'wkt',
                 'pPoint_start', 'pPoint_end']

        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)

        sJson = json.dumps(obj,
                          sort_keys=True,
                          indent=4,
                          ensure_ascii=True,
                          cls=EdgeClassEncoder)
        return sJson

    def set_edge_id(self, lEdgeID: int) -> None:
        """
        Set the edge ID.

        Args:
            lEdgeID (int): New edge ID

        Raises:
            TypeError: If lEdgeID is not an integer
        """
        if not isinstance(lEdgeID, (int, np.integer)):
            raise TypeError(f"Edge ID must be an integer, got {type(lEdgeID)}")
        self.lEdgeID = int(lEdgeID)

    def set_edge_index(self, lEdgeIndex: int) -> None:
        """
        Set the edge array index.

        Args:
            lEdgeIndex (int): New edge index

        Raises:
            TypeError: If lEdgeIndex is not an integer
        """
        if not isinstance(lEdgeIndex, (int, np.integer)):
            raise TypeError(f"Edge index must be an integer, got {type(lEdgeIndex)}")
        self.lEdgeIndex = int(lEdgeIndex)

    def set_upstream_index(self, lIndex_upstream: int) -> None:
        """
        Set the upstream edge index.

        Args:
            lIndex_upstream (int): New upstream edge index

        Raises:
            TypeError: If lIndex_upstream is not an integer
        """
        if not isinstance(lIndex_upstream, (int, np.integer)):
            raise TypeError(f"Upstream index must be an integer, got {type(lIndex_upstream)}")
        self.lIndex_upstream = int(lIndex_upstream)

    def set_downstream_index(self, lIndex_downstream: int) -> None:
        """
        Set the downstream edge index.

        Args:
            lIndex_downstream (int): New downstream edge index

        Raises:
            TypeError: If lIndex_downstream is not an integer
        """
        if not isinstance(lIndex_downstream, (int, np.integer)):
            raise TypeError(f"Downstream index must be an integer, got {type(lIndex_downstream)}")
        self.lIndex_downstream = int(lIndex_downstream)

    def is_valid(self) -> bool:
        """
        Check if the edge has valid attributes.

        An edge is considered valid if:
        - It has non-zero length
        - Both vertices are valid
        - It has a valid ID or index

        Returns:
            bool: True if edge has valid attributes
        """
        has_valid_length = self.dLength > 0
        has_valid_vertices = (hasattr(self, 'pVertex_start') and
                            hasattr(self, 'pVertex_end') and
                            self.pVertex_start.is_valid() and
                            self.pVertex_end.is_valid())
        has_valid_id = self.lEdgeID > 0 or self.lEdgeIndex >= 0
        return has_valid_length and has_valid_vertices and has_valid_id

    def copy(self) -> 'pyedge':
        """
        Create a deep copy of the edge.

        Returns:
            pyedge: A new edge object with the same attributes
        """
        new_edge = pyedge(self.pVertex_start.copy(), self.pVertex_end.copy())
        new_edge.lEdgeID = self.lEdgeID
        new_edge.lEdgeIndex = self.lEdgeIndex
        new_edge.lIndex_upstream = self.lIndex_upstream
        new_edge.lIndex_downstream = self.lIndex_downstream
        return new_edge

    def reverse(self) -> 'pyedge':
        """
        Create a reversed copy of the edge.

        Returns a new edge with start and end vertices swapped,
        and upstream/downstream indices swapped.

        Returns:
            pyedge: A new edge object with reversed direction
        """
        reversed_edge = pyedge(self.pVertex_end.copy(), self.pVertex_start.copy())
        reversed_edge.lEdgeID = self.lEdgeID
        reversed_edge.lEdgeIndex = self.lEdgeIndex
        reversed_edge.lIndex_upstream = self.lIndex_downstream
        reversed_edge.lIndex_downstream = self.lIndex_upstream
        return reversed_edge

    def get_midpoint(self) -> pyvertex:
        """
        Calculate and return the midpoint of the edge.

        Returns:
            pyvertex: A new vertex at the midpoint of the edge
        """
        mid_lon = (self.pVertex_start.dLongitude_degree + self.pVertex_end.dLongitude_degree) / 2.0
        mid_lat = (self.pVertex_start.dLatitude_degree + self.pVertex_end.dLatitude_degree) / 2.0

        param = {
            'dLongitude_degree': mid_lon,
            'dLatitude_degree': mid_lat
        }
        return pyvertex(param)

    def is_connected_to(self, other: 'pyedge') -> bool:
        """
        Check if this edge is connected to another edge.

        Two edges are connected if they share a common vertex.

        Args:
            other (pyedge): Another edge to check connection with

        Returns:
            bool: True if edges share a common vertex
        """
        if not isinstance(other, pyedge):
            return False

        return (self.pVertex_end == other.pVertex_start or
                self.pVertex_end == other.pVertex_end or
                self.pVertex_start == other.pVertex_start or
                self.pVertex_start == other.pVertex_end)

    def is_upstream_of(self, other: 'pyedge') -> bool:
        """
        Check if this edge is directly upstream of another edge.

        Args:
            other (pyedge): Another edge to check

        Returns:
            bool: True if this edge's end vertex matches other's start vertex
        """
        if not isinstance(other, pyedge):
            return False
        return self.pVertex_end == other.pVertex_start

    def is_downstream_of(self, other: 'pyedge') -> bool:
        """
        Check if this edge is directly downstream of another edge.

        Args:
            other (pyedge): Another edge to check

        Returns:
            bool: True if this edge's start vertex matches other's end vertex
        """
        if not isinstance(other, pyedge):
            return False
        return self.pVertex_start == other.pVertex_end
