"""
PyVertex class for representing vertices in flowline networks.

This module extends the pypoint class from pyearth to add flowline-specific
vertex attributes and functionality.
"""

import json
from json import JSONEncoder
from typing import Dict, Any, Optional
import numpy as np
from pyearth.toolbox.mesh.point import pypoint


class VertexClassEncoder(JSONEncoder):
    """
    Custom JSON encoder for pyvertex objects.

    Handles numpy data types and converts them to native Python types
    for JSON serialization.
    """

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, pyvertex):
            return json.loads(obj.tojson())
        return JSONEncoder.default(self, obj)


class pyvertex(pypoint):
    """
    Vertex class for flowline network representation.

    Extends the pypoint class to include vertex-specific attributes used
    in flowline topology and stream network analysis. A vertex represents
    a point in the flowline network that may be shared between multiple
    flowlines or serve as a confluence point.

    Attributes:
        lVertexIndex (int): Array index for the vertex (default: -1)
        lVertexID (int): Unique identifier for the vertex (default: 1)
        lFlowlineID (int): Associated flowline ID, primarily used for
                          intersection operations (default: -1)

    Inherits from pypoint:
        dLongitude_degree (float): Longitude in degrees
        dLatitude_degree (float): Latitude in degrees
        dX_meter (float): X coordinate in meters
        dY_meter (float): Y coordinate in meters
        dZ_meter (float): Z coordinate in meters
        dElevation (float): Elevation value

    Args:
        aParameter (dict): Dictionary containing vertex parameters.
                          Required keys: 'dLongitude_degree', 'dLatitude_degree'
                          Optional keys: 'lVertexIndex', 'lVertexID', 'lFlowlineID',
                                       'x', 'y', 'z', 'dElevation'

    Example:
        >>> vertex_params = {
        ...     'dLongitude_degree': -77.0,
        ...     'dLatitude_degree': 38.0,
        ...     'lVertexID': 5
        ... }
        >>> vertex = pyvertex(vertex_params)
        >>> print(vertex)
        pyvertex(ID=5, Lon=-77.0, Lat=38.0)
    """

    def __init__(self, aParameter: Dict[str, Any]) -> None:
        """
        Initialize a vertex object.

        Args:
            aParameter (dict): Dictionary containing vertex parameters.
                             Must include 'dLongitude_degree' and 'dLatitude_degree'.

        Raises:
            ValueError: If required parameters are missing
            TypeError: If aParameter is not a dictionary
        """
        if not isinstance(aParameter, dict):
            raise TypeError(f"aParameter must be a dictionary, got {type(aParameter)}")

        # Initialize parent class (pypoint)
        super().__init__(aParameter)

        # Vertex-specific attributes with defaults
        # lVertexIndex: Used for array indexing in computational operations
        self.lVertexIndex = int(aParameter.get("lVertexIndex", -1))

        # lVertexID: Unique identifier for the vertex in the network
        self.lVertexID = int(aParameter.get("lVertexID", 1))

        # lFlowlineID: Associated flowline ID (used during intersection operations)
        self.lFlowlineID = int(aParameter.get("lFlowlineID", -1))

    def __repr__(self) -> str:
        """
        Return a detailed string representation of the vertex.

        Returns:
            str: Detailed representation including vertex ID and coordinates
        """
        return (
            f"pyvertex(ID={self.lVertexID}, Index={self.lVertexIndex}, "
            f"Lon={self.dLongitude_degree:.6f}, Lat={self.dLatitude_degree:.6f}, "
            f"FlowlineID={self.lFlowlineID})"
        )

    def __str__(self) -> str:
        """
        Return a concise string representation of the vertex.

        Returns:
            str: Concise representation with ID and coordinates
        """
        return f"pyvertex(ID={self.lVertexID}, Lon={self.dLongitude_degree:.6f}, Lat={self.dLatitude_degree:.6f})"

    def __hash__(self) -> int:
        """
        Return hash value for the vertex.

        Uses parent class hash which is based on rounded coordinates.
        This allows vertices to be used in sets and as dictionary keys.

        Returns:
            int: Hash value based on coordinates
        """
        return super().__hash__()

    def __eq__(self, other: Any) -> bool:
        """
        Check if two vertices are equal based on their coordinates.

        Vertices are considered equal if their coordinates match within
        a precision threshold, regardless of their IDs or indices.

        Args:
            other: Another object to compare with

        Returns:
            bool: True if vertices have the same coordinates
        """
        if not isinstance(other, pyvertex):
            return NotImplemented
        return super().__eq__(other)

    def tojson(self) -> str:
        """
        Convert vertex object to a JSON string.

        Serializes all vertex attributes except internal computed values
        like radians. Uses the custom VertexClassEncoder to handle
        numpy data types.

        Returns:
            str: JSON string representation of the vertex

        Example:
            >>> vertex = pyvertex({'dLongitude_degree': -77.0, 'dLatitude_degree': 38.0})
            >>> json_str = vertex.tojson()
        """
        aSkip = ["dLongitude_radian", "dLatitude_radian", "wkt"]

        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)

        sJson = json.dumps(
            obj, sort_keys=True, indent=4, ensure_ascii=True, cls=VertexClassEncoder
        )
        return sJson

    def set_vertex_id(self, lVertexID: int) -> None:
        """
        Set the vertex ID.

        Args:
            lVertexID (int): New vertex ID

        Raises:
            TypeError: If lVertexID is not an integer
        """
        if not isinstance(lVertexID, (int, np.integer)):
            raise TypeError(f"Vertex ID must be an integer, got {type(lVertexID)}")
        self.lVertexID = int(lVertexID)

    def set_vertex_index(self, lVertexIndex: int) -> None:
        """
        Set the vertex array index.

        Args:
            lVertexIndex (int): New vertex index

        Raises:
            TypeError: If lVertexIndex is not an integer
        """
        if not isinstance(lVertexIndex, (int, np.integer)):
            raise TypeError(
                f"Vertex index must be an integer, got {type(lVertexIndex)}"
            )
        self.lVertexIndex = int(lVertexIndex)

    def set_flowline_id(self, lFlowlineID: int) -> None:
        """
        Set the associated flowline ID.

        Args:
            lFlowlineID (int): New flowline ID

        Raises:
            TypeError: If lFlowlineID is not an integer
        """
        if not isinstance(lFlowlineID, (int, np.integer)):
            raise TypeError(f"Flowline ID must be an integer, got {type(lFlowlineID)}")
        self.lFlowlineID = int(lFlowlineID)

    def is_valid(self) -> bool:
        """
        Check if the vertex has valid attributes.

        Returns:
            bool: True if vertex has valid coordinates and at least one valid ID
        """
        has_valid_coords = (
            -180 <= self.dLongitude_degree <= 180 and -90 <= self.dLatitude_degree <= 90
        )
        has_valid_id = self.lVertexID > 0 or self.lVertexIndex >= 0
        return has_valid_coords and has_valid_id

    def copy(self) -> "pyvertex":
        """
        Create a deep copy of the vertex.

        Returns:
            pyvertex: A new vertex object with the same attributes
        """
        param = {
            "dLongitude_degree": self.dLongitude_degree,
            "dLatitude_degree": self.dLatitude_degree,
            "x": self.dX_meter,
            "y": self.dY_meter,
            "z": self.dZ_meter,
            "dElevation": self.dElevation,
            "lVertexIndex": self.lVertexIndex,
            "lVertexID": self.lVertexID,
            "lFlowlineID": self.lFlowlineID,
        }
        return pyvertex(param)
