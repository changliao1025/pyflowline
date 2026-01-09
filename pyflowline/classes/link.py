"""
PyCellLink class for representing connections between mesh cells in flowline networks.

This module defines the pycelllink class which represents the topological
connections between mesh cells through shared edges. Cell links are essential
for defining the connectivity and adjacency relationships in mesh networks.
"""

import json
from json import JSONEncoder
from typing import Optional, Any, Tuple
import numpy as np

from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.edge import pyedge
from pyflowline.classes.flowline import pyflowline
from pyflowline.classes.meshcell import pymeshcell


class LinkClassEncoder(JSONEncoder):
    """
    Custom JSON encoder for pycelllink objects.

    Handles numpy data types, mesh cell objects, and other complex types,
    converting them to native Python types for JSON serialization.
    """

    def default(self, obj):
        if isinstance(obj, (np.integer, np.int32, np.int64)):
            return int(obj)
        if isinstance(obj, (np.float32, np.float64)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, list):
            return obj
        if pyedge and isinstance(obj, pyedge):
            return getattr(obj, "lEdgeID", str(obj))
        if pyvertex and isinstance(obj, pyvertex):
            return getattr(obj, "lVertexID", str(obj))
        if pyflowline and isinstance(obj, pyflowline):
            return getattr(obj, "lFlowlineID", str(obj))
        if pymeshcell and isinstance(obj, pymeshcell):
            return getattr(obj, "lCellID", str(obj))
        if isinstance(obj, pycelllink):
            return obj.lLinkID
        return JSONEncoder.default(self, obj)


class pycelllink:
    """
    Cell link class for representing connections between mesh cells.

    Represents the topological connection between two mesh cells through a
    shared edge. Cell links are fundamental for defining mesh connectivity,
    neighbor relationships, and flow routing in mesh-based models.

    A cell link connects exactly two mesh cells (start and end) and is
    associated with the edge that forms the boundary between them. This
    representation is essential for:
    - Mesh topology analysis
    - Flow routing between cells
    - Neighbor traversal algorithms
    - Network connectivity validation

    Attributes:
        lLinkIndex (int): Sequential index of the link in a collection
        lLinkID (int): Unique identifier for the link
        pCell_start (pymeshcell): Starting/source mesh cell
        pCell_end (pymeshcell): Ending/target mesh cell
        pEdge_link (pyedge): Edge object that connects the two cells
        dLink_length (float): Length of the connecting edge
        dLink_width (float): Width of the connection (if applicable)
        iLink_type (int): Type of link (internal, boundary, etc.)
        bLink_active (bool): Whether the link is active/valid

    Args:
        pCell_start_in (pymeshcell): The starting/source mesh cell object
        pCell_end_in (pymeshcell): The ending/target mesh cell object
        pEdge_link_in (pyedge): The edge object that links the two cells

    Raises:
        TypeError: If input parameters are not of the expected types
        ValueError: If cells or edge are None, or if cells are the same

    Example:
        >>> cell1 = pysquare(lon1, lat1, edges1, vertices1)
        >>> cell2 = pysquare(lon2, lat2, edges2, vertices2)
        >>> shared_edge = pyedge(vertex1, vertex2)
        >>> link = pycelllink(cell1, cell2, shared_edge)
        >>> link.set_link_id(100)
        >>> print(link)
        pycelllink(ID=100, Start=1, End=2, Length=1234.56m)
    """

    def __init__(self, pCell_start_in, pCell_end_in, pEdge_link_in) -> None:
        """
        Initialize a cell link object.

        Creates a connection between two mesh cells through a shared edge.
        The link represents the topological relationship and can store
        additional properties like length, type, and activity status.

        Args:
            pCell_start_in: The starting/source mesh cell object
            pCell_end_in: The ending/target mesh cell object
            pEdge_link_in: The edge object that connects the two cells

        Raises:
            TypeError: If input parameters are not of the expected types
            ValueError: If cells or edge are None, or if cells are identical
        """
        # Input validation
        if pCell_start_in is None:
            raise ValueError("Starting cell cannot be None")
        if pCell_end_in is None:
            raise ValueError("Ending cell cannot be None")
        if pEdge_link_in is None:
            raise ValueError("Linking edge cannot be None")

        # Check that cells are different (a cell shouldn't link to itself)
        if (
            hasattr(pCell_start_in, "lCellID")
            and hasattr(pCell_end_in, "lCellID")
            and pCell_start_in.lCellID == pCell_end_in.lCellID
            and pCell_start_in.lCellID != -1
        ):
            raise ValueError("Cannot create link: start and end cells are the same")

        # Validate cell types if possible
        if pymeshcell:
            if not isinstance(pCell_start_in, pymeshcell):
                raise TypeError(
                    f"Starting cell must be a mesh cell type, got {type(pCell_start_in)}"
                )
            if not isinstance(pCell_end_in, pymeshcell):
                raise TypeError(
                    f"Ending cell must be a mesh cell type, got {type(pCell_end_in)}"
                )

        # Initialize link attributes
        self.lLinkIndex: int = 0
        self.lLinkID: int = 0

        # Core link components
        self.pCell_start = pCell_start_in
        self.pCell_end = pCell_end_in
        self.pEdge_link = pEdge_link_in

        # Additional link properties
        self.dLink_length: float = -1.0  # Length of the connecting edge
        self.dLink_width: float = -1.0  # Width of the connection
        self.iLink_type: int = 0  # Type of link (0=internal, 1=boundary, etc.)
        self.bLink_active: bool = True  # Whether the link is active

        # Calculate link properties if possible
        self._calculate_link_properties()

    def __repr__(self) -> str:
        """
        Return a detailed string representation of the cell link.

        Returns:
            str: Detailed representation including link ID, cell IDs, and properties
        """
        start_id = (
            getattr(self.pCell_start, "lCellID", "Unknown")
            if self.pCell_start
            else "None"
        )
        end_id = (
            getattr(self.pCell_end, "lCellID", "Unknown") if self.pCell_end else "None"
        )

        return (
            f"pycelllink(ID={self.lLinkID}, Index={self.lLinkIndex}, "
            f"Start={start_id}, End={end_id}, Length={self.dLink_length:.2f}m, "
            f"Type={self.iLink_type}, Active={self.bLink_active})"
        )

    def __str__(self) -> str:
        """
        Return a concise string representation of the cell link.

        Returns:
            str: Concise representation with link ID and connected cell IDs
        """
        start_id = (
            getattr(self.pCell_start, "lCellID", "Unknown")
            if self.pCell_start
            else "None"
        )
        end_id = (
            getattr(self.pCell_end, "lCellID", "Unknown") if self.pCell_end else "None"
        )

        return (
            f"pycelllink(ID={self.lLinkID}, Start={start_id}, End={end_id}, "
            f"Length={self.dLink_length:.2f}m)"
        )

    def __hash__(self) -> int:
        """
        Return hash value for the cell link.

        Uses the combination of start cell, end cell, and edge for hashing.
        This allows links to be used in sets and as dictionary keys.

        Returns:
            int: Hash value based on link components
        """
        start_hash = hash(id(self.pCell_start)) if self.pCell_start else 0
        end_hash = hash(id(self.pCell_end)) if self.pCell_end else 0
        edge_hash = hash(id(self.pEdge_link)) if self.pEdge_link else 0

        return hash((start_hash, end_hash, edge_hash))

    def __eq__(self, other: Any) -> bool:
        """
        Check if two cell links are equal.

        Links are considered equal if they connect the same cells through
        the same edge (in either direction).

        Args:
            other: Another object to compare with

        Returns:
            bool: True if links are equivalent
        """
        if not isinstance(other, pycelllink):
            return NotImplemented

        # Check if same cells and edge (forward direction)
        same_forward = (
            self.pCell_start == other.pCell_start
            and self.pCell_end == other.pCell_end
            and self.pEdge_link == other.pEdge_link
        )

        # Check if same cells and edge (reverse direction)
        same_reverse = (
            self.pCell_start == other.pCell_end
            and self.pCell_end == other.pCell_start
            and self.pEdge_link == other.pEdge_link
        )

        return same_forward or same_reverse

    def _calculate_link_properties(self) -> None:
        """
        Calculate link properties such as length and width.

        Attempts to extract geometric properties from the connecting edge
        and adjacent cells.
        """
        try:
            # Calculate link length from edge
            if self.pEdge_link and hasattr(self.pEdge_link, "calculate_length"):
                self.dLink_length = self.pEdge_link.calculate_length()
            elif self.pEdge_link and hasattr(self.pEdge_link, "dLength"):
                self.dLink_length = self.pEdge_link.dLength

            # Try to determine link type based on cell properties
            if (
                self.pCell_start
                and self.pCell_end
                and hasattr(self.pCell_start, "iFlag_coast")
                and hasattr(self.pCell_end, "iFlag_coast")
            ):

                # Check if this is a boundary link (connects to coast/boundary)
                if self.pCell_start.iFlag_coast == 1 or self.pCell_end.iFlag_coast == 1:
                    self.iLink_type = 1  # Boundary link
                else:
                    self.iLink_type = 0  # Internal link

        except Exception:
            # If calculation fails, keep default values
            pass

    def set_link_id(self, lLinkID: int) -> None:
        """
        Set the unique identifier for the link.

        Args:
            lLinkID (int): New link ID

        Raises:
            TypeError: If lLinkID is not an integer
        """
        if not isinstance(lLinkID, (int, np.integer)):
            raise TypeError(f"Link ID must be an integer, got {type(lLinkID)}")

        self.lLinkID = int(lLinkID)

    def set_link_index(self, lLinkIndex: int) -> None:
        """
        Set the sequential index for the link.

        Args:
            lLinkIndex (int): New link index

        Raises:
            TypeError: If lLinkIndex is not an integer
        """
        if not isinstance(lLinkIndex, (int, np.integer)):
            raise TypeError(f"Link index must be an integer, got {type(lLinkIndex)}")

        self.lLinkIndex = int(lLinkIndex)

    def set_link_type(self, iLink_type: int) -> None:
        """
        Set the type classification for the link.

        Args:
            iLink_type (int): Link type (0=internal, 1=boundary, 2=special, etc.)

        Raises:
            TypeError: If iLink_type is not an integer
            ValueError: If iLink_type is negative
        """
        if not isinstance(iLink_type, (int, np.integer)):
            raise TypeError(f"Link type must be an integer, got {type(iLink_type)}")
        if iLink_type < 0:
            raise ValueError(f"Link type must be non-negative, got {iLink_type}")

        self.iLink_type = int(iLink_type)

    def set_active(self, bActive: bool) -> None:
        """
        Set the active status of the link.

        Args:
            bActive (bool): True to activate link, False to deactivate

        Raises:
            TypeError: If bActive is not a boolean
        """
        if not isinstance(bActive, bool):
            raise TypeError(f"Active status must be a boolean, got {type(bActive)}")

        self.bLink_active = bActive

    def reverse_direction(self) -> None:
        """
        Reverse the direction of the link.

        Swaps the start and end cells, effectively reversing the link direction
        while maintaining the same connection.
        """
        temp_cell = self.pCell_start
        self.pCell_start = self.pCell_end
        self.pCell_end = temp_cell

    def get_connected_cells(self) -> Tuple[Any, Any]:
        """
        Get the two cells connected by this link.

        Returns:
            Tuple[Any, Any]: (start_cell, end_cell)
        """
        return (self.pCell_start, self.pCell_end)

    def get_other_cell(self, cell) -> Optional[Any]:
        """
        Get the cell on the other end of the link from a given cell.

        Args:
            cell: The reference cell

        Returns:
            Optional[Any]: The other connected cell, or None if reference cell not found

        Raises:
            ValueError: If cell is None
        """
        if cell is None:
            raise ValueError("Reference cell cannot be None")

        if self.pCell_start == cell:
            return self.pCell_end
        elif self.pCell_end == cell:
            return self.pCell_start
        else:
            return None

    def contains_cell(self, cell) -> bool:
        """
        Check if the link connects to a specific cell.

        Args:
            cell: The cell to check for

        Returns:
            bool: True if the cell is connected by this link

        Raises:
            ValueError: If cell is None
        """
        if cell is None:
            raise ValueError("Cell to check cannot be None")

        return self.pCell_start == cell or self.pCell_end == cell

    def get_link_properties(self) -> dict:
        """
        Get comprehensive properties of the link.

        Returns:
            dict: Dictionary containing all link properties
        """
        properties = {
            "link_id": self.lLinkID,
            "link_index": self.lLinkIndex,
            "link_type": self.iLink_type,
            "active": self.bLink_active,
            "length": self.dLink_length,
            "width": self.dLink_width,
            "start_cell_id": (
                getattr(self.pCell_start, "lCellID", None) if self.pCell_start else None
            ),
            "end_cell_id": (
                getattr(self.pCell_end, "lCellID", None) if self.pCell_end else None
            ),
            "edge_id": (
                getattr(self.pEdge_link, "lEdgeID", None) if self.pEdge_link else None
            ),
        }

        return properties

    def is_boundary_link(self) -> bool:
        """
        Check if this is a boundary link.

        A boundary link connects to cells that are on the mesh boundary
        or have special boundary properties.

        Returns:
            bool: True if this is a boundary link
        """
        return self.iLink_type == 1

    def is_valid(self) -> bool:
        """
        Check if the link has valid properties.

        A link is considered valid if:
        - Both connected cells exist
        - The connecting edge exists
        - The link is active
        - Geometric properties are reasonable

        Returns:
            bool: True if link is valid
        """
        has_valid_cells = self.pCell_start is not None and self.pCell_end is not None
        has_valid_edge = self.pEdge_link is not None
        has_reasonable_length = self.dLink_length >= 0

        return (
            has_valid_cells
            and has_valid_edge
            and has_reasonable_length
            and self.bLink_active
        )

    def calculate_flow_direction(self) -> int:
        """
        Determine flow direction based on cell properties.

        Uses cell elevation or other properties to determine the likely
        flow direction across the link.

        Returns:
            int: Flow direction (1=start to end, -1=end to start, 0=unknown)
        """
        try:
            # Try to use elevation data if available
            if hasattr(self.pCell_start, "dElevation_mean") and hasattr(
                self.pCell_end, "dElevation_mean"
            ):

                elev_start = self.pCell_start.dElevation_mean
                elev_end = self.pCell_end.dElevation_mean

                if elev_start > elev_end + 1e-6:  # Small tolerance
                    return 1  # Flow from start to end (downhill)
                elif elev_end > elev_start + 1e-6:
                    return -1  # Flow from end to start (downhill)
                else:
                    return 0  # Same elevation or unknown

            # Try to use stream order if available
            if hasattr(self.pCell_start, "iStream_order_burned") and hasattr(
                self.pCell_end, "iStream_order_burned"
            ):

                order_start = self.pCell_start.iStream_order_burned
                order_end = self.pCell_end.iStream_order_burned

                if order_start < order_end:
                    return 1  # Flow from lower to higher order
                elif order_end < order_start:
                    return -1  # Flow from lower to higher order

        except Exception:
            pass

        return 0  # Unknown direction

    def copy(self) -> "pycelllink":
        """
        Create a deep copy of the cell link.

        Returns:
            pycelllink: A new cell link object with the same properties
        """
        new_link = pycelllink(self.pCell_start, self.pCell_end, self.pEdge_link)

        # Copy all attributes
        new_link.lLinkIndex = self.lLinkIndex
        new_link.lLinkID = self.lLinkID
        new_link.dLink_length = self.dLink_length
        new_link.dLink_width = self.dLink_width
        new_link.iLink_type = self.iLink_type
        new_link.bLink_active = self.bLink_active

        return new_link

    def tojson(self) -> str:
        """
        Convert cell link object to a JSON string.

        Serializes all link attributes including connected cells and edge
        information. Uses the custom LinkClassEncoder to handle complex objects.

        Returns:
            str: JSON string representation of the cell link

        Example:
            >>> link = pycelllink(cell1, cell2, edge)
            >>> json_str = link.tojson()
        """
        obj = self.__dict__.copy()

        # Create a clean copy for serialization
        obj_clean = {}
        for key, value in obj.items():
            if key in ["pCell_start", "pCell_end", "pEdge_link"]:
                # Replace with IDs for serialization
                if hasattr(value, "lCellID"):
                    obj_clean[key + "_id"] = value.lCellID
                elif hasattr(value, "lEdgeID"):
                    obj_clean[key + "_id"] = value.lEdgeID
                else:
                    obj_clean[key + "_id"] = str(value) if value else None
            else:
                obj_clean[key] = value

        sJson = json.dumps(
            obj_clean, sort_keys=True, indent=4, ensure_ascii=True, cls=LinkClassEncoder
        )
        return sJson
