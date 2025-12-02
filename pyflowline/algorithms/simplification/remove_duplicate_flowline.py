
"""
Improved duplicate removal for pyflowline objects based on geometry only.
"""

from pyflowline.algorithms.simplification.fast_remove_duplicates import remove_duplicate_flowlines_fast

def remove_duplicate_flowline(aFlowline_in):
    """
    Remove duplicate flowlines based purely on vertex geometry.

    This function ignores all IDs and focuses only on the geometric properties
    of flowlines by comparing all vertices in the flowline path.

    Args:
        aFlowline_in: List of pyflowline objects

    Returns:
        List of unique pyflowline objects (based on geometry only)
    """
    return remove_duplicate_flowlines_fast(aFlowline_in)

