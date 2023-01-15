import numpy as np
from math import cos, sin, sqrt, acos
def longlat_to_3d(lonr, latr):
    """Convert a point given latitude and longitude in radians to
    3-dimensional space, assuming a sphere radius of one."""
    a = cos(latr) * cos(lonr)
    b = cos(latr) * sin(lonr)
    c = sin(latr)
    return np.array((a,b,c))