cimport cython
from libc.math cimport sin, cos, arcsin, sqrt
""" Low-level function for pyflowline
"""
# Authors: Chang Liao

@cython.boundscheck(False)  # deactivate bnds checking
cpdef calculate_distance_based_on_lon_lat(dLon1_in, dLat1_in, dLon2_in, dLat2_in):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    cdef pi , r    
    cdef double lon1, lat1, lon2, lat2
    cdef double dlon, dlat, a, b, c
    #constant
    pi = 3.14159265358979323846264338327
    # Radius of earth in kilometers. Use 3956 for miles 
    r = 6378137.0 
    lon1 = dLon1_in /180.0 * pi
    lat1 = dLat1_in /180.0 * pi
    lon2 = dLon2_in /180.0 * pi
    lat2 = dLat2_in /180.0 * pi
    
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)*sin(dlat/2) + cos(lat1) * cos(lat2) * sin(dlon/2)*sin(dlon/2)
    b = 2 * arcsin(sqrt(a)) 
    c = b * r
    return c

@cython.boundscheck(False)  # deactivate bnds checking
cpdef convert_360_to_180(dLongitude_in):
    """[This function is modified from
    http://www.idlcoyote.com/map_tips/lonconvert.html]

    Args:
        dLongitude_in ([type]): [description]

    Returns:
        [type]: [description]
    """ 
    cdef a, dLongitude_out
    a = int(dLongitude_in /180.0)
    dLongitude_out = dLongitude_in - a * 360.0
    return dLongitude_out