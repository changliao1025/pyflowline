import numpy as np
cimport cython
from libcpp.vector cimport vector
from libc.math cimport sin, cos, asin,acos, sqrt, abs
""" Low-level function for pyflowline
"""
# Authors: Chang Liao

#constant

cdef double pi = 3.14159265358979323846264338327
cdef double dRadius = 6378137.0 

@cython.boundscheck(False)  # deactivate bnds checking
cpdef  calculate_distance_based_on_lon_lat(double dLongitude_degree1_in, double dLatitude_degree1_in, double dLongitude_degree2_in, double dLatitude_degree2_in):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
   
    cdef double dLongitude_radian1_in, dLatitude_radian1_in
    cdef double dLongitude_radian2_in, dLatitude_radian2_in
    cdef double dLongtitude_diff, dLatitude_diff, a, b, c
   
    # Radius of earth in kilometers. Use 3956 for miles 
   
    dLongitude_radian1_in = dLongitude_degree1_in /180.0 * pi
    dLatitude_radian1_in = dLatitude_degree1_in /180.0 * pi
    dLongitude_radian2_in = dLongitude_degree2_in /180.0 * pi
    dLatitude_radian2_in = dLatitude_degree2_in /180.0 * pi
    
    # haversine formula 
    dLongtitude_diff = dLongitude_radian2_in - dLongitude_radian1_in 
    dLatitude_diff = dLatitude_radian2_in - dLatitude_radian1_in 
    a = sin(dLatitude_diff/2)*sin(dLatitude_diff/2) + cos(dLatitude_radian1_in) * cos(dLatitude_radian2_in) * sin(dLongtitude_diff/2)*sin(dLongtitude_diff/2)
    b = 2 * asin(sqrt(a)) 
    c = b * dRadius
    return c

@cython.boundscheck(False)  # deactivate bnds checking
cpdef  convert_360_to_180(double dLongitude_in):
    """[This function is modified from
    http://www.idlcoyote.com/map_tips/lonconvert.html]

    Args:
        dLongitude_in ([type]): [description]

    Returns:
        [type]: [description]
    """ 
    cdef int a 
    cdef double dLongitude_out
    a = int(dLongitude_in /180.0)
    dLongitude_out = dLongitude_in - a * 360.0
    return dLongitude_out

@cython.boundscheck(False)  # deactivate bnds checking
cpdef  find_vertex_in_list(list aVertex_in, pVertex_in, double dThreshold_in = 1.0E-6):
    """[find the index of a vertex in a list]

    Args:
        aVertex_in ([type]): [description]
        pVertex_in ([type]): [description]

    Returns:
        [type]: [description]
    """
    cdef int i
    cdef int iFlag_exist = 0
    cdef int lIndex= -1
    cdef int nVertex
    cdef double dDistance
    nVertex= len(aVertex_in)

    if nVertex > 0 :

        for i in range(nVertex):
            pVertex = aVertex_in[i]
            dDistance = pVertex.calculate_distance(pVertex_in)
            #if pVertex == pVertex_in:
            if dDistance < dThreshold_in:
                iFlag_exist = 1      
                lIndex = i 
                break                
            else:
                pass

        pass        
        
    else:
        pass
    
    return iFlag_exist, lIndex

@cython.boundscheck(False)  # deactivate bnds checking
cpdef find_vertex_on_edge(list aVertex_in, pEdge_in):
    #
    cdef int iFlag_exist = 0
    cdef int nVertex, npoint
    cdef vector[int] aIndex ,aIndex_order
    cdef vector[double] aDistance 
  
    nVertex= len(aVertex_in)
    npoint = 0    
    if nVertex > 0 :
        for i in range( nVertex):
            pVertex = aVertex_in[i]
            iFlag_overlap, dDistance, diff = pEdge_in.check_vertex_on_edge(pVertex)
            if iFlag_overlap == 1:                
                iFlag_exist = 1                   
                aDistance.push_back(dDistance)
                aIndex.push_back(i)
                npoint = npoint + 1          
            else:                
                if  diff < 1.0:
                    iFlag_overlap = pEdge_in.check_vertex_on_edge(pVertex) 
                pass

        #re-order 
        if iFlag_exist == 1 :
            x = np.array(aDistance)
            b = np.argsort(x)
            c = np.array(aIndex)
            d= c[b]
            aIndex_order = list(d)    

    else:
        pass
    
    return iFlag_exist, npoint , aIndex_order

@cython.boundscheck(False)  # deactivate bnds checking
cpdef add_unique_vertex(list aVertex_in, pVertex_in, double dThreshold_in = 1.0E-6):
    """[add a vertex to a list if it is not already included]

    Args:
        aVertex_in ([type]): [description]
        pVertex_in ([type]): [description]

    Returns:
        [type]: [description]
    """
    cdef int iFlag_exist
    cdef int nVertex
    cdef int dummy
    iFlag_exist = 0
    nVertex = len(aVertex_in)     

    iFlag_exist, dummy =  find_vertex_in_list(aVertex_in, pVertex_in, dThreshold_in)

    if iFlag_exist == 1:
        pass
    else:
        aVertex_in.append(pVertex_in)
        pass

    return aVertex_in, iFlag_exist


@cython.boundscheck(False)  # deactivate bnds checking
#cdef angle_between_vectors_coordinates(double *u, double *v):
cpdef angle_between_vectors_coordinates(double x1, double y1, double z1, double x2, double y2, double z2):
    """Return the angle between two vectors in any dimension space,
    in degrees.
    """
    cdef double a, b, c, d, e, f 
    a = x1*x2 + y1*y2 + z1*z2
    b = sqrt( x1*x1 + y1*y1 + z1*z1   ) 
    c = sqrt( x2*x2 + y2*y2 + z2*z2   ) 
    d = a / (b* c)
    if d > 1:
        d = 1
    if d < -1:
        d = -1
    e = acos(d)
    f = e / pi * 180.0
    return f 

@cython.boundscheck(False)  # deactivate bnds checking
cpdef (double, double, double) longlat_to_3d(dLongitude_degree_in, dLatitude_degree_in):
    """Convert a point given latitude and longitude in radians to
    3-dimensional space, assuming a sphere radius of one."""

    #cdef double a, b, c
    cdef double x, y, z
    cdef double dLongitude_degree, dLatitude_degree

    dLongitude_radian =  dLongitude_degree_in / 180.0 * pi
    dLatitude_radian = dLatitude_degree_in / 180.0 * pi

    x = cos(dLatitude_radian) * cos(dLongitude_radian)
    y = cos(dLatitude_radian) * sin(dLongitude_radian)
    z = sin(dLatitude_radian)
 
    return x, y, z


@cython.boundscheck(False)  # deactivate bnds checking
cpdef calculate_angle_betwen_vertex(dLongitude_degree1_in, dLatitude_degree1_in, dLongitude_degree2_in, dLatitude_degree2_in, dLongitude_degree3_in, dLatitude_degree3_in):
    #all in degree

    cdef double angle3deg
    cdef double x1, y1, z1
    cdef double x2, y2, z2 
    cdef double x3, y3, z3 
    cdef double x4, y4, z4
    cdef double x5, y5, z5
    
    # The points in 3D space
    x1, y1, z1 = longlat_to_3d(dLongitude_degree1_in, dLatitude_degree1_in)
    x2, y2, z2 = longlat_to_3d(dLongitude_degree2_in, dLatitude_degree2_in)
    x3, y3, z3 = longlat_to_3d(dLongitude_degree3_in, dLatitude_degree3_in)
    # Vectors in 3D space
    
    x4 = x1 - x2
    y4 = y1 - y2
    z4 = z1 - z2
    #c3vec[i] = aCoordinate3[i] - aCoordinate2[i]    
    x5 = x3 - x2
    y5 = y3 - y2
    z5 = z3 - z2

    angle3deg = angle_between_vectors_coordinates( x4, y4, z4, x5, y5, z5)
    return  angle3deg

@cython.boundscheck(False)  # deactivate bnds checking
cpdef calculate_distance_to_plane(dLongitude_degree1_in, dLatitude_degree1_in, dLongitude_degree2_in, dLatitude_degree2_in, dLongitude_degree3_in, dLatitude_degree3_in):    
    cdef double x1, y1, z1
    cdef double x2, y2, z2 
    cdef double x3, y3, z3 
    cdef double distance
    cdef b, c  
    # The points in 3D space
    x1, y1, z1 = longlat_to_3d(dLongitude_degree1_in, dLatitude_degree1_in)
    x2, y2, z2 = longlat_to_3d(dLongitude_degree2_in, dLatitude_degree2_in)
    x3, y3, z3 = longlat_to_3d(dLongitude_degree3_in, dLatitude_degree3_in)
    #The formula is x+b*y+c*z=0    
    c = (-x1*y3 + x3* y1)/( z1*y3 - z3*y1 )
    b = (-x1*z3 + x3 * z1 ) / (y1 * z3 - y3*z1)
    distance = abs(  x2 + b * y2 + c * z2 )
    return distance

