cimport cython
from libcpp.vector cimport vector
from libcpp.algorithm cimport sort
from libcpp.utility cimport pair
from libc.math cimport M_PI, sin, cos, asin,acos, sqrt, abs
from tinyr import RTree
from cython.operator cimport dereference as deref

""" Low-level function for pyflowline
"""
# Authors: Chang Liao

#constant

cdef double dRadius = 6378137.0

@cython.boundscheck(False)  # deactivate bnds checking
cpdef  calculate_distance_based_on_longitude_latitude(double dLongitude_degree1_in, double dLatitude_degree1_in, double dLongitude_degree2_in, double dLatitude_degree2_in):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians

    cdef double dLongitude_radian1_in, dLatitude_radian1_in
    cdef double dLongitude_radian2_in, dLatitude_radian2_in
    cdef double dLongtitude_diff, dLatitude_diff, a, b, c

    # Radius of earth in kilometers. Use 3956 for miles

    dLongitude_radian1_in = dLongitude_degree1_in /180.0 * M_PI
    dLatitude_radian1_in = dLatitude_degree1_in /180.0 * M_PI
    dLongitude_radian2_in = dLongitude_degree2_in /180.0 * M_PI
    dLatitude_radian2_in = dLatitude_degree2_in /180.0 * M_PI

    # haversine formula
    dLongtitude_diff = dLongitude_radian2_in - dLongitude_radian1_in
    dLatitude_diff = dLatitude_radian2_in - dLatitude_radian1_in
    a = sin(dLatitude_diff/2)*sin(dLatitude_diff/2) + cos(dLatitude_radian1_in) * cos(dLatitude_radian2_in) * sin(dLongtitude_diff/2)*sin(dLongtitude_diff/2)
    b = 2 * asin(sqrt(a))
    c = b * dRadius
    return c

@cython.boundscheck(False)  # deactivate bnds checking
@cython.wraparound(False)  # deactivate negative indexing
cpdef calculate_distance_based_on_longitude_latitude_array(double[:] dLongitude_degree1_in, double[:] dLatitude_degree1_in,
                                                         double[:] dLongitude_degree2_in, double[:] dLatitude_degree2_in):
    """
    Calculate the great circle distance between arrays of points
    on the earth (specified in decimal degrees)

    Args:
        dLongitude_degree1_in: Array of longitude values for first points (in degrees)
        dLatitude_degree1_in: Array of latitude values for first points (in degrees)
        dLongitude_degree2_in: Array of longitude values for second points (in degrees)
        dLatitude_degree2_in: Array of latitude values for second points (in degrees)

    Returns:
        numpy.ndarray: Array of distances in meters
    """
    import numpy as np

    cdef int i
    cdef int n = dLongitude_degree1_in.shape[0]
    cdef double[:] result = np.zeros(n, dtype=np.float64)

    cdef double dLongitude_radian1, dLatitude_radian1
    cdef double dLongitude_radian2, dLatitude_radian2
    cdef double dLongtitude_diff, dLatitude_diff, a, b, c

    for i in range(n):
        # Convert decimal degrees to radians
        dLongitude_radian1 = dLongitude_degree1_in[i] / 180.0 * M_PI
        dLatitude_radian1 = dLatitude_degree1_in[i] / 180.0 * M_PI
        dLongitude_radian2 = dLongitude_degree2_in[i] / 180.0 * M_PI
        dLatitude_radian2 = dLatitude_degree2_in[i] / 180.0 * M_PI

        # Haversine formula
        dLongtitude_diff = dLongitude_radian2 - dLongitude_radian1
        dLatitude_diff = dLatitude_radian2 - dLatitude_radian1
        a = sin(dLatitude_diff/2)*sin(dLatitude_diff/2) + cos(dLatitude_radian1) * cos(dLatitude_radian2) * sin(dLongtitude_diff/2)*sin(dLongtitude_diff/2)
        b = 2 * asin(sqrt(a))
        c = b * dRadius

        result[i] = c

    return np.asarray(result)

# For compatibility, also provide a version that works with NumPy arrays directly
@cython.boundscheck(False)
@cython.wraparound(False)
def calculate_distance_based_on_longitude_latitude_numpy(dLongitude_degree1_in, dLatitude_degree1_in,
                                                      dLongitude_degree2_in, dLatitude_degree2_in):
    """
    Calculate the great circle distance between arrays of points using NumPy arrays

    Args:
        dLongitude_degree1_in: NumPy array of longitude values for first points (in degrees)
        dLatitude_degree1_in: NumPy array of latitude values for first points (in degrees)
        dLongitude_degree2_in: NumPy array of longitude values for second points (in degrees)
        dLatitude_degree2_in: NumPy array of latitude values for second points (in degrees)

    Returns:
        numpy.ndarray: Array of distances in meters
    """
    import numpy as np

    # Convert decimal degrees to radians
    dLongitude_radian1_in = dLongitude_degree1_in / 180.0 * M_PI
    dLatitude_radian1_in = dLatitude_degree1_in / 180.0 * M_PI
    dLongitude_radian2_in = dLongitude_degree2_in / 180.0 * M_PI
    dLatitude_radian2_in = dLatitude_degree2_in / 180.0 * M_PI

    # Haversine formula
    dLongtitude_diff = dLongitude_radian2_in - dLongitude_radian1_in
    dLatitude_diff = dLatitude_radian2_in - dLatitude_radian1_in
    a = np.sin(dLatitude_diff/2)**2 + np.cos(dLatitude_radian1_in) * np.cos(dLatitude_radian2_in) * np.sin(dLongtitude_diff/2)**2
    c = 2 * np.arcsin(np.sqrt(a)) * dRadius

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
    cdef double x, y, left, right, bottom, top
    nVertex= len(aVertex_in)

    if nVertex > 0 :
        index_vertex = RTree( max_cap=5, min_cap=2)
        for i in range(nVertex):
            x = aVertex_in[i].dLongitude_degree
            y = aVertex_in[i].dLatitude_degree
            left =   x - 1E-5
            right =  x + 1E-5
            bottom = y - 1E-5
            top =    y + 1E-5
            pBound= (left, bottom, right, top)
            index_vertex.insert(i, pBound)

        aIntersect = list(index_vertex.search_surrounding([pVertex_in.dLongitude_degree, pVertex_in.dLatitude_degree]))

        for k in aIntersect:
            pVertex = aVertex_in[k]
            dDistance = pVertex.calculate_distance(pVertex_in)
            #if pVertex == pVertex_in:
            if dDistance < dThreshold_in:
                iFlag_exist = 1
                lIndex = k
                break
            else:
                pass
        pass

    return iFlag_exist, lIndex

cdef int compare_pairs(pair[double, int] a, pair[double, int] b):
    return a.first < b.first

@cython.boundscheck(False)  # deactivate bnds checking
cpdef find_vertex_on_edge(list aVertex_in, pEdge_in):
    #
    cdef int iFlag_exist = 0
    cdef int nVertex, npoint
    cdef vector[int] aIndex, aIndex_order
    cdef vector[double] aDistance
    cdef double x, y, left, right, bottom, top
    cdef vector[pair[double, int]] distance_index_pairs

    nVertex= len(aVertex_in)
    npoint = 0
    if nVertex > 0 :
        index_vertex = RTree(max_cap=5, min_cap=2)
        for i in range(nVertex):
            lID = i
            x = aVertex_in[i].dLongitude_degree
            y = aVertex_in[i].dLatitude_degree
            left = x - 1E-5
            right= x + 1E-5
            bottom= y -1E-5
            top=    y + 1E-5
            pBound= (left, bottom, right, top)
            index_vertex.insert(lID, pBound)  #
            pass
        #now the new vertex
        pVertex_start = pEdge_in.pVertex_start
        pVertex_end = pEdge_in.pVertex_end
        x1=pVertex_start.dLongitude_degree
        y1=pVertex_start.dLatitude_degree
        x2=pVertex_end.dLongitude_degree
        y2=pVertex_end.dLatitude_degree
        left = min(x1, x2)
        right = max(x1, x2)
        bottom = min(y1, y2)
        top = max(y1, y2)
        pBound= (left, bottom, right, top)
        aIntersect = list(index_vertex.search(pBound))
        for k in aIntersect:
            pVertex = aVertex_in[k]
            iFlag_overlap, dDistance, diff = pEdge_in.check_vertex_on_edge(pVertex)
            if iFlag_overlap == 1:
                iFlag_exist = 1
                aDistance.push_back(dDistance)
                aIndex.push_back(k)
                npoint = npoint + 1
            else:
                if diff < 1.0:
                    iFlag_overlap = pEdge_in.check_vertex_on_edge(pVertex)
                pass

        #re-order
        if iFlag_exist == 1 :
            # Create a vector of pairs (distance, index)
            for i in range(aDistance.size()):
                distance_index_pairs.push_back((aDistance[i], aIndex[i]))

            # Sort the vector of pairs by the first element (distance)
            sort(distance_index_pairs.begin(), distance_index_pairs.end(), compare_pairs)

            # Extract the sorted indices into a Python list
            aIndex_order = [distance_index_pairs[i].second for i in range(distance_index_pairs.size())]
            pass

    else:
        pass

    return iFlag_exist, npoint, aIndex_order

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
    if pVertex_in is None:
        raise ValueError("Input vertex is None.")

    iFlag_exist, dummy =  find_vertex_in_list(aVertex_in, pVertex_in, dThreshold_in)

    if iFlag_exist == 1:
        pass
    else:
        aVertex_in.append(pVertex_in)
        pass

    return aVertex_in, iFlag_exist

@cython.boundscheck(False)  # deactivate bnds checking
cpdef angle_between_vectors_coordinates(double x1, double y1, double z1, double x2, double y2, double z2):
    """Return the angle between two vectors in any dimension space,
    in degrees.
    """
    cdef double a, b, c, d, e, f
    a = x1*x2 + y1*y2 + z1*z2
    b = sqrt( x1*x1 + y1*y1 + z1*z1   )
    c = sqrt( x2*x2 + y2*y2 + z2*z2   )
    if b == 0 or c == 0:
        raise ValueError("Zero-length vector encountered. Cannot calculate angle.")
    d = a / (b* c)
    if d > 1:
        d = 1
    if d < -1:
        d = -1
    e = acos(d)
    f = e / M_PI * 180.0
    return f

@cython.boundscheck(False)  # deactivate bnds checking
cpdef (double, double, double) longlat_to_3d(dLongitude_degree_in, dLatitude_degree_in):
    """Convert a point given latitude and longitude in radians to
    3-dimensional space, assuming a sphere radius of one."""

    #cdef double a, b, c
    cdef double x, y, z
    cdef double dLongitude_degree, dLatitude_degree

    dLongitude_radian =  dLongitude_degree_in / 180.0 * M_PI
    dLatitude_radian = dLatitude_degree_in / 180.0 * M_PI

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

@cython.boundscheck(False)  # Disable bounds checking for performance
@cython.wraparound(False)   # Disable negative indexing for performance
def calculate_distance_to_plane(double dLongitude_degree1_in, double dLatitude_degree1_in,
                                    double dLongitude_degree2_in, double dLatitude_degree2_in,
                                    double dLongitude_degree3_in, double dLatitude_degree3_in):
    """
    Calculate the distance of a point to a plane defined by three points in 3D space.
    """

    cdef double x1, y1, z1
    cdef double x2, y2, z2
    cdef double x3, y3, z3
    cdef double v1_x, v1_y, v1_z
    cdef double v2_x, v2_y, v2_z
    cdef double normal_x, normal_y, normal_z
    cdef double A, B, C, D
    cdef double distance

    # Convert the three points to 3D coordinates
    x1, y1, z1 = longlat_to_3d(dLongitude_degree1_in, dLatitude_degree1_in)
    x2, y2, z2 = longlat_to_3d(dLongitude_degree2_in, dLatitude_degree2_in)
    x3, y3, z3 = longlat_to_3d(dLongitude_degree3_in, dLatitude_degree3_in)

    # Calculate two vectors on the plane
    v1_x = x2 - x1
    v1_y = y2 - y1
    v1_z = z2 - z1

    v2_x = x3 - x1
    v2_y = y3 - y1
    v2_z = z3 - z1

    # Compute the normal vector using the cross product
    normal_x = v1_y * v2_z - v1_z * v2_y
    normal_y = v1_z * v2_x - v1_x * v2_z
    normal_z = v1_x * v2_y - v1_y * v2_x

    # Check if the normal vector is zero (points are collinear)
    if abs(normal_x) < 1e-10 and abs(normal_y) < 1e-10 and abs(normal_z) < 1e-10:
        #raise ValueError("The three points are collinear in 3D space. A plane cannot be defined.")
        return 0.0  # Return zero distance if points are collinear

    # Calculate the plane equation coefficients (A, B, C, D)
    A = normal_x
    B = normal_y
    C = normal_z
    D = -(A * x1 + B * y1 + C * z1)

    # Calculate the distance of the second point to the plane
    distance = abs(A * x2 + B * y2 + C * z2 + D) / sqrt(A**2 + B**2 + C**2)

    return distance

@cython.boundscheck(False)  # deactivate bnds checking
cpdef calculate_distance_to_plane_old(dLongitude_degree1_in, dLatitude_degree1_in, dLongitude_degree2_in, dLatitude_degree2_in, dLongitude_degree3_in, dLatitude_degree3_in):
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
    # Check for zero denominators
    denominator_c = z1 * y3 - z3 * y1
    denominator_b = y1 * z3 - y3 * z1

    if denominator_c == 0 or denominator_b == 0:
        print(z1, y3, z3, y1)
        print(y1, z3, y3, z1)
        print(dLongitude_degree1_in, dLatitude_degree1_in)
        print(dLongitude_degree2_in, dLatitude_degree2_in)
        print(dLongitude_degree3_in, dLatitude_degree3_in)
        raise ValueError("Division by zero encountered in the calculation of coefficients 'b' and 'c'. Check input points.")

    # Calculate coefficients
    c = (-x1 * y3 + x3 * y1) / denominator_c
    b = (-x1 * z3 + x3 * z1) / denominator_b
    distance = abs(  x2 + b * y2 + c * z2 )
    return distance

