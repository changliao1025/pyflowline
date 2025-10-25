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

@cython.boundscheck(False)  # deactivate bnds checking
cpdef find_vertex_in_list(list aVertex_in, pVertex_in, double dThreshold_in = 1.0E-6):
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


