from rtree.index import Index as RTreeindex
def find_vertex_in_list(aVertex_in, pVertex_in, dThreshold_in=1.0E-6):
    """[find the index of a vertex in a list]

    Args:
        aVertex_in ([type]): [description]
        pVertex_in ([type]): [description]

    Returns:
        [type]: [description]
    """

    iFlag_exist = 0
    lIndex= -1
    nVertex= len(aVertex_in)
    index_vertex = RTreeindex()
    for i in range(nVertex):
        lID = i
        x=aVertex_in[i].dLongitude_degree
        y=aVertex_in[i].dLatitude_degree
        left   = x - 1E-5
        right  = x + 1E-5
        bottom = y - 1E-5
        top    = y + 1E-5
        pBound = (left, bottom, right, top)
        index_vertex.insert(lID, pBound)  #
        pass
    #now the new vertex
    x = pVertex_in.dLongitude_degree
    y = pVertex_in.dLatitude_degree
    delta = 1E-5  # or your desired search tolerance
    bbox = (x - delta, y - delta, x + delta, y + delta)
    aIntersect = list(index_vertex.intersection(bbox))
    #aIntersect = list(index_vertex.search_surrounding([pVertex_in.dLongitude_degree, pVertex_in.dLatitude_degree]))
    for k in aIntersect:
        pVertex = aVertex_in[k]
        #dDistance = pVertex.calculate_distance(pVertex_in)
        if pVertex == pVertex_in: #if dDistance < dThreshold_in:
            iFlag_exist = 1
            lIndex = k
            break
        else:
            pass
    

    return iFlag_exist, lIndex