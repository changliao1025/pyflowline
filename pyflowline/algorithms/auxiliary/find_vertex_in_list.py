
import importlib
import numpy as np
iFlag_cython = importlib.util.find_spec("cython") 
if iFlag_cython is not None:
    from pyflowline.external.tinyr.tinyr.tinyr import RTree
    iFlag_use_rtree = 1
else:
    iFlag_use_rtree =0
    pass

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

    if iFlag_use_rtree == 1:

        #can we use rtree here?
        #index_vertex = rtree.index.Index()
        index_vertex = RTree(max_cap=5, min_cap=2)
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
        x=pVertex_in.dLongitude_degree
        y=pVertex_in.dLatitude_degree   
        left = x - 1E-5
        right = x + 1E-5
        bottom = y-1E-5
        top =    y+1E-5
        pBound= (left, bottom, right, top)
        aIntersect = list(index_vertex.search(pBound))

        for k in aIntersect:
            pVertex = aVertex_in[k]
            #dDistance = pVertex.calculate_distance(pVertex_in)
            if pVertex == pVertex_in: #if dDistance < dThreshold_in:
            
                iFlag_exist = 1      
                lIndex = k
                break                
            else:
                pass
        pass
    else:

        if nVertex > 0 :
            for i in np.arange( nVertex):
                pVertex = aVertex_in[i]
                #dDistance = pVertex.calculate_distance(pVertex_in)
                if pVertex == pVertex_in: #if dDistance < dThreshold_in:
                
                    iFlag_exist = 1      
                    lIndex = i 
                    break                
                else:
                    pass

            pass        

        else:
            pass
    
    return iFlag_exist, lIndex