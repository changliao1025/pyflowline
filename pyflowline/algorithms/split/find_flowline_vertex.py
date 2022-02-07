
from pyflowline.algorithms.auxiliary.find_index_in_list import add_unique_vertex

def find_flowline_vertex(aFlowline_in):
   
    
    nFlowline = len(aFlowline_in)
 
    aVertex=list()
    for i in range(0, nFlowline):      
        pFlowline = aFlowline_in[i]
        pVertex_start = pFlowline.pVertex_start
        pVertex_end = pFlowline.pVertex_end
        aVertex, dummy = add_unique_vertex(aVertex, pVertex_start)
        aVertex, dummy = add_unique_vertex(aVertex, pVertex_end)
        pass


    aVertex_middle=list()
    for i in range(0, nFlowline):      
        pFlowline = aFlowline_in[i]
        aVertex2 = pFlowline.aVertex
        nVertex = pFlowline.nVertex
        for j in range(1, nVertex-1):   
            pVertex= aVertex2[j]
            aVertex_middle.append( pVertex)

    #aVertex_middle2 = list()
    #for elem in aVertex_middle:
    #    if aVertex_middle.count(elem) > 1:
    #        aVertex_middle2, dummy = add_unique_vertex(aVertex_middle2, elem)            
    #    else:
    #        pass

    #aVertex_middle2 =  [x for n, x in enumerate(aVertex_middle) if x in aVertex_middle[:n]] #[item for item, count in collections.Counter(aVertex_middle).items() if count > 1]
    
    #merge
    #aVertex = list(set().union(aVertex, aVertex_middle2))


    return aVertex
