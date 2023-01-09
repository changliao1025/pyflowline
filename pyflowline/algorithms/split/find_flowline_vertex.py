
#from pyflowline.algorithms.auxiliary.find_index_in_list import add_unique_vertex
from pyflowline.algorithms.cython.kernel import add_unique_vertex
def find_flowline_vertex(aFlowline_in, dThreshold_in=1.0E-6): 
    nFlowline = len(aFlowline_in) 
    aVertex=list()
    for i in range(0, nFlowline):      
        pFlowline = aFlowline_in[i]
        pVertex_start = pFlowline.pVertex_start
        pVertex_end = pFlowline.pVertex_end
        aVertex, dummy = add_unique_vertex(aVertex, pVertex_start,dThreshold_in)
        aVertex, dummy = add_unique_vertex(aVertex, pVertex_end, dThreshold_in)
        pass

    return aVertex
