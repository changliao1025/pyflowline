import numpy as np
from pystream.shared.vertex import pyvertex

from pystream.find_vertex_in_list import find_vertex_in_list
def add_unique_vertex(aVertex, pVertex_in):
    iFlag_exist = 0
    nVertex = len(aVertex)
    
    #for i in np.arange( nVertex):
    #    pVertex = aVertex[i]
    #    if pVertex == pVertex_in:
    #        iFlag_exist = 1
    #        pass
    #    else:
    #        pass
    #    pass        

    iFlag_exist, dummy =  find_vertex_in_list(aVertex, pVertex_in)

    if iFlag_exist == 1:
        pass
    else:
        #add it into the dic
        aVertex.append(pVertex_in)
        pass

    return aVertex, iFlag_exist

