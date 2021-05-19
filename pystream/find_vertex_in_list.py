
import numpy as np

def find_vertex_in_list(aVertex, pVertex_in):

    iFlag_exist = 0
    lIndex= -1
    nVertex= len(aVertex)

    if nVertex > 0 :
        for i in np.arange( nVertex):
            pVertex = aVertex[i]
            if pVertex == pVertex_in:
                iFlag_exist = 1      
                lIndex = i 
                pass
            else:
                pass

        pass        
        
    else:
        pass
    
    return iFlag_exist, lIndex
   