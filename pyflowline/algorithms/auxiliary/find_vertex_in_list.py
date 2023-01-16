import numpy as np

def find_vertex_in_list(aVertex_in, pVertex_in):
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

    if nVertex > 0 :
        for i in np.arange( nVertex):
            pVertex = aVertex_in[i]
            if pVertex == pVertex_in:
                iFlag_exist = 1      
                lIndex = i 
                break                
            else:
                pass

        pass        
        
    else:
        pass
    
    return iFlag_exist, lIndex