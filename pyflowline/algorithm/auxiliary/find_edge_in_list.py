import numpy as np
def find_edge_in_list(aEdge, pEdge_in):
    iFlag_exist = 0
    lIndex= -1
    nEdge= len(aEdge)

    if nEdge > 0 :
        for i in np.arange( nEdge):
            pEdge = aEdge[i]
            if pEdge == pEdge_in:
                iFlag_exist = 1      
                lIndex = i 
                break                
            else:
                pass

        pass        
        
    else:
        pass
    
    return iFlag_exist, lIndex
