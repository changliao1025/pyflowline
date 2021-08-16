import numpy as np
def find_flowline_in_list(aFlowline, pFlowline_in):
    iFlag_exist = 0
    lIndex= -1
    nFlowline= len(aFlowline)

    if nFlowline > 0 :
        for i in np.arange( nFlowline):
            pFlowline = aFlowline[i]
            if pFlowline == pFlowline_in:
                iFlag_exist = 1      
                lIndex = i 
                break                
            else:
                pass

        pass        
        
    else:
        pass
    
    return iFlag_exist, lIndex
