import sys, os
import numpy as np


from pystream.check_head_water import check_head_water

def remove_small_river(aFlowline_in, dThreshold):

    nFlowline = len(aFlowline_in)
    aFlowline_out=list()        

    lID = 0

    
    
    for i in range(nFlowline):
        pFlowline = aFlowline_in[i]      
        pVertex_start = pFlowline.pVertex_start
        pVertex_end = pFlowline.pVertex_end
        
        dLength = pFlowline.calculate_length()
        
        if check_head_water(aFlowline_in, pVertex_start)==1:
            if dLength > dThreshold :
                pFlowline.lIndex = lID
                aFlowline_out.append(pFlowline)
                lID = lID + 1 
                pass
            else:
                #print('small')
                pass
        else:        
            pFlowline.lIndex = lID
            aFlowline_out.append(pFlowline)
            lID = lID +1       
            pass
        
        
        pass
    

    return aFlowline_out