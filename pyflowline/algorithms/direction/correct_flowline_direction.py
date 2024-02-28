import sys
import numpy as np 
from pyflowline.algorithms.auxiliary.check_head_water import check_head_water
sys.setrecursionlimit(100000)
lFlowlineIndex=0
def correct_flowline_direction(aFlowline_in, pVertex_outlet_in):    
    """_summary_ This function should expect the flowline may not be ordered, so the stream order info is not available.

    Args:
        aFlowline_in (_type_): _description_
        pVertex_outlet_in (_type_): _description_

    Returns:
        _type_: List of flowline, ordered from outlet to headwater
    """
    
    aFlowline_out= list()     
    #we have to go reversely    
    aFlag_process=None
    global lFlowlineIndex    
    nFlowline = len(aFlowline_in)
    # Create sets for faster lookup
    pVertex_start_in_set = {flowline.pVertex_start for flowline in aFlowline_in}    
    aFlag_process=np.full(nFlowline, 0, dtype =int)   
    dDiatance_min = float('inf')
    for i in range(nFlowline):        
        pFlowline = aFlowline_in[i]                
        pVertex_end = pFlowline.pVertex_end
        dDiatance = pVertex_end.calculate_distance(pVertex_outlet_in)
        if dDiatance < dDiatance_min:
            dDiatance_min = dDiatance
            lIndex_outlet = i                  
        
    lFlowlineIndex = 0 
    pFlowline = aFlowline_in[lIndex_outlet]
    pFlowline.lFlowlineIndex = lFlowlineIndex
    pFlowline.iFlag_dam = 1 # this one is not dam, but it should be preserved
    aFlowline_out.append(pFlowline)
    pVertex_start = pFlowline.pVertex_start
    pVertex_end = pFlowline.pVertex_end
    lFlowlineIndex = lFlowlineIndex + 1    
    #we might find more than 1 upstream    
    def find_upstream_flowline( pVertex_end_in):
        nUpstream = 0 
        aUpstream=list()
        aReverse=list()
        for i in range(nFlowline):
            pFlowline = aFlowline_in[i]
            pVerter_start = pFlowline.pVertex_start
            pVerter_end = pFlowline.pVertex_end
            #if pVerter_end == pVertex_start_in  and pVerter_start!=pVertex_end_in:
            if pVerter_end in pVertex_start_in_set and pVerter_start != pVertex_end_in:
                if aFlag_process[i] != 1:
                    nUpstream = nUpstream + 1
                    aUpstream.append(i)
                    aReverse.append(0)
                    aFlag_process[i] = 1
                    pass
            else:
                #if pVerter_start == pVertex_start_in and pVerter_end !=pVertex_end_in :
                if pVerter_start in pVertex_start_in_set and pVerter_end != pVertex_end_in:
                    if aFlag_process[i] != 1:
                        nUpstream = nUpstream + 1
                        aUpstream.append(i)
                        aReverse.append(1)
                        aFlag_process[i] = 1
                        pass
                pass

            pass
        return nUpstream, aUpstream, aReverse
    
    
    def tag_upstream(pVertex_start_in, pVertex_end_in):
        
        global lFlowlineIndex
        if(check_head_water(aFlowline_in, pVertex_start_in)==1):            
            pass
        else:
            nUpstream, aUpstream, aReverse = find_upstream_flowline( pVertex_end_in)
            if nUpstream > 0:                
                for j in range(nUpstream):
                    pFlowline = aFlowline_in[ aUpstream[j] ] 
                    if (aReverse[j]==1):             
                        pFlowline.reverse()                    
                        pass
                    else:
                        pass
                    pFlowline.lFlowlineIndex = lFlowlineIndex
                    aFlowline_out.append(pFlowline)
                    lFlowlineIndex = lFlowlineIndex + 1
                    tag_upstream(  pFlowline.pVertex_start, pFlowline.pVertex_end  )            

                pass
            else:
                pass

    tag_upstream(pVertex_start, pVertex_end)
    
    return aFlowline_out