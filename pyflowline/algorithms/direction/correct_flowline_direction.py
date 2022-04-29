import numpy as np 
from pyflowline.algorithms.auxiliary.check_head_water import check_head_water
lID=0
aFlag_process=None

def correct_flowline_direction(aFlowline_in, pVertex_outlet_in):        
    #we have to go reversely    
    aFlowline_out= list()   
    global lID    
    global aFlag_process
    nFlowline = len(aFlowline_in)
    aFlag_process=np.full(nFlowline, 0, dtype =int)
    iFlag_first = 1
    for i in range(nFlowline):
        pFlowline = aFlowline_in[i]
        iFlag_dam = pFlowline.iFlag_dam
        iStream_order = pFlowline.iStream_order
        pVertex_start = pFlowline.pVertex_start
        pVertex_end = pFlowline.pVertex_end
        dDiatance = pVertex_end.calculate_distance( pVertex_outlet_in )   
        if iFlag_first ==1:
            dDiatance_min = dDiatance                
            lIndex_outlet = i            
            iFlag_first=0    
        else:            
            if  dDiatance < dDiatance_min:
                dDiatance_min = dDiatance
                #found it
                lIndex_outlet = i      
                pass    
            else:
          
                pass
        
    lID = 0 
    pFlowline = aFlowline_in[lIndex_outlet]
    pFlowline.lIndex = lID
    pFlowline.iFlag_dam = 1 # this one is not dam, but it should be preserved
    aFlowline_out.append(pFlowline)
    pVertex_start = pFlowline.pVertex_start
    pVertex_end = pFlowline.pVertex_end
    lID = lID + 1    

    #we might find more than 1 upstream    
    def find_upstream_flowline(pVertex_start_in, pVertex_end_in):
        nUpstream = 0 
        aUpstream=list()
        aReverse=list()
        for i in range(nFlowline):
            pFlowline = aFlowline_in[i]
            pVerter_start = pFlowline.pVertex_start
            pVerter_end = pFlowline.pVertex_end
            if pVerter_end == pVertex_start_in  and pVerter_start!=pVertex_end_in:
                if aFlag_process[i] != 1:
                    nUpstream = nUpstream + 1
                    aUpstream.append(i)
                    aReverse.append(0)
                    aFlag_process[i] = 1
                    pass
            else:
                if pVerter_start == pVertex_start_in and pVerter_end !=pVertex_end_in :
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
        if(check_head_water(aFlowline_in, pVertex_start_in)==1):            
            pass
        else:
            nUpstream, aUpstream, aReverse = find_upstream_flowline(pVertex_start_in, pVertex_end_in)
            if nUpstream > 0:
                global lID
                for j in range(nUpstream):
                    pFlowline = aFlowline_in[ aUpstream[j] ] 
                    if (aReverse[j]==1):             
                        pFlowline.reverse()                    
                        pass
                    else:
                        pass
                    pFlowline.lIndex = lID
                    aFlowline_out.append(pFlowline)
                    lID = lID + 1
                    tag_upstream(  pFlowline.pVertex_start, pFlowline.pVertex_end  )            

                pass
            else:
                pass

    tag_upstream(pVertex_start, pVertex_end)
    
    return aFlowline_out