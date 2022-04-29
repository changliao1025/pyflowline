
import numpy as np 
from pyflowline.algorithms.auxiliary.check_head_water import check_head_water
def define_stream_order(aFlowline_in):
    nFlowline = len(aFlowline_in)
    aFlowline_out = list()
    
    if nFlowline == 0 :
        print ('data incomplete')
    else:       
        aStream_order = np.full(nFlowline, 0, dtype=int)    
        for i in range(nFlowline):
            pFlowline = aFlowline_in[i]
            pVertex_start=pFlowline.pVertex_start
            if check_head_water(aFlowline_in, pVertex_start)==1:     
                aStream_order[i] = 1
                pass       
     
        while aStream_order[0] == 0:
            for  i in range(nFlowline):
                if aStream_order[i] !=0:
                    continue

                pFlowline = aFlowline_in[i]
                iFlag_upstream_done = 1
                pVertex_start=pFlowline.pVertex_start
                pVertex_end=pFlowline.pVertex_end
                aStrord=list()
                for  j in range(nFlowline):
                    pFlowline2 = aFlowline_in[j]
                    pVertex_start2=pFlowline2.pVertex_start
                    pVertex_end2=pFlowline2.pVertex_end
                    if pVertex_start == pVertex_end2:
                        if aStream_order[j] == 0:
                            iFlag_upstream_done=0
                            break
                        else:
                            aStrord.append( aStream_order[j]  )
                
                if iFlag_upstream_done==1:
                    #get unique value
                    dummy = np.array(aStrord)
                    dummy1 = np.unique(dummy)

                    if len(dummy1) == 1: #all upstreams have the same order
                        aStream_order[i] = aStrord[0] + 1
                    else:
                        aStream_order[i] = np.max(dummy)    

        for i in range(nFlowline):
            pFlowline = aFlowline_in[i]      
            pFlowline.iStream_order =    aStream_order[i]   
            aFlowline_out.append(pFlowline)

    return aFlowline_out, aStream_order