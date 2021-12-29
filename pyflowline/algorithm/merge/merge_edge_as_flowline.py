
import numpy as np
from pyflowline.algorithm.auxiliary.find_index_in_list import find_edge_in_list
lID = 0
def merge_edge_as_flowline(aEdge_in, pVertex_outlet_in):
    iFlag_first=1
    nEdge = len(aEdge_in)

    def merge_edge_reach(lIndex_in, pVertex_start_in, pVertex_end_in):
        global lID
        pEdge = aEdge_in[lIndex_in]
        pVertex_current = pVertex_start_in
        
        
            
        for j in range(0, nEdge):      
            pEdge2 = aEdge_in[j]                
            pVertex_start = pEdge2.pVertex_start
            pVertex_end = pEdge2.pVertex_end
            if pVertex_end == pVertex_current:
                #this is the upstream, merge them 
                pEdge = pEdge.merge_upstream(pEdge2)
                pVertex_current = pVertex_start
                #print(j)
                break
            else:
                pass

        #save 
        pEdge.lIndex = lID
        aFlowline_out.append(pEdge)
        
        lID = lID + 1        
        #go to next 
        if find_vertex_in_list(aVertex_headwater, pVertex_current)[0] ==1: 
            pass
        else:
            #it must be confluence
            if find_vertex_in_list(aVertex_confluence, pVertex_current)[0] ==1: 
                for j in range(0, nFlowline):                      
                    pFlowline3 = aEdge_in[j]                
                    pVertex_start = pFlowline3.pVertex_start
                    pVertex_end = pFlowline3.pVertex_end
                    if pVertex_end == pVertex_current:
                        merge_flowline_reach(j, pVertex_start, pVertex_end)
                        pass
            else:
                print('something is wrong?')
                pass

    for i in range(nEdge):        
        pEdge = aEdge_in[i]                
        pVertex_start = pEdge.pVertex_start
        pVertex_end = pEdge.pVertex_end
        dDiatance = pVertex_end.calculate_distance( pVertex_outlet_in)
        

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
                #print(dDiatance)
                pass
            
    pEdge = aEdge_in[lIndex_outlet]            
    pVertex_start = pEdge.pVertex_start
    pVertex_end = pEdge.pVertex_end      

    merge_edge_reach(lIndex_outlet, pVertex_start, pVertex_end) 
    pass