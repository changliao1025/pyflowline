
import numpy as np
from pyflowline.algorithms.auxiliary.find_index_in_list import find_vertex_in_list

lID = 0
def merge_flowline(aFlowline_in, aVertex_in, \
    pVertex_outlet_in, \
        aIndex_headwater_in,\
            aIndex_middle_in, \
                aIndex_confluence_in ):

    nVertex=len(aVertex_in)
    nFlowline = len(aFlowline_in)
    aFlowline_out=list()               
    aVertex = np.array(aVertex_in)
    aIndex_headwater = np.array(aIndex_headwater_in)
    aIndex_middle = np.array(aIndex_middle_in)
    aIndex_confluence = np.array(aIndex_confluence_in)
    if aIndex_middle.size == 0:
        return aFlowline_in
    
    aVertex_headwater=aVertex[aIndex_headwater]    
    aVertex_middle=aVertex[aIndex_middle]   
    if aIndex_confluence.size > 0:        
        aVertex_confluence=aVertex[aIndex_confluence]    
    
    def merge_flowline_reach(lIndex_in, pVertex_start_in, pVertex_end_in):
        global lID
        pFlowline = aFlowline_in[lIndex_in]
        pVertex_current = pVertex_start_in
        dummy = lIndex_in
        while (find_vertex_in_list(aVertex_middle, pVertex_current)[0] ==1):            
            for j in range(0, nFlowline):      
                pFlowline2 = aFlowline_in[j]                
                pVertex_start = pFlowline2.pVertex_start
                pVertex_end = pFlowline2.pVertex_end
                if pVertex_end == pVertex_current:
                    pFlowline = pFlowline.merge_upstream(pFlowline2)
                    pVertex_current = pVertex_start
                    dummy = j
                    break
                else:
                    pass

        pFlowline.lIndex = lID
        aFlowline_out.append(pFlowline)        
        lID = lID + 1        
        #go to next 
        if find_vertex_in_list(aVertex_headwater, pVertex_current)[0] ==1: 
            return
        else:
            #confluence
            if find_vertex_in_list(aVertex_confluence, pVertex_current)[0] ==1: 
                for k in range(0, nFlowline):                      
                    pFlowline3 = aFlowline_in[k]                
                    pVertex_start = pFlowline3.pVertex_start
                    pVertex_end = pFlowline3.pVertex_end
                    if pVertex_end == pVertex_current:
                        merge_flowline_reach(k, pVertex_start, pVertex_end)
                                
    
    iFlag_first=1
    for i in range(nFlowline):        
        pFlowline = aFlowline_in[i]                
        pVertex_start = pFlowline.pVertex_start
        pVertex_end = pFlowline.pVertex_end
        dDiatance = pVertex_end.calculate_distance( pVertex_outlet_in)
        if iFlag_first ==1:
            dDiatance_min = dDiatance                
            lIndex_outlet = i            
            iFlag_first=0    
        else:            
            if  dDiatance < dDiatance_min:
                dDiatance_min = dDiatance           
                lIndex_outlet = i                
                pass    
            else:
                pass
            
    pFlowline = aFlowline_in[lIndex_outlet]            
    pVertex_start = pFlowline.pVertex_start
    pVertex_end = pFlowline.pVertex_end     
    merge_flowline_reach(lIndex_outlet, pVertex_start, pVertex_end)   
    return aFlowline_out
