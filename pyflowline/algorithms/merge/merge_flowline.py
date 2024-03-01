
import numpy as np
import importlib
iFlag_cython = importlib.util.find_spec("cython") 
if iFlag_cython is not None:
    from pyflowline.algorithms.cython.kernel import find_vertex_in_list
else:
    from pyflowline.algorithms.auxiliary.find_vertex_in_list import find_vertex_in_list

lID = 0

def merge_flowline(aFlowline_in, 
                   aVertex_in, 
                    pVertex_outlet_in, 
                    aIndex_headwater_in,
                    aIndex_middle_in, 
                    aIndex_confluence_in ):
    """Merge the flowline feature to its simplist format

    Args:
        aFlowline_in (_type_): _description_
        aVertex_in (_type_): _description_
        pVertex_outlet_in (_type_): _description_
        aIndex_headwater_in (_type_): _description_
        aIndex_middle_in (_type_): _description_
        aIndex_confluence_in (_type_): _description_

    Returns:
        List: The flowline are strictly ordered from outlet to headwater
    """

    #nVertex=len(aVertex_in)
    nFlowline = len(aFlowline_in)
    aFlowline_out=list()               
    aVertex = np.array(aVertex_in)
    aIndex_headwater = np.array(aIndex_headwater_in)
    aIndex_middle = np.array(aIndex_middle_in)
    aIndex_confluence = np.array(aIndex_confluence_in)
    if aIndex_middle.size == 0:
        return aFlowline_in 

    #convert to set
    aVertex_headwater_set = set(aVertex[aIndex_headwater])  
    aVertex_middle_set = set(aVertex[aIndex_middle])  

    if aIndex_confluence.size > 0:        
        iFlag_confluence = 1  
        aVertex_confluence_set = set(aVertex[aIndex_confluence])  
    else:
        iFlag_confluence = 0
    
    # Build the dictionary
    vertex_to_flowline = {flowline.pVertex_end: flowline for flowline in aFlowline_in}  

    def merge_flowline_reach(lIndex_in, pVertex_start_in):
        global lID
        pFlowline = aFlowline_in[lIndex_in]
        iSegment = pFlowline.iStream_segment
        pVertex_current = pVertex_start_in
        
        #while (find_vertex_in_list(aVertex_middle.tolist(), pVertex_current)[0] ==1):            
        while (pVertex_current in aVertex_middle_set): 
            if pVertex_current in vertex_to_flowline:
                pFlowline2 = vertex_to_flowline[pVertex_current]
                pFlowline = pFlowline.merge_upstream(pFlowline2)
                pVertex_current = pFlowline2.pVertex_start

            #old method    
            #for j in range(0, nFlowline):      
            #    pFlowline2 = aFlowline_in[j]                
            #    pVertex_start = pFlowline2.pVertex_start
            #    pVertex_end = pFlowline2.pVertex_end
            #    if pVertex_end == pVertex_current:
            #        pFlowline = pFlowline.merge_upstream(pFlowline2)
            #        pVertex_current = pVertex_start                    
            #        break                
               
        #go to next 
        #if find_vertex_in_list(aVertex_headwater.tolist(), pVertex_current)[0] ==1: 
        if pVertex_current in aVertex_headwater_set:
            pFlowline.iStream_segment = iSegment      
            pFlowline.iStream_order = 1           
            pFlowline.lFlowlineIndex = lID
            aFlowline_out.append(pFlowline)        
            lID = lID + 1 
            return
        else:
            pFlowline.iStream_segment = iSegment              
            pFlowline.lFlowlineIndex = lID
            aFlowline_out.append(pFlowline)        
            lID = lID + 1 
            #confluence
            #if find_vertex_in_list(aVertex_confluence.tolist(), pVertex_current)[0] ==1: 
            if pVertex_current in aVertex_confluence_set:
                for k in range(0, nFlowline):                      
                    pFlowline3 = aFlowline_in[k]                
                    pVertex_start = pFlowline3.pVertex_start
                    pVertex_end = pFlowline3.pVertex_end
                    if pVertex_end == pVertex_current:
                        merge_flowline_reach(k, pVertex_start)
                               
        
    #find outlet
    dDiatance_min = float('inf')
    for i in range(nFlowline):        
        pFlowline = aFlowline_in[i]                
        pVertex_end = pFlowline.pVertex_end
        dDiatance = pVertex_end.calculate_distance(pVertex_outlet_in)
        if dDiatance < dDiatance_min:
            dDiatance_min = dDiatance
            lIndex_outlet = i
            
    pFlowline = aFlowline_in[lIndex_outlet]            
    pVertex_start = pFlowline.pVertex_start
    pVertex_end = pFlowline.pVertex_end   

    

    #now start from outlet
    if iFlag_confluence == 1:
        #check whether outlet is a confluence
        #if  (find_vertex_in_list(aVertex_confluence.tolist(), pVertex_end)[0] ==1):
        if pVertex_end in aVertex_confluence_set:
            #if pVertex_end in end_vertex_to_flowline:
            #    merge_flowline_reach(end_vertex_to_flowline[pVertex_end], pVertex_start_dummy)      
            for i in range(nFlowline):
                pFlowline = aFlowline_in[i]  
                pVertex_start_dummy = pFlowline.pVertex_start
                pVertex_end_dummy = pFlowline.pVertex_end
                if pVertex_end == pVertex_end_dummy:
                    #this is 
                    merge_flowline_reach(i, pVertex_start_dummy)  
            
        else:
            merge_flowline_reach(lIndex_outlet, pVertex_start)   
    else:
        merge_flowline_reach(lIndex_outlet, pVertex_start)   
    return aFlowline_out
