
import numpy as np
import importlib
iFlag_cython = importlib.util.find_spec("cython") 

from pyflowline.algorithms.split.find_flowline_vertex import find_flowline_vertex
if iFlag_cython is not None:
    from pyflowline.algorithms.cython.kernel import add_unique_vertex
    from pyflowline.algorithms.cython.kernel import find_vertex_in_list
else:
    from pyflowline.algorithms.auxiliary.find_index_in_list import add_unique_vertex
    from pyflowline.algorithms.auxiliary.find_vertex_in_list import find_vertex_in_list

def find_flowline_confluence(aFlowline_in, pVertex_outlet_in):    
    """_summary_

    Args:
        aFlowline_in (_type_): _description_
        pVertex_outlet_in (_type_): _description_

    Returns:
        List: _description_
    """

    nFlowline = len(aFlowline_in) 
    aVertex=list()
    aIndex_headwater=list()
    aIndex_confluence=list()
    aIndex_middle =list()
    lIndex_outlet = -1

    aVertex = find_flowline_vertex(aFlowline_in)

    vertex_to_index = {vertex: index for index, vertex in enumerate(aVertex)}


    nVertex=len(aVertex)
    aConnectivity  = np.full(  nVertex , 0, dtype=int )
    #iFlag_first=1
    #for i in range(nVertex):        
    #    pVertex = aVertex[i]
    #    dDiatance = pVertex.calculate_distance( pVertex_outlet_in) 
    #    if iFlag_first ==1:
    #        dDiatance_min = dDiatance                
    #        lIndex_outlet = i            
    #        pVertex_outlet_out =  pVertex
    #        iFlag_first=0
    #    else:
    #        if  dDiatance < dDiatance_min:
    #            dDiatance_min = dDiatance
    #            #found it
    #            lIndex_outlet = i       
    #            pVertex_outlet_out =  pVertex
    #            pass    
    #        else:
    #            pass

    distances = [vertex.calculate_distance(pVertex_outlet_in) for vertex in aVertex]
    lIndex_outlet = np.argmin(distances)
    pVertex_outlet_out = aVertex[lIndex_outlet]

    #a outlet can be a confluence, so we won't set it     
    aConnectivity[lIndex_outlet] =0 

    #for i in range(0, nFlowline):      
    #    pFlowline = aFlowline_in[i]
    #    pVertex_start = pFlowline.pVertex_start
    #    pVertex_end = pFlowline.pVertex_end
    #    
    #    iFlag_exist, lIndex =  find_vertex_in_list(aVertex, pVertex_end)  
    #    if lIndex != lIndex_outlet:
    #        #if check_head_water(aFlowline_in, pVertex_start)==1:
    #        if pFlowline.iStream_order == 1:
    #            iFlag_exist, lIndex =  find_vertex_in_list(aVertex, pVertex_start)            
    #            aConnectivity[lIndex] = 1
    #            iFlag_exist, lIndex =  find_vertex_in_list(aVertex, pVertex_end)            
    #            aConnectivity[lIndex] = aConnectivity[lIndex] + 1        
    #        else:
    #            #middle point, or confluence
    #            iFlag_exist, lIndex =  find_vertex_in_list(aVertex, pVertex_start)
    #            if iFlag_exist ==1:
    #                aConnectivity[lIndex] = aConnectivity[lIndex] + 1
#
    #            iFlag_exist, lIndex =  find_vertex_in_list(aVertex, pVertex_end)
    #            if iFlag_exist ==1:
    #                aConnectivity[lIndex] = aConnectivity[lIndex] + 1
    #    else:
    #        #it is outlet, but an outlet can be a confluence as well
    #        iFlag_exist, lIndex =  find_vertex_in_list(aVertex, pVertex_start)
    #        if iFlag_exist ==1:
    #            aConnectivity[lIndex] = aConnectivity[lIndex] + 1
    #        iFlag_exist, lIndex =  find_vertex_in_list(aVertex, pVertex_end)
    #        if iFlag_exist ==1:
    #            aConnectivity[lIndex] = aConnectivity[lIndex] + 1
    #        pass
    for pFlowline in aFlowline_in:
        pVertex_start = pFlowline.pVertex_start
        pVertex_end = pFlowline.pVertex_end
        lIndex_end = vertex_to_index.get(pVertex_end)

        if lIndex_end != lIndex_outlet:
            if pFlowline.iStream_order == 1:
                lIndex_start = vertex_to_index.get(pVertex_start)
                aConnectivity[lIndex_start] = 1
                aConnectivity[lIndex_end] += 1
            else:
                lIndex_start = vertex_to_index.get(pVertex_start)
                if lIndex_start is not None:
                    aConnectivity[lIndex_start] += 1
                if lIndex_end is not None:
                    aConnectivity[lIndex_end] += 1
        else:
            lIndex_start = vertex_to_index.get(pVertex_start)
            if lIndex_start is not None:
                aConnectivity[lIndex_start] += 1
            if lIndex_end is not None:
                aConnectivity[lIndex_end] += 1
    #reset outlet
    
    #for i in range(0, nVertex): 
    #    if i == lIndex_outlet:
    #        if (aConnectivity[i] >1):
    #            aIndex_confluence.append(i)
    #            pass
    #        pass
    #    else:
    #        if (aConnectivity[i] >=3):
    #            aIndex_confluence.append(i)
    #        else:
    #            if (aConnectivity[i] ==1):
    #                aIndex_headwater.append(i)
    #            else:
    #                if (aConnectivity[i] ==2):  
    #                    aIndex_middle.append(i)
    #                else:
    #                    #this is outlet
    #                    pass

    aIndex_headwater = [i for i, x in enumerate(aConnectivity) if x == 1 and i != lIndex_outlet]
    aIndex_middle = [i for i, x in enumerate(aConnectivity) if x == 2 and i != lIndex_outlet]
    aIndex_confluence = [i for i, x in enumerate(aConnectivity) if x >= 3 or (i == lIndex_outlet and x > 1)]
                    
    return aVertex, lIndex_outlet, aIndex_headwater, aIndex_middle, aIndex_confluence, aConnectivity, pVertex_outlet_out
