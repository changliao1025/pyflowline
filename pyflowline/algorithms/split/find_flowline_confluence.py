
import numpy as np

from pyflowline.algorithms.auxiliary.check_head_water import check_head_water


import importlib
iFlag_cython = importlib.util.find_spec("cython") 
if iFlag_cython is not None:
    from pyflowline.algorithms.cython.kernel import add_unique_vertex
    from pyflowline.algorithms.cython.kernel import find_vertex_in_list
else:
    from pyflowline.algorithms.auxiliary.find_index_in_list import add_unique_vertex
    from pyflowline.algorithms.auxiliary.find_vertex_in_list import find_vertex_in_list

def find_flowline_confluence(aFlowline_in, pVertex_outlet_in):    
    nFlowline = len(aFlowline_in) 
    aVertex=list()
    aIndex_headwater=list()
    aIndex_confluence=list()
    aIndex_middle =list()
    lIndex_outlet = -1
    for i in range(0, nFlowline):      
        pFlowline = aFlowline_in[i]
        pVertex_start = pFlowline.pVertex_start
        pVertex_end = pFlowline.pVertex_end
        aVertex, dummy = add_unique_vertex(aVertex, pVertex_start)
        aVertex, dummy = add_unique_vertex(aVertex, pVertex_end)
        pass

    nVertex=len(aVertex)
    aConnectivity  = np.full(  nVertex , 0, dtype=int )
    iFlag_first=1
    for i in range(nVertex):        
        pVertex = aVertex[i]
        dDiatance = pVertex.calculate_distance( pVertex_outlet_in) 
        if iFlag_first ==1:
            dDiatance_min = dDiatance                
            lIndex_outlet = i            
            pVertex_outlet_out =  pVertex
            iFlag_first=0
        else:
            if  dDiatance < dDiatance_min:
                dDiatance_min = dDiatance
                #found it
                lIndex_outlet = i       
                pVertex_outlet_out =  pVertex
                pass    
            else:
                pass

    #a outlet can be a confluence, so we won't set it     
    #aConnectivity[lIndex_outlet] =0 
    for i in range(0, nFlowline):      
        pFlowline = aFlowline_in[i]
        pVertex_start = pFlowline.pVertex_start
        pVertex_end = pFlowline.pVertex_end
        

        iFlag_exist, lIndex =  find_vertex_in_list(aVertex, pVertex_end)  
        if lIndex != lIndex_outlet:

            if check_head_water(aFlowline_in, pVertex_start)==1:
                iFlag_exist, lIndex =  find_vertex_in_list(aVertex, pVertex_start)            
                aConnectivity[lIndex] = 1
                iFlag_exist, lIndex =  find_vertex_in_list(aVertex, pVertex_end)            
                aConnectivity[lIndex] = aConnectivity[lIndex] + 1        
            else:
                #middle point, or confluence
                iFlag_exist, lIndex =  find_vertex_in_list(aVertex, pVertex_start)
                if iFlag_exist ==1:
                    aConnectivity[lIndex] = aConnectivity[lIndex] + 1

                iFlag_exist, lIndex =  find_vertex_in_list(aVertex, pVertex_end)
                if iFlag_exist ==1:
                    aConnectivity[lIndex] = aConnectivity[lIndex] + 1
        else:
            #it is outlet, but an outlet can be a confluence as well
            iFlag_exist, lIndex =  find_vertex_in_list(aVertex, pVertex_start)
            if iFlag_exist ==1:
                aConnectivity[lIndex] = aConnectivity[lIndex] + 1
            iFlag_exist, lIndex =  find_vertex_in_list(aVertex, pVertex_end)
            if iFlag_exist ==1:
                aConnectivity[lIndex] = aConnectivity[lIndex] + 1
            pass

    #reset outlet
    #aConnectivity[lIndex_outlet] = 0
    for i in range(0, nVertex): 
        if i == lIndex_outlet:
            if (aConnectivity[i] >1):
                aIndex_confluence.append(i)
                pass
            pass
        else:
            if (aConnectivity[i] >=3):
                aIndex_confluence.append(i)
            else:
                if (aConnectivity[i] ==1):
                    aIndex_headwater.append(i)
                else:
                    if (aConnectivity[i] ==2):  
                        aIndex_middle.append(i)
                    else:
                        #this is outlet
                        pass


    return aVertex, lIndex_outlet, aIndex_headwater, aIndex_middle, aIndex_confluence, aConnectivity, pVertex_outlet_out
