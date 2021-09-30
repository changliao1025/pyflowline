import os, sys
from osgeo import ogr, osr, gdal, gdalconst
import numpy as np

from pyflowline.algorithm.auxiliary.find_vertex_in_list import find_vertex_in_list
lID = 0
def merge_flowline(aFlowline_in, aVertex, \
    pVertex_outlet_in, \
        aIndex_headwater,\
            aIndex_middle, \
                aIndex_confluence ):

    nVertex=len(aVertex)
    nFlowline = len(aFlowline_in)
    aFlowline_out=list()           
    
    aVertex = np.array(aVertex)
    aIndex_headwater = np.array(aIndex_headwater)
    aIndex_middle = np.array(aIndex_middle)
    aIndex_confluence = np.array(aIndex_confluence)

    if aIndex_middle.size == 0:
        return aFlowline_in
    
    aVertex_headwater=aVertex[aIndex_headwater]    
    aVertex_middle=aVertex[aIndex_middle]   

    if aIndex_confluence.size > 0:        
        aVertex_confluence=aVertex[aIndex_confluence]    
    
    def merge_flowline_reach(lIndex, pVertex_start_in, pVertex_end_in):
        global lID
        pFlowline = aFlowline_in[lIndex]
        pVertex_current = pVertex_start_in
        
        while (find_vertex_in_list(aVertex_middle, pVertex_current)[0] ==1):
            
            for j in range(0, nFlowline):      
                pFlowline2 = aFlowline_in[j]                
                pVertex_start = pFlowline2.pVertex_start
                pVertex_end = pFlowline2.pVertex_end
                if pVertex_end == pVertex_current:
                    #this is the upstream, merge them 
                    pFlowline = pFlowline.merge_upstream(pFlowline2)
                    pVertex_current = pVertex_start
                    #print(j)
                    break
                else:
                    pass

        #save 
        pFlowline.lIndex = lID
        aFlowline_out.append(pFlowline)
        
        lID = lID + 1        
        #go to next 
        if find_vertex_in_list(aVertex_headwater, pVertex_current)[0] ==1: 
            pass
        else:
            #it must be confluence
            if find_vertex_in_list(aVertex_confluence, pVertex_current)[0] ==1: 
                for j in range(0, nFlowline):                      
                    pFlowline3 = aFlowline_in[j]                
                    pVertex_start = pFlowline3.pVertex_start
                    pVertex_end = pFlowline3.pVertex_end
                    if pVertex_end == pVertex_current:
                        merge_flowline_reach(j, pVertex_start, pVertex_end)
                        pass
            else:
                print('something is wrong?')
                pass
            
    
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
                #found it
                lIndex_outlet = i                
                pass    
            else:
                #print(dDiatance)
                pass
            
    pFlowline = aFlowline_in[lIndex_outlet]            
    pVertex_start = pFlowline.pVertex_start
    pVertex_end = pFlowline.pVertex_end      

    merge_flowline_reach(lIndex_outlet, pVertex_start, pVertex_end)   



    return aFlowline_out
