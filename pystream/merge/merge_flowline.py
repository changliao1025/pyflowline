import os, sys
from osgeo import ogr, osr, gdal, gdalconst
import numpy as np

from pystream.find_vertex_in_list import find_vertex_in_list
lID = 0
def merge_flowline(aFlowline_in, aVertex, \
    pVertex_outlet_in, \
        aIndex_headwater,\
            aIndex_middle, \
                aIndex_confluence ):
    nVertex=len(aVertex)
    nFlowline = len(aFlowline_in)
    aFlowline_out=list()           
    
    #aVertex_headwater=aVertex[aIndex_headwater]
    aVertex_headwater=aVertex[ aIndex_headwater[0]:aIndex_headwater[len(aIndex_headwater)-1] ]
    #aVertex_middle=aVertex[aIndex_middle]
    aVertex_middle=aVertex[ aIndex_middle[0]:aIndex_middle[len(aIndex_middle)-1] ]
    #aVertex_confluence=aVertex[aIndex_confluence]
    aVertex_confluence=aVertex[ aIndex_confluence[0]:aIndex_confluence[len(aIndex_confluence)-1] ]
    
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
                    print(j)
                    break
                else:
                    pass

        #save 
        pFlowline.lIndex = lID
        aFlowline_out.append(pFlowline)
        
        lID = lID + 1        
        #go to next 
        if find_vertex_in_list(aVertex_headwater, pVertex_current)[0] !=1: 
            pass
        else:
            #it must be confluence
            for j in range(0, nFlowline):                      
                pFlowline3 = aFlowline_in[j]                
                pVertex_start = pFlowline3.pVertex_start
                pVertex_end = pFlowline3.pVertex_end
                if pVertex_end == pVertex_current:
                    merge_flowline_reach(j, pVertex_start, pVertex_end)
                    pass
            
    
    
    for i in range(nFlowline):        
        pFlowline = aFlowline_in[i]                
        pVertex_start = pFlowline.pVertex_start
        pVertex_end = pFlowline.pVertex_end
        dDiatance = pVertex_end.calculate_distance( pVertex_outlet_in)
        if  dDiatance < 100.0:                
            lIndex_outlet = i
            break
            pass
            

    merge_flowline_reach(lIndex_outlet, pVertex_start, pVertex_end)   



    return 
