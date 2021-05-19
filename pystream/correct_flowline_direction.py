import os, sys
from pystream.shared.flowline import pyflowline
import numpy as np 
from osgeo import ogr, osr, gdal, gdalconst
from shapely.geometry import Point, LineString
from shapely.ops import split
from shapely.wkt import loads

lID=0
aFlag_process=None

def correct_flowline_direction(aFlowline_in, pVertex_outlet):
        
    #we have to go reversely    
    aFlowline_out= list()   

    global lID    
    global aFlag_process

    nFlowline = len(aFlowline_in)
    aFlag_process=np.full(nFlowline, 0, dtype =int)

    for i in range(nFlowline):
        pFlowline = aFlowline_in[i]
        pVertex_start = pFlowline.pVertex_start
        pVertex_end = pFlowline.pVertex_end
        dDiatance = pVertex_end.calculate_distance( pVertex_outlet)
        if  dDiatance < 100.0:
            pFlowline.lIndex = lID
            aFlowline_out.append(pFlowline)
            lID = lID +1
            break
            pass    
        else:
            #print(dDiatance)
            pass
            
        pass     

    def check_head_water(pVertex_start_in):
        iFlag_head_water = -1
        iCount = 0
        for i in range(nFlowline):
            pFlowline = aFlowline_in[i]
            pVerter_start = pFlowline.pVertex_start
            pVerter_end = pFlowline.pVertex_end
            if pVertex_start_in == pVerter_end:
                iCount = iCount + 1
                pass

            if  pVertex_start_in == pVerter_start:
                iCount = iCount + 1
                pass

            pass
        if iCount == 1:
            iFlag_head_water=1
            
        return iFlag_head_water

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
                if aFlag_process[i] !=1:
                    nUpstream = nUpstream + 1
                    aUpstream.append(i)
                    aReverse.append(0)
                    aFlag_process[i] = 1
                    pass
            else:
                if pVerter_start == pVertex_start_in and pVerter_end !=pVertex_end_in :
                    if aFlag_process[i] !=1:
                        nUpstream = nUpstream + 1
                        aUpstream.append(i)
                        aReverse.append(1)
                        aFlag_process[i] = 1
                        pass
                pass

            pass

        print(aUpstream)
        return nUpstream, aUpstream, aReverse
    
    
    def tag_upstream(pVertex_start_in, pVertex_end_in):
        if(check_head_water(pVertex_start_in)==1):            
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