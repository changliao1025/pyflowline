import sys, os
import numpy as np
from osgeo import gdal, osr, ogr, gdalconst


def remove_flowline_loop(aFlowline_in):
    aFlowline_out=list()

    nFlowline = len(aFlowline_in)
        
    
    def find_paralle_stream(i, pVertex_start_in):        
        ndownstream=0
        aDownstream=list()
        for j in range(nFlowline):
            pFlowline = aFlowline_in[j]
            pVertex_start = pFlowline.pVertex_start
            pVertex_end = pFlowline.pVertex_end
            if pVertex_start == pVertex_start_in and i !=j :
                ndownstream= ndownstream+1
                aDownstream.append(j)
                pass
                
        return ndownstream, aDownstream
    lID=0
    aFlag = np.full(nFlowline, 0, dtype=int)
    for i in range(nFlowline):
        pFlowline = aFlowline_in[i]      
        pVertex_start = pFlowline.pVertex_start
        pVertex_end = pFlowline.pVertex_end
        
        ndownstream , aDownstream = find_paralle_stream(i, pVertex_start)
        if ndownstream == 0:
            if aFlag[i] !=1:
                pFlowline.lIndex = lID
                aFlowline_out.append(pFlowline)
                lID = lID + 1
                aFlag[i]=1
            pass
        else:                     
            #more than one, so we only take the current one
            if(ndownstream>0):                
                if aFlag[i] !=1:
                    pFlowline.lIndex = lID
                    aFlowline_out.append(pFlowline)
                    lID = lID + 1
                    aFlag[i]=1
                    pass

                #set all to treated
                for k in range(ndownstream):
                    aFlag[  aDownstream[k]] = 1
            pass
        
        pass 

    return aFlowline_out