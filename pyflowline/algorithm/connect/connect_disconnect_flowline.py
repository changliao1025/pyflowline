
import os, sys
from pyflowline.shared.edge import pyedge
import numpy as np
from osgeo import ogr, osr, gdal, gdalconst


    
def connect_disconnect_flowline(aFlowline_in, aVertex, aThreshold):
    
    aFlowline_out = list()
    nFlowline = len(aFlowline_in)
    nvertex = len(aVertex)
    for j in range(nFlowline):
        pFlowline = aFlowline_in[j]
        pVertex_start = pFlowline.pVertex_start
        pVertex_end = pFlowline.pVertex_end
        iFlag_found = 0
        for i in range(nvertex):
            if iFlag_found ==1:
                break
                
            pVertex=aVertex[i]
            dDistance1 = pVertex.calculate_distance(pVertex_start)
            dDistance2 = pVertex.calculate_distance(pVertex_end)
            
            if dDistance1 < aThreshold[i]:                
                pFlowline.aVertex.insert(0, pVertex)
                pEdge = pyedge(pVertex, pFlowline.aEdge[0].pVertex_start)
                pFlowline.aEdge.insert(0, pEdge)
                #update count
                pFlowline.nVertex = pFlowline.nVertex + 1
                pFlowline.nEdge = pFlowline.nEdge + 1                
                pFlowline.pVertex_start = pVertex
                iFlag_found = 1
                break
                
            else:
                if dDistance2 < aThreshold[i]:                    
                    nvertex2 = pFlowline.nVertex
                    nedge2 = pFlowline.nEdge
                    pFlowline.aVertex.insert(nvertex2, pVertex)
                    pEdge = pyedge(pFlowline.aEdge[nedge2-1].pVertex_end, pVertex )
                    pFlowline.aEdge.insert(nedge2, pEdge)
                    #update count
                    pFlowline.nVertex = pFlowline.nVertex + 1
                    pFlowline.nEdge = pFlowline.nEdge + 1
                    pFlowline.pVertex_end = pVertex
                    iFlag_found = 1
                    break
                    
                else:
                    pass

            pass
    
        aFlowline_out.append(pFlowline)
        pass

    return aFlowline_out

           

            

  