import os, sys
from osgeo import ogr, osr, gdal, gdalconst
import numpy as np

from pystream.add_unique_vertex import add_unique_vertex

def find_flowline_vertex(aFlowline_in):
   
    
    nFlowline = len(aFlowline_in)
    #build dictionary
    aVertex=list()
    for i in range(0, nFlowline):      
        pFlowline = aFlowline_in[i]
        pVertex_start = pFlowline.pVertex_start
        pVertex_end = pFlowline.pVertex_end
        aVertex, dummy = add_unique_vertex(aVertex, pVertex_start)
        aVertex, dummy = add_unique_vertex(aVertex, pVertex_end)
        pass

    return aVertex
