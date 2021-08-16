import os, sys
from osgeo import ogr, osr, gdal, gdalconst
import numpy as np
from pystream.shared.vertex import pyvertex
from pystream.shared.edge import pyedge
from pystream.shared.flowline import pyflowline

from pystream.algorithm.auxiliary.find_vertex_in_list import find_vertex_in_list

def split_flowline_to_reach(aFlowline_in):
    aFlowline_out = list()
    nFlowline = len(aFlowline_in)

    for i in range(nFlowline):
        pFlowline = aFlowline_in[i]
        iStream_order = pFlowline.iStream_order
        nVertex = pFlowline.nVertex
        nEdge = pFlowline.nVertex
        iPart = 0
        
        

        for j in range(nEdge):
            pEdge = pFlowline.aEdge[j]
            aEdge = list()
            aEdge.append(pEdge)
            aEdge = pyflowline(aEdge)
            pFlowline1 = pyflowline(aEdge)      
            aFlowline_out.append(pFlowline1)
    
            pass



        pass   

    return aFlowline_out

