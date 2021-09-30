import os, sys
from osgeo import ogr, osr, gdal, gdalconst
import numpy as np
from pyflowline.shared.vertex import pyvertex
from pyflowline.shared.edge import pyedge
from pyflowline.shared.flowline import pyflowline

from pyflowline.algorithm.auxiliary.find_vertex_in_list import find_vertex_in_list

def split_flowline(aFlowline_in, aVertex_in):
    aFlowline_out = list()
    nFlowline = len(aFlowline_in)

    for i in range(nFlowline):
        pFlowline = aFlowline_in[i]
        iStream_order = pFlowline.iStream_order
        iFlag_dam = pFlowline.iFlag_dam
        nVertex = pFlowline.nVertex

        iPart = 0
        
        aVertex  = list()

        for j in range(nVertex):
            pVertex = pFlowline.aVertex[j]
            iFlag_exist, lIndex = find_vertex_in_list( aVertex_in,  pVertex)
            if iFlag_exist == 1:
                iPart = iPart + 1
                aVertex.append(j)
                pass
            else:
                pass
        if iPart < 2:
            print('Something is wrong')
            pass

        if iPart ==2:
            aFlowline_out.append(pFlowline)
            pass
        else:

            nLine = iPart-1
            for k in range(nLine):
                t = aVertex[k]
                s = aVertex[k+1]
                aEdge=list()
                for l in range(t,s):
                    pVertex0 = pFlowline.aVertex[l]
                    pVertex1 = pFlowline.aVertex[l+1]
                    pEdge = pyedge(pVertex0, pVertex1)
                    aEdge.append(pEdge)
                    pass

                pFlowline1 = pyflowline(aEdge)
                pFlowline1.iStream_order = iStream_order
                pFlowline1.iFlag_dam = iFlag_dam
                aFlowline_out.append(pFlowline1)
                pass
    
            pass



        pass   

    return aFlowline_out

