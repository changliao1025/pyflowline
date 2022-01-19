import os, sys
from telnetlib import IP
import numpy as np
from osgeo import ogr, osr, gdal, gdalconst
from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.edge import pyedge
from pyflowline.classes.flowline import pyflowline

from pyflowline.algorithms.auxiliary.find_index_in_list import find_vertex_in_list, find_vertex_on_edge


def split_flowline(aFlowline_in, aVertex_in, iFlag_intersect = None):
    aFlowline_out = list()
    nFlowline = len(aFlowline_in)

    for i in range(nFlowline):
        pFlowline = aFlowline_in[i]
        iStream_order = pFlowline.iStream_order
        iStream_segment = pFlowline.iStream_segment
        
        iFlag_dam = pFlowline.iFlag_dam
        nVertex = pFlowline.nVertex
        nEdge= pFlowline.nEdge

        iPart = 0
        
        aVertex  = list()
        aVertex_all = list()
        for j in range(nEdge):
            pEdge=pFlowline.aEdge[j]
            pVertex = pEdge.pVertex_start
            aVertex_all.append(pVertex)
            iFlag_exist, lIndex = find_vertex_in_list( aVertex_in,  pVertex)
            if iFlag_exist == 1:                
                iPart = iPart + 1
                aVertex.append(pVertex)
                
            if iFlag_intersect is not None:
                iFlag_exist, npoint, aIndex = find_vertex_on_edge( aVertex_in,  pEdge)
                #they must be ordered
                if iFlag_exist==1:
                    for m in range(npoint):
                        pVertex_dummy = aVertex_in[aIndex[m]]
                    
                        iPart = iPart + 1
                        aVertex.append(pVertex_dummy)
                        aVertex_all.append(pVertex_dummy)
                    pass

        #the last vertex
        pVertex = pFlowline.pVertex_end
        aVertex_all.append(pVertex)
        iFlag_exist, lIndex = find_vertex_in_list( aVertex_in,  pVertex)
        if iFlag_exist == 1:
            iPart = iPart + 1
            aVertex.append(pVertex)
            pass


        if iPart == 0 :
            print('Something is wrong')
            pass
        else:
            if iPart ==1:
                print('Something is wrong')
                pass
            else:
                if iPart ==2:
                    aFlowline_out.append(pFlowline)
                    pass
                else: #greater than 2

                    nLine = iPart-1
                    #rebuild index
                    aVertex_index =list()
                    for m in range(iPart):
                        pVertex= aVertex[m]
                        iFlag_exist, lIndex = find_vertex_in_list( aVertex_all,  pVertex)
                        if iFlag_exist ==1:
                            aVertex_index.append(lIndex)
                    
                    #find duplicate


                    for k in range(nLine):
                        t = aVertex_index[k]
                        s = aVertex_index[k+1]
                        if s!=t:
                            aEdge=list()
                            for l in range(t,s):
                                pVertex0 = aVertex_all[l]  #pFlowline.  aVertex[l]
                                pVertex1 = aVertex_all[l+1]  #pFlowline.    aVertex[l+1]
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

