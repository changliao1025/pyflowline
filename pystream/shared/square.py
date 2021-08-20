from abc import ABCMeta, abstractmethod
import numpy as np
from osgeo import gdal, osr, ogr
from pystream.shared.vertex import pyvertex
from pystream.shared.edge import pyedge
from pystream.shared.cell import pycell

class pysquare(pycell):
    #lIndex=0 
    nFlowline=0
    dLength=0.0
    dArea=0.0
    dx_center=0.0
    dy_center=0.0
    aEdge=None
    aVertex=None
    aFlowline=None
    lCellID  = -1
    aNeighbor=None #the global ID of all neighbors
    nNeighbor=-1

    def __init__(self, aEdge, aVertex):    
        nEdge = len(aEdge)
        if nEdge != 4:
            pass
        else:
            
                
            self.aEdge = aEdge
            self.aVertex = aVertex #the first one and last one are the same
            self.nEdge = 4
            self.nVertex = 4

            dx=0.0
            dy=0.0
            for i in range(4):
                dx = dx + aVertex[i].dx
                dy = dy + aVertex[i].dy
                pass

            self.dx_center = dx/4.0
            self.dy_center = dy/4.0

            pass
        pass
    
    def has_this_edge(self, pEdge_in):
        iFlag_found = 0
        for pEdge in self.aEdge:
            if pEdge.is_overlap(pEdge_in):
                iFlag_found =1 
                break
            else:
                pass       
        
        return iFlag_found

    def which_edge_cross_this_vertex(self, pVertex_in):
        iFlag_found = 0
        pEdge_out = None
        for pEdge in self.aEdge:
            if( pEdge.check_vertex_on_edge(pVertex_in) ==1 ):
                iFlag_found =1
                pEdge_out = pEdge
                break

            else:
                pass

        return iFlag_found, pEdge_out

    def calculate_cell_area(self):
        dLength_edge = self.dLength

        dArea = dLength_edge * dLength_edge 

        self.dArea = dArea
        return dArea
    def calculate_edge_length(self):
        dArea = self.dArea
        dLength_edge = np.sqrt(   dArea   )
        self.dLength = dLength_edge
        return dLength_edge

    def share_edge(self, other):
        iFlag_share = 0
        for pEdge in self.aEdge:
            for pEdge2 in other.aEdge:
                if pEdge.is_overlap(pEdge2) ==1 :
                    iFlag_share = 1 
                    break


        return iFlag_share

