
import numpy as np
from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.edge import pyedge
from pyflowline.classes.cell import pycell

class pytin(pycell):    
    nFlowline=0
    nVertex =0 
    nEdge=0
    dLength=0.0
    dArea=0.0
    dX_center_meter=0.0
    dY_center_meter=0.0
    aEdge=None
    aVertex=None
    aFlowline=None
    lCellID  = -1
    aNeighbor=None #the global ID of all neighbors
    nNeighbor=-1
    def __init__(self, aEdge,aVertex, dLon, dLat):       
        nEdge = len(aEdge)
        if nEdge !=3:
            pass
        else:           
            self.aEdge = aEdge
            self.aVertex = aVertex #the first one and last one are the same
            self.nEdge = len(aEdge)
            self.nVertex = len(aVertex) - 1
            self.dLongitude_center = dLon
            self.dLatitude_center = dLat
            pVertex = dict()        
            pVertex['lon'] =self.dLongitude_center
            pVertex['lat'] =self.dLatitude_center           
            self.pVertex_center = pyvertex(pVertex)
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
            iFlag, dummy ,diff = pEdge.check_vertex_on_edge(pVertex_in)
            if( iFlag ==1 ):
                iFlag_found =1
                pEdge_out = pEdge
                break
            else:
                pass

        return iFlag_found, pEdge_out
    
    def calculate_cell_area(self):           
        self.dArea = 0.0
        return self.dArea

    def calculate_edge_length(self):        
        self.dLength_edge =0.0
        return self.dLength_edge
    
    def share_edge(self, other):
        iFlag_share = 0
        for pEdge in self.aEdge:
            for pEdge2 in other.aEdge:
                if pEdge.is_overlap(pEdge2) ==1 :
                    iFlag_share = 1 
                    break

        return iFlag_share
