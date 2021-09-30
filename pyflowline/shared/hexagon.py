from abc import ABCMeta, abstractmethod
import numpy as np
import json
from osgeo import gdal, osr, ogr
from pyflowline.shared.vertex import pyvertex
from pyflowline.shared.edge import pyedge
from pyflowline.shared.cell import pycell
from pyearth.gis.location.calculate_polygon_area import calculate_polygon_area
import numpy as np
import json
from json import JSONEncoder
class NumpyArrayEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return JSONEncoder.default(self, obj)

class pyhexagon(pycell):
    #lIndex=-1  
    
    nFlowline=0
    dLength=0.0
    dArea=0.0
    dx_center=0.0
    dy_center=0.0

    dLon_center=0.0
    dLat_center=0.0
    aEdge=None
    aVertex=None
    aFlowline=None
    lCellID  = -1
    aNeighbor=None #the global ID of all neighbors

    pVertex_center = None

    def __init__(self, aEdge, aVertex):    

        nEdge = len(aEdge)
        if nEdge != 6:
            pass
        else:
            
                
            self.aEdge = aEdge
            self.aVertex = aVertex #the first one and last one are the same
            self.nEdge = 6
            self.nVertex = 6

            dLon=0.0
            dLat=0.0
            for i in range(self.nVertex):
                dLon = dLon + aVertex[i].dLongitude
                dLat = dLat + aVertex[i].dLatitude
                pass

            self.dLon_center = dLon/self.nVertex
            self.dLat_center = dLat/self.nVertex
            pVertex = dict()        
            pVertex['lon'] =self.dLon_center
            pVertex['lat'] =self.dLat_center           
            self.pVertex_center = pyvertex(pVertex)

            self.lCellID_downstream_burned=-1
            self.iStream_order_burned=-1
            self.iStream_segment_burned=-1
            self.dElevation=-9999.0

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
        #dLength_edge = self.dLength

        #dLength_edge = np.sqrt(  2.0 * dArea / (3.0* np.sqrt(3.0))  )
        #dArea = dLength_edge * dLength_edge * (3.0* np.sqrt(3.0)) /2.0

        lons=list()
        lats=list()
        
        for i in range(self.nVertex):
            
            lons.append( self.aVertex[i].dLongitude )
            lats.append( self.aVertex[i].dLatitude )


        self.dArea = calculate_polygon_area(lats, lons)

        
        return self.dArea

    def calculate_edge_length(self):
        dArea = self.dArea
        dLength_edge = np.sqrt(  2.0 * dArea / (3.0* np.sqrt(3.0))  )
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
    
    def export_to_json(self):
        sJson = json.dumps(self.__dict__, ensure_ascii=False, indent=4, cls=NumpyArrayEncoder) 
        return sJson

