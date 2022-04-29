from abc import ABCMeta, abstractmethod
import json
from json import JSONEncoder
import numpy as np
from pyflowline.classes.vertex import pyvertex
from pyflowline.algorithms.auxiliary.gdal_functions import  calculate_angle_betwen_vertex, \
    calculate_polygon_area, calculate_distance_to_plane

class EdgeClassEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.float32):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, list):
            pass  
        if isinstance(obj, pyvertex):
            return json.loads(obj.tojson()) #lVertexID
          
            
        return JSONEncoder.default(self, obj)

class pyedge(object):
    __metaclass__ = ABCMeta 
    lEdgeID=-1
    lIndex=-1
    pVertex_start = None
    pVertex_end = None
    dLength=0.0    
    lIndex_upstream=-1
    lIndex_downstream=-1

    def __init__(self, pVertex_start_in, pVertex_end_in):
        if pVertex_start_in == pVertex_end_in:
            print('The two vertices are the same')    
        else:
            self.pVertex_start = pVertex_start_in
            self.pVertex_end = pVertex_end_in
            self.dLength = self.calculate_length()

        return

    def calculate_length(self):
        dLength =0.0
        dLength = self.pVertex_start.calculate_distance( self.pVertex_end)
        self.dLength= dLength
        return dLength

    def check_shared_vertex(self, other):
        '''
        check whether two edges are sharing the same vertex
        '''
        iFlag_shared =-1
        v0 = self.pVertex_start
        v1 = self.pVertex_end
        v2 = other.pVertex_start
        v3 = other.pVertex_end
        if v0 == v2 or v0 == v2 or v1==v2 or v1==v3:
            iFlag_shared =1
        else:
            iFlag_shared=0

        return iFlag_shared
    
    def check_upstream(self, other):
        iFlag_upstream =-1
        v0 = self.pVertex_start
        v1 = self.pVertex_end
        v2 = other.pVertex_start
        v3 = other.pVertex_end
        if v0 == v3:
            iFlag_upstream =1
        else:
            iFlag_upstream=0

        return iFlag_upstream

    def check_downstream(self, other):
        iFlag_downstream =-1
        v0 = self.pVertex_start
        v1 = self.pVertex_end
        v2 = other.pVertex_start
        v3 = other.pVertex_end
        if v1 == v2:
            iFlag_downstream =1
        else:
            iFlag_downstream=0

        return iFlag_downstream

    def reverse(self):
        v0 = self.pVertex_start
        v1 = self.pVertex_end
        self.pVertex_start = v1
        self.pVertex_end = v0
    
    def is_overlap(self, pEdge_in):
        iFlag_overlap = 0
        pVertex_start1 = self.pVertex_start
        pVertex_end1 = self.pVertex_end
        pVertex_start2 = pEdge_in.pVertex_start
        pVertex_end2 = pEdge_in.pVertex_end

        if pVertex_start1 == pVertex_start2 and pVertex_end1 == pVertex_end2:
            iFlag_overlap = 1
        else:
            if  pVertex_start1 == pVertex_end2 and pVertex_end1 == pVertex_start2:
                iFlag_overlap = 1
            else:
                iFlag_overlap = 0

        return iFlag_overlap

    def check_vertex_on_edge(self, pVertex_in):
        iFlag =0 
        dDistance = -1
        dDistance_plane = 9999
        pVertex_start = self.pVertex_start
        pVertex_end = self.pVertex_end
        self.dLength = pVertex_start.calculate_distance(pVertex_end)
        if pVertex_in != pVertex_start and pVertex_in!=pVertex_end:
            d1 = pVertex_start.calculate_distance(pVertex_in)            
            d2 = pVertex_end.calculate_distance(pVertex_in)  
            d3 = d1 + d2 - self.dLength
            angle3deg = calculate_angle_betwen_vertex(\
                 pVertex_start.dLongitude_degree, pVertex_start.dLatitude_degree,\
                 pVertex_in.dLongitude_degree, pVertex_in.dLatitude_degree,\
                 pVertex_end.dLongitude_degree,pVertex_end.dLatitude_degree)

            dDistance_plane = calculate_distance_to_plane(\
                 pVertex_start.dLongitude_degree, pVertex_start.dLatitude_degree,\
                 pVertex_in.dLongitude_degree, pVertex_in.dLatitude_degree,\
                 pVertex_end.dLongitude_degree,pVertex_end.dLatitude_degree)
            lons = [pVertex_start.dLongitude_degree,pVertex_in.dLongitude_degree,pVertex_end.dLongitude_degree]
            lats = [pVertex_start.dLatitude_degree, pVertex_in.dLatitude_degree, pVertex_end.dLatitude_degree]
            dArea = calculate_polygon_area(lons, lats)

            if  angle3deg > 178 and d3 < 1.0: #care
                iFlag = 1
                dDistance = d1
            else:
                iFlag = 0
            
        else:
                iFlag = 0 

        return iFlag, dDistance, dDistance_plane
    
    def __eq__(self, other):                
        iFlag_overlap = self.is_overlap(other)  
        return iFlag_overlap

    def __ne__(self, other):
        return not self.__eq__(other)

    def tojson(self):

        obj = self.__dict__.copy()
        
        sJson = json.dumps(obj, \
            sort_keys=True, \
                indent = 4, \
                    ensure_ascii=True, \
                        cls=EdgeClassEncoder)
        return sJson

