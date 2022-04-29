from abc import ABCMeta
import numpy as np
import json
from json import JSONEncoder
from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.flowline import pyflowline
from pyflowline.algorithms.auxiliary.gdal_functions import calculate_distance_based_on_lon_lat
from pyflowline.algorithms.auxiliary.gdal_functions import  calculate_angle_betwen_vertex

class ConfluenceClassEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()        
        if isinstance(obj, np.float32):
            return float(obj)        
        if isinstance(obj, list):
            pass  
        if isinstance(obj, pyvertex):
            return json.loads(obj.tojson())       
        if isinstance(obj, pyflowline):
            return obj.lFlowlineID
        
        return JSONEncoder.default(self, obj)

class pyconfluence():
    __metaclass__ = ABCMeta  
    lIndex=-1 
    lConfluenceID=-1    
    pVertex_confluence=None    
    pFlowline_downstream=None
    aFlowline_upstream=list()
    dAngle_upstream=0.0    
    def __init__(self, pVertex_center, aFlowline_upstream_in, pFlowline_downstream_in):        
        try:     
            self.pVertex_confluence      = pVertex_center       
            self.aFlowline_upstream      = aFlowline_upstream_in  
            self.pFlowline_downstream      = pFlowline_downstream_in  
            
        except:
            print('Initialization of confluence failed!')
        
        return
    
    def calculate_branching_angle(self):        
        #normally there are 2 edges meet at confluence        
        if len(self.aFlowline_upstream)==2:
            pFlowline1 = self.aFlowline_upstream[0]
            pFlowline2 = self.aFlowline_upstream[1]
            nedge1 = pFlowline1.nEdge
            nedge2 = pFlowline2.nEdge                
            x1 = pFlowline1.aEdge[nedge1-1].pVertex_start.dLongitude_degree
            y1 = pFlowline1.aEdge[nedge1-1].pVertex_start.dLatitude_degree
            x2 = self.pVertex_confluence.dLongitude_degree
            y2 = self.pVertex_confluence.dLatitude_degree
            x3 = pFlowline2.aEdge[nedge2-1].pVertex_start.dLongitude_degree
            y3 = pFlowline2.aEdge[nedge2-1].pVertex_start.dLatitude_degree
            self.dAngle_upstream = calculate_angle_betwen_vertex(x1, y1, x2, y2, x3, y3)
        else:
            print('multiple upstream')
            print(len(self.aFlowline_upstream))
            self.dAngle_upstream=0.0

        return self.dAngle_upstream    
    
    def tojson(self):
        aSkip = ['aFlowline_upstream']
        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)
        sJson = json.dumps(obj,  \
                sort_keys=True, \
                indent = 4, \
                ensure_ascii=True, \
                cls=ConfluenceClassEncoder)
        return sJson