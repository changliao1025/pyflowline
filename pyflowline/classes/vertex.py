from abc import ABCMeta
import numpy as np
import json
from json import JSONEncoder
from pyflowline.algorithms.auxiliary.gdal_functions import calculate_distance_based_on_lon_lat

class VertexClassEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()        
        return JSONEncoder.default(self, obj)

class pyvertex(object):
    __metaclass__ = ABCMeta  
    lIndex=-1 #this index will be used for array
    lVertexID=-1
    lFlowlineID = -1  #we use this id only for intersect
    dX_meter=-9999
    dY_meter=-9999
    dZ_meter=-9999
    dLongitude_degree=0.0
    dLatitude_degree=0.0
    dLongitude_radian=0.0
    dLatitude_radian=0.0
    dElevation=0.0    
    def __init__(self, aParameter):
        if 'x' in aParameter:            
            self.dX_meter             = float(aParameter['x'])
        
        if 'y' in aParameter:            
            self.dY_meter             = float(aParameter['y'])
        
        if 'z' in aParameter:            
            self.dZ_meter             = float(aParameter['z'])
        
        #longitude and latitude are always required     
        try:     
            self.dLongitude_degree      = float(aParameter['dLongitude_degree'])                 
            self.dLatitude_degree       = float(aParameter['dLatitude_degree'])
        except:
            print('Initialization of vertex failed!')
        
        return
    
    def __eq__(self, other):
        iFlag = -1
        
        c = self.calculate_distance(other)
        if( c < 1.0E-6 ): #be careful
            iFlag = 1
        else:
            iFlag = 0       

        return iFlag

    def __ne__(self, other):
        return not self.__eq__(other)
    
    def calculate_distance(self, other):                
        lon1 = self.dLongitude_degree
        lat1 = self.dLatitude_degree    
        lon2 = other.dLongitude_degree
        lat2 = other.dLatitude_degree
        dDistance = calculate_distance_based_on_lon_lat(lon1,  lat1, lon2, lat2)        
        return dDistance
    
    def tojson(self):
        sJson = json.dumps(self.__dict__, \
                sort_keys=True, \
                indent = 4, \
                ensure_ascii=True, \
                cls=VertexClassEncoder)
        return sJson

  


