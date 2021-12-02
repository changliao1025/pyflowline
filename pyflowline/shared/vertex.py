from abc import ABCMeta, abstractmethod
import numpy as np
import json
from json import JSONEncoder

from pyearth.gis.location.calculate_distance_based_on_lon_lat import calculate_distance_based_on_lon_lat

class VertexClassEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        
        return JSONEncoder.default(self, obj)

class pyvertex(object):
    __metaclass__ = ABCMeta  
    lVertexID=-1
    dx=-9999
    dy=-9999
    dz=-9999
    dLongitude_degree=0.0
    dLatitude_degree=0.0
    dElevation=0.0
    lIndex=-1 #this index will be used for array
    def __init__(self, aParameter):
        if 'x' in aParameter:            
            self.dx             = float(aParameter['x'])
        
        if 'y' in aParameter:            
            self.dy             = float(aParameter['y'])
        
        if 'z' in aParameter:            
            self.dz             = float(aParameter['z'])
        
        #if 'lon' in aParameter:            
        self.dLongitude_degree      = float(aParameter['lon'])
        
        #if 'lat' in aParameter:            
        self.dLatitude_degree       = float(aParameter['lat'])

        
        return
    
    def __eq__(self, other):
        iFlag = -1
        
        c = self.calculate_distance(other)
        if( c < 0.0001 ): #be careful
            iFlag = 1
        else:
            iFlag = 0       

        return iFlag

    def __ne__(self, other):
        return not self.__eq__(other)
    
    def calculate_distance(self, other):
        
        x1 = self.dx
        y1 = self.dy
        z1 = self.dz
        lon1 = self.dLongitude_degree
        lat1 = self.dLatitude_degree
        
        x2 = other.dx
        y2 = other.dy
        z2 = other.dz
        lon2 = other.dLongitude_degree
        lat2 = other.dLatitude_degree

        #if x1!=-9999 and x2!=-9999:
        #    a = (x1-x2) * (x1-x2)
        #    b = (y1-y2) * (y1-y2)
        #    c = np.sqrt(a+b)
        #else:
        #    #use latitude longitude
        #    pass

        c= calculate_distance_based_on_lon_lat(lon1,  lat1, lon2, lat2)

        dDistance = c

        return dDistance
    
    def tojson(self):
        sJson = json.dumps(self.__dict__, \
            sort_keys=True, \
                indent = 4, \
                    ensure_ascii=True, \
                        cls=VertexClassEncoder)
        return sJson

  


