from abc import ABCMeta, abstractmethod
import numpy as np

class pyvertex(object):
    __metaclass__ = ABCMeta  
    dx=0.0
    dy=0.0
    dz=0.0
    dLongitude=0.0
    dLatitude=0.0
    dElevation=0.0
    lIndex=-1 #this index will be used for array
    def __init__(self, aParameter):
        if 'x' in aParameter:            
            self.dx             = float(aParameter['x'])
        
        if 'y' in aParameter:            
            self.dy             = float(aParameter['y'])
        
        if 'z' in aParameter:            
            self.dz             = float(aParameter['z'])
        
        if 'lon' in aParameter:            
            self.dLongitude      = float(aParameter['lon'])
        
        if 'lat' in aParameter:            
            self.dLatitude       = float(aParameter['lat'])

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
        
        x1 = self.dx
        y1 = self.dy
        x2 = other.dx
        y2 = other.dy
        a = (x1-x2) * (x1-x2)
        b = (y1-y2) * (y1-y2)
        c = np.sqrt(a+b)

        dDistance = c

        return dDistance

  


