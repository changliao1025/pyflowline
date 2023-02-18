
import json
from json import JSONEncoder
import importlib
import numpy as np

iFlag_cython = importlib.util.find_spec("cython") 
if iFlag_cython is not None:
    from pyflowline.algorithms.cython.kernel import calculate_distance_based_on_lon_lat
else:
    from pyflowline.algorithms.auxiliary.gdal_functions import calculate_distance_based_on_lon_lat

class VertexClassEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()        
        return JSONEncoder.default(self, obj)

class pyvertex(object):
    """
    The vertex class

    Args:
        object (_type_): None

    Returns:
        pyvertex: A vertex object
    """
 
    
  
    lIndex=-1 
    """this index will be used for array - class variable"""
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
        """
        Initilize a vertex object

        Args:           
            aParameter (dict): A dictionary parameters
        """

        if 'x' in aParameter:            
            self.dX_meter             = float(aParameter['x'])
        
        if 'y' in aParameter:            
            self.dY_meter             = float(aParameter['y'])
        
        if 'z' in aParameter:            
            self.dZ_meter             = float(aParameter['z'])
        
        #longitude and latitude are always required     
        try:     
            self.dLongitude_degree      = float(aParameter['dLongitude_degree'])      
            """dLongitude_degree - object variable"""           
            self.dLatitude_degree       = float(aParameter['dLatitude_degree'])
            """dLatitude_degree - object variable"""

            self.dLongitude_radian = np.radians(self.dLongitude_degree)
            self.dLatitude_radian = np.radians(self.dLatitude_degree)
            
        except:
            print('Initialization of vertex failed!')
        
        return

    def toNvector(self): 
        """
        Note: replicated in LatLon_NvectorEllipsoidal

        Returns:
            pynvector: A nvector object
        """
        
        a = self.dLatitude_radian
        b = self.dLongitude_radian
        c = np.sin(a)
        e = np.cos(a)
        d = np.sin(b) 
        f = np.cos(b)
        #// right-handed vector: x -> 0°E,0°N; y -> 90°E,0°N, z -> 90°N
        x = e * f
        y = e * d
        z = c
        point =dict()
        point['x'] = x
        point['y'] = y
        point['z'] = z
        pNvector = pynvector(point)
        return pNvector
    
    def __eq__(self, other):
        """
        Check whether two vertices are equivalent

        Args:
            other (pyvertex): The other vertex

        Returns:
            int: 1 if equivalent, 0 if not
        """
        iFlag = -1
        dThreshold_in = 1.0E-6        
        c = self.calculate_distance(other)
        if( c <= dThreshold_in ): #be careful
            #print(self.dLongitude_degree ,self.dLatitude_degree , other.dLongitude_degree, other.dLatitude_degree)
            iFlag = 1
        else:
            iFlag = 0       

        return iFlag

    def __ne__(self, other):
        """
        Check whether two vertices are equivalent

        Args:
            other (pyvertex): The other vertex

        Returns:
            int: 0 if equivalent, 1 if not
        """
        return not self.__eq__(other)
    
    def calculate_distance(self, other):
        """
        Calculate the distance between two vertices

        Args:
            other (pyvertex): The other vertex

        Returns:
            float: The great circle distance
        """
        dDistance = 0.0
        lon1 = self.dLongitude_degree
        lat1 = self.dLatitude_degree    
        lon2 = other.dLongitude_degree
        lat2 = other.dLatitude_degree
        dDistance = calculate_distance_based_on_lon_lat(lon1, lat1, lon2, lat2)        
        return dDistance
    
    def tojson(self):
        """
        Convert a vecter object to a json string

        Returns:
            json str: A json string
        """
        aSkip = ['dLongitude_radian', \
                'dLatitude_radian']

        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)
            
        #sJson = json.dumps(self.__dict__, \
        sJson = json.dumps(obj, \
                sort_keys=True, \
                indent = 4, \
                ensure_ascii=True, \
                cls=VertexClassEncoder)
        return sJson

  
class pynvector(object):
    """
    The vector class

    Args:
        object (_type_): _description_

    Returns:
        _type_: _description_
    """
    dX=-9999
    dY=-9999
    dZ=-9999

    def __init__(self, aParameter):
        if 'x' in aParameter:            
            self.dX           = float(aParameter['x'])
        
        if 'y' in aParameter:            
            self.dY            = float(aParameter['y'])
        
        if 'z' in aParameter:            
            self.dZ             = float(aParameter['z'])
        self.dLength = self.calculate_length()
        self.dX = self.dX/self.dLength
        self.dY = self.dY /self.dLength
        self.dZ= self.dZ /self.dLength

        return
    def calculate_length(self):
        self.dLength = np.sqrt( self.dX * self.dX + self.dY * self.dY + self.dZ * self.dZ )
        return self.dLength

    def plus(self, other):
        point =dict()
        point['x'] = self.dX + other.dX
        point['y'] = self.dY + other.dY
        point['z'] = self.dZ + other.dZ

        pnv = pynvector( point )

        return pnv

    def toLatLon(self):

        x = self.dX
        y = self.dY
        z = self.dZ

        a = np.arctan2(z, np.sqrt(x*x + y*y))
        b = np.arctan2(y, x)

        point0= dict()   
        point0['dLongitude_degree'] = np.degrees(b)
        point0['dLatitude_degree'] = np.degrees(a)

        pv = pyvertex(point0)
        return pv
        

