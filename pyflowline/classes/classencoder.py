import numpy as np
import json
from json import JSONEncoder

from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.edge import pyedge
from pyflowline.classes.flowline import pyflowline

from pyflowline.classes.hexagon import pyhexagon
from pyflowline.classes.square import pysquare
from pyflowline.classes.latlon import pylatlon
from pyflowline.classes.mpas import pympas
from pyflowline.classes.basin import pybasin

class ClassEncoder(JSONEncoder):
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
        if isinstance(obj, pyedge):
            return obj.lEdgeID        
        if isinstance(obj, pyflowline):
            return obj.lFlowlineID
        if isinstance(obj, pyhexagon):
            return obj.lCellID
        if isinstance(obj, pysquare):
            return obj.lCellID
        if isinstance(obj, pylatlon):
            return obj.lCellID
        if isinstance(obj, pympas):
            return obj.lCellID
        
        if isinstance(obj, pybasin):
            return json.loads(obj.tojson())    
            
        return JSONEncoder.default(self, obj)