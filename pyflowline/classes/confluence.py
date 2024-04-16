import importlib.util
import json
from json import JSONEncoder
import numpy as np
from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.flowline import pyflowline

iFlag_cython = importlib.util.find_spec("cython") 
if iFlag_cython is not None:
    from pyflowline.algorithms.cython.kernel import calculate_angle_betwen_vertex
else:
    from pyearth.gis.geometry.calculate_angle_betwen_vertex import  calculate_angle_betwen_vertex

class ConfluenceClassEncoder(JSONEncoder):
    """Confluence Class Encoder

    Args:
        JSONEncoder (_type_): _description_
    """
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
    """
    The pyconfluence class
    Returns:
        object: A confluence object
    """

    lIndex=-1 
    lConfluenceID=-1    
    pVertex_confluence=None    
    pFlowline_downstream=None
    aFlowline_upstream=list()
    dAngle_upstream=0.0    

    def __init__(self, pVertex_center, aFlowline_upstream_in, pFlowline_downstream_in):
        """
        Initialize a pyconfluence object

        Args:
            pVertex_center (pyvertex): The center vertex
            aFlowline_upstream_in (list [pyflowline]): A list of upstream flowlines
            pFlowline_downstream_in (pyflowline): The downstream flowline
        """
        try:     
            self.pVertex_confluence      = pVertex_center       
            self.aFlowline_upstream      = aFlowline_upstream_in  
            self.pFlowline_downstream    = pFlowline_downstream_in  
            
        except:
            print('Initialization of confluence failed!')
        
        return
    
    def calculate_branching_angle(self):
        """
        Calcualte the confluence branching angle (https://www.pnas.org/doi/10.1073/pnas.1215218109)

        Returns:
            float: The branching angle in degree
        """
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
        """
        Convert a pyconfluence object to json

        Returns:
            json str: A json string
        """
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