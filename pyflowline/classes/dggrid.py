import json
from json import JSONEncoder
import numpy as np
from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.edge import pyedge
from pyflowline.classes.cell import pycell
from pyflowline.classes.flowline import pyflowline
from pyflowline.external.pyearth.gis.gdal.gdal_functions import calculate_polygon_area

class DggridClassEncoder(JSONEncoder):
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
        if isinstance(obj, pydggrid):
            return obj.lCellID              
        
            
        return JSONEncoder.default(self, obj)

class pydggrid(pycell):
    """
    The hexagon cell class

    Args:
        pycell (obj): None

    Returns:
        pydggrid: A dggrid object
    """
       
    lCellID  = -1
    nFlowline=0
    nVertex =0 
    nEdge=0  
    dLength=0.0
    dArea=0.0    
    dX_center_meter=0.0
    dY_center_meter=0.0
    dz_center=0.0
    dLongitude_center_degree=0.0
    dLatitude_center_degree=0.0
    dElevation_mean=-9999
    dLength_flowline=0.0
    iFlag_intersected=-1
    iFlag_coast = 0
    lCellID_downstream_burned=-1
    iStream_order_burned=-1
    iStream_segment_burned=-1
    aEdge=None
    aEdgeID=None
    aVertex=None
    aVertexID=None
    pVertex_center = None
    aFlowline=None   
    nNeighbor=-1
    nNeighbor_land=-1
    nNeighbor_ocean=-1
    nNeighbor_land_virtual = -1
    aNeighbor_land_virtual = None
    aNeighbor=None #the global ID of all neighbors
    aNeighbor_land=None #the global ID of all neighbors
    aNeighbor_ocean=None #the global ID of all neighbors
    aNeighbor_distance = None

    def __init__(self, dLon, dLat, aEdge, aVertex):
        """
        Initilize a hexagon cell object

        Args:
            dLon (float): The longitude of center 
            dLat (float): The latitude of center 
            aEdge (list [pyedge]): A list of edges that define the hexagon
            aVertex (list [pyvertex]): A list of vertices the define the hexagon
        """
        nEdge = len(aEdge)
                             
        self.aEdge = aEdge
        self.aVertex = aVertex #the first one and last one are the same
        self.nEdge = nEdge
        self.nVertex = nEdge
        self.dLongitude_center_degree = dLon
        self.dLatitude_center_degree = dLat            
        pVertex = dict()        
        pVertex['dLongitude_degree'] =self.dLongitude_center_degree
        pVertex['dLatitude_degree'] =self.dLatitude_center_degree           
        self.pVertex_center = pyvertex(pVertex)
        self.lCellID_downstream_burned=-1
        self.iStream_order_burned=-1
        self.iStream_segment_burned=-1
        self.dElevation_mean=-9999.0
        pass
      
    
    def has_this_edge(self, pEdge_in):
        """
        Check whether the hexagon contains an edge

        Args:
            pEdge_in (pyedge): The edge to be checked

        Returns:
            int: 1 if found, 0 if not
        """
        iFlag_found = 0
        for pEdge in self.aEdge:
            if pEdge.is_overlap(pEdge_in):
                iFlag_found =1 
                break
            else:
                pass       
        
        return iFlag_found

    def which_edge_cross_this_vertex(self, pVertex_in):
        """
        Find which edge overlap with a vertex

        Args:
            pVertex_in (pyvertex): The vertex to be checked

        Returns:
            tuple [int, pyedge]: 1 if found, with the edge object; 0 if not found
        """
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
        """
        Calculate the area of the hexagon cell

        Returns:
            float: The area in m2
        """

        lons=list()
        lats=list()        
        for i in range(self.nVertex):            
            lons.append( self.aVertex[i].dLongitude_degree )
            lats.append( self.aVertex[i].dLatitude_degree )

        self.dArea = calculate_polygon_area( lons, lats)        
        return self.dArea

    def calculate_edge_length(self):
        """
        Calculate the effective length of the hexagon cell

        Returns:
            float: The effective length
        """
        dArea = self.dArea
        dLength_edge = np.sqrt(  2.0 * dArea / (3.0* np.sqrt(3.0))  )
        self.dLength = dLength_edge
        return dLength_edge
    
    def share_edge(self, other):
        """
        Check whether a hexagon shares an edge with another hexagon

        Args:
            other (pydggrid): The other hexagon

        Returns:
            int: 1 if share, 0 if not
        """
        iFlag_share = 0
        for pEdge in self.aEdge:
            for pEdge2 in other.aEdge:
                if pEdge.is_overlap(pEdge2) ==1 :
                    iFlag_share = 1 
                    break

        return iFlag_share
    
    def tojson(self):
        """
        Convert a hexagon object to a json string

        Returns:
            json str: A json string
        """
        aSkip = ['aEdge', \
                'aFlowline']
                
        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)
        sJson = json.dumps(obj, \
            sort_keys=True, \
            indent = 4, \
            ensure_ascii=True, \
                cls=DggridClassEncoder) 
        return sJson

    

