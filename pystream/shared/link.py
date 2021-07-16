from abc import ABCMeta, abstractmethod

from osgeo import gdal, osr, ogr

from pystream.shared.vertex import pyvertex
from pystream.shared.edge import pyedge
from pystream.shared.cell import pycell
from pystream.shared.flowline import pyflowline
from pystream.shared.mpas import pympas

import json
import numpy as np
from json import JSONEncoder
class LinkClassEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, pyedge):
            return obj.lEdgeID
        if isinstance(obj, pyvertex):
            return obj.lVertexID
        if isinstance(obj, pyflowline):
            return obj.lFlowlineID
        if isinstance(obj, pympas):
            return obj.lCellID    
        return JSONEncoder.default(self, obj)

class pyhexagonlink(object):
    lIndex=0
    
    dLength=0.0
    pHexagon_start=None
    pHexagon_end=None
    pEdge_link=None

    def __init__(self, pHexagon_start_in, pHexagon_end_in, pEdge_link_in):   
        self.pHexagon_start = pHexagon_start_in
        self.pHexagon_end = pHexagon_end_in
        self.pEdge_link = pEdge_link_in
        return

class pympaslink(object):
    lIndex=0
    
    dLength=0.0
    pMpas_start=None
    pMpas_end=None
    pEdge_link=None

    def __init__(self, pMpas_start_in, pMpas_end_in, pEdge_link_in):   
        self.pMpas_start = pMpas_start_in
        self.pMpas_end = pMpas_end_in
        self.pEdge_link = pEdge_link_in
        return

    def export_to_json(self):
        #sJson = json.dumps(self.__dict__, f, ensure_ascii=False, indent=4, cls=CellClassEncoder) 
        sJson = json.dumps(self.__dict__, ensure_ascii=True, indent=4, cls=LinkClassEncoder)
        return sJson
   
    

