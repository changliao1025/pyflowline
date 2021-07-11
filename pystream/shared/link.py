from abc import ABCMeta, abstractmethod

from osgeo import gdal, osr, ogr

from pystream.shared.vertex import pyvertex
from pystream.shared.edge import pyedge
from pystream.shared.cell import pycell


class pyhexagonlink(object):
    lIndex=0
    
    dLength=0.0
    pHexagon_start=None
    pHexagon_end=None
    pFlowline_link=None

    def __init__(self, pHexagon_start_in, pHexagon_end_in, pFlowline_link_in):   
        self.pHexagon_start = pHexagon_start_in
        self.pHexagon_end = pHexagon_end_in
        self.pFlowline_link = pFlowline_link_in
        return
    
   
    

