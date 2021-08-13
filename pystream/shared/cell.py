import abc
from abc import ABCMeta, abstractmethod
import numpy as np
import json
from osgeo import gdal, osr, ogr
from pystream.shared.vertex import pyvertex
from pystream.shared.edge import pyedge

class pycell(metaclass=ABCMeta):
    
    @abstractmethod
    def __init__(self, aEdge):    
        pass

    @abstractmethod
    def calculate_cell_area(self, aEdge):

        pass

