import abc
from abc import ABCMeta, abstractmethod

from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.edge import pyedge

import enum
# Using enum class create enumerations
class celltype(enum.Enum):
   hexagon = 1
   square = 2
   latlon = 3
   mpas = 4
   tin = 5

class pycell(metaclass=ABCMeta):
    
    @abstractmethod
    def __init__(self, aEdge):    
        pass

    @abstractmethod
    def calculate_cell_area(self, aEdge):
        pass

