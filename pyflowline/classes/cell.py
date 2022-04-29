import abc
from abc import ABCMeta, abstractmethod

from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.edge import pyedge

class pycell(metaclass=ABCMeta):
    
    @abstractmethod
    def __init__(self, aEdge):    
        pass

    @abstractmethod
    def calculate_cell_area(self, aEdge):
        pass

