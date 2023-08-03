import enum
from abc import ABCMeta, abstractmethod

# Using enum class create enumerations
class celltype(enum.Enum):
    """_summary_

    Args:
        enum (_type_): For enumeration feature
    """
    hexagon = 1
    square = 2
    latlon = 3
    mpas = 4
    dggrid = 5
    tin = 6

class pycell(metaclass=ABCMeta):
    
    @abstractmethod
    def __init__(self, aEdge):    
        pass

    @abstractmethod
    def calculate_cell_area(self, aEdge):
        pass

