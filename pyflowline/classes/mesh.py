from abc import abstractmethod
class mesh(object):
    """Abstract base class for
    """
    
    dLatitude_bot = -90
    dLatitude_top = 90
    dLongitude_left = -180
    dLongitude_right = 180
    sFilename_mesh=''   
    @abstractmethod
    def __init__(self):

        pass
    @abstractmethod
    def calculate_mesh_area(self):

        pass