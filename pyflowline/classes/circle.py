import os
import json
from json import JSONEncoder
import importlib.util
import numpy as np
from osgeo import ogr, gdal, osr
from geographiclib.geodesic import Geodesic

from gcsbuffer.classes.vertex import pyvertex
from gcsbuffer.formats.save_vertex_as_polygon import save_vertex_as_polygon
from gcsbuffer.algorithms.split_by_length import split_edge_by_length
from gcsbuffer.algorithms.find_minimal_enclosing_polygon import find_minimal_enclosing_polygon


iFlag_cython = importlib.util.find_spec("cython")
if iFlag_cython is not None:
    from pyflowline.algorithms.cython.kernel import calculate_angle_betwen_vertex
    from pyflowline.algorithms.cython.kernel import calculate_distance_to_plane
else:
    from pyearth.gis.geometry.calculate_angle_betwen_vertex import  calculate_angle_betwen_vertex
    from pyearth.gis.geometry.calculate_distance_to_plane import calculate_distance_to_plane

class CircleClassEncoder(JSONEncoder):
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
            return json.loads(obj.tojson())


        return JSONEncoder.default(self, obj)

class pycircle(object):
    """The pyedge class

    Args:
        object (object): None

    Returns:
        pyedge: A edge object
    """

    pVertex_center = None
    aVertex_circle = None
    
    def __init__(self, pVertex_center_in, aVertex_circle_in):
        """
        Initilize a pyedge object

        Args:
            pVertex_start_in (pyvertex): The starting vertex
            pVertex_end_in (pyvertex): The ending vertex
        """
        self.pVertex_center = pVertex_center_in
        self.aVertex_circle = aVertex_circle_in


        return