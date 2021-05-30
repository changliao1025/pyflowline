import os
import json
from pystream.shared.edge import pyedge
from osgeo import ogr, osr, gdal, gdalconst
import numpy as np

from shapely.geometry import Point, LineString, MultiLineString
from shapely.wkt import loads

from pystream.shared.vertex import pyvertex
from pystream.shared.flowline import pyflowline
from pystream.shared.hexagon import pyhexagon

def convert_coordinates_to_hexagon(aCoordinates):
    npoint = len(aCoordinates)
    
    aVertex=list()
    for i in range(npoint):
        x = aCoordinates[i][0]
        y = aCoordinates[i][1]
        dummy = dict()
        dummy['x'] =x
        dummy['y'] =y
        pVertex = pyvertex(dummy)
        aVertex.append(pVertex)
        
    aEdge=list()
    for j in range(npoint-1):
        pEdge = pyedge( aVertex[j], aVertex[j+1] )
        aEdge.append(pEdge)
    
    pHexagon = pyhexagon( aEdge, aVertex)
    
    return pHexagon