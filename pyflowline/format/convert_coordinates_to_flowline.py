import os
import json
from pyflowline.classes.edge import pyedge
from osgeo import ogr, osr, gdal, gdalconst
import numpy as np
from pyflowline.algorithm.auxiliary.reproject_coordinates import reproject_coordinates





from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.flowline import pyflowline

def convert_gcs_coordinates_to_flowline(aCoordinates):
    npoint = len(aCoordinates)
    
    aVertex=list()
    for i in range(npoint):
        x = aCoordinates[i][0]
        y = aCoordinates[i][1]
        dummy = dict()
        dummy['lon'] =x
        dummy['lat'] =y
        pVertex = pyvertex(dummy)
        aVertex.append(pVertex)
        
    aEdge=list()
    for j in range(npoint-1):
        pEdge = pyedge( aVertex[j], aVertex[j+1] )
        aEdge.append(pEdge)
    
    pLine = pyflowline( aEdge)
    
    return pLine

def convert_pcs_coordinates_to_flowline(aCoordinates, pSpatialRef_shapefile_in):
    npoint = len(aCoordinates)
    
    aVertex=list()
    for i in range(npoint):
        x = aCoordinates[i][0]
        y = aCoordinates[i][1]
        dummy = dict()
        dummy['x'] =x
        dummy['y'] =y

        lon, lat = reproject_coordinates(x, y, pSpatialRef_shapefile_in)
        dummy['lon'] = lon
        dummy['lat'] = lat
        pVertex = pyvertex(dummy)
        aVertex.append(pVertex)
        
    aEdge=list()
    for j in range(npoint-1):
        pEdge = pyedge( aVertex[j], aVertex[j+1] )
        aEdge.append(pEdge)
    
    pLine = pyflowline( aEdge)
    
    return pLine