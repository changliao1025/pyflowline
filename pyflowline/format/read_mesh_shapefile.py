import os
import json
from pyflowline.classes.edge import pyedge
from osgeo import ogr, osr, gdal, gdalconst
import numpy as np

from shapely.geometry import Point, LineString, MultiLineString
from shapely.wkt import loads

from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.flowline import pyflowline



def read_mesh_shapefile(sFilename_mesh_in):
    """
    convert a shpefile to json format.
    This function should be used for stream flowline only.
    """

    aMesh=list()

    pDriver_json = ogr.GetDriverByName('GeoJSON')
    #pDriver_shapefile = ogr.GetDriverByName('ESRI Shapefile')
   
    pDataset_shapefile = pDriver_json.Open(sFilename_mesh_in, gdal.GA_ReadOnly)
    pLayer_shapefile = pDataset_shapefile.GetLayer(0)
    pSpatialRef_shapefile = pLayer_shapefile.GetSpatialRef()

    #we also need to spatial reference

    return aMesh, pSpatialRef_shapefile