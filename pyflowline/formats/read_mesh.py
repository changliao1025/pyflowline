import os
import json
import numpy as np
from osgeo import ogr, osr, gdal, gdalconst
from shapely.geometry import Point, LineString, MultiLineString
from shapely.wkt import loads

from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.edge import pyedge
from pyflowline.classes.flowline import pyflowline



def read_mesh_json(sFilename_mesh_in):
    """
    convert a shpefile to json format.
    This function should be used for stream flowline only.
    """

    aMesh_out=list()
    pDriver_json = ogr.GetDriverByName('GeoJSON') 
   
    pDataset_shapefile = pDriver_json.Open(sFilename_mesh_in, gdal.GA_ReadOnly)
    pLayer_shapefile = pDataset_shapefile.GetLayer(0)
    pSpatial_reference_out = pLayer_shapefile.GetSpatialRef()

    #we also need to spatial reference

    return aMesh_out, pSpatial_reference_out



def read_mesh_shapefile(sFilename_mesh_in):
    """
    convert a shpefile to json format.
    This function should be used for stream flowline only.
    """
    aMesh_out=list()

    
    pDriver_shapefile = ogr.GetDriverByName('ESRI Shapefile')
   
    pDataset_shapefile = pDriver_shapefile.Open(sFilename_mesh_in, gdal.GA_ReadOnly)
    pLayer_shapefile = pDataset_shapefile.GetLayer(0)
    pSpatial_reference_out = pLayer_shapefile.GetSpatialRef()

    #we also need to spatial reference

    return aMesh_out, pSpatial_reference_out