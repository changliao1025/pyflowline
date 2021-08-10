import os
import json
from pystream.shared.edge import pyedge
from osgeo import ogr, osr, gdal, gdalconst
import numpy as np

from shapely.geometry import Point, LineString, MultiLineString
from shapely.wkt import loads

from pystream.shared.vertex import pyvertex
from pystream.shared.flowline import pyflowline



def read_mesh_shapefile(sFilename_shapefile_in):
    """
    convert a shpefile to json format.
    This function should be used for stream flowline only.
    """

    aMesh=list()

    pDriver_json = ogr.GetDriverByName('GeoJSON')
    pDriver_shapefile = ogr.GetDriverByName('ESRI Shapefile')
   
    pDataset_shapefile = pDriver_shapefile.Open(sFilename_shapefile_in, gdal.GA_ReadOnly)
    pLayer_shapefile = pDataset_shapefile.GetLayer(0)
    pSpatialRef_shapefile = pLayer_shapefile.GetSpatialRef()

    
    
        
        
    
    #we also need to spatial reference

    return aMesh, pSpatialRef_shapefile