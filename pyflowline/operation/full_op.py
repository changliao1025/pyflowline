import os, sys
import numpy as np
import osgeo
from osgeo import ogr, osr, gdal, gdalconst
from pyearth.gis.gdal.gdal_function import obtain_raster_metadata
from pyearth.gis.gdal.gdal_function import reproject_coordinates
from pyearth.gis.projection.degree_to_meter import degree_to_meter

from pyflowline.operation.create_mesh_op import create_mesh_op
from pyflowline.operation.preprocess_flowline_op import preprocess_flowline_op
from pyflowline.operation.intersect_flowline_with_mesh_with_postprocess_op import intersect_flowline_with_mesh_with_postprocess_op



def full_op(opyflowline_in):
    iFlag_simplification = opyflowline_in.iFlag_simplification
    iFlag_create_mesh = opyflowline_in.iFlag_create_mesh
    iFlag_intersect = opyflowline_in.iFlag_intersect
    
    if iFlag_simplification == 1:
        preprocess_flowline_op(opyflowline_in)

    if iFlag_create_mesh==1:
        create_mesh_op(opyflowline_in)

    if iFlag_intersect ==1:
        
        intersect_flowline_with_mesh_with_postprocess_op(opyflowline_in)


    