import os, sys
import numpy as np
import osgeo
from osgeo import ogr, osr, gdal, gdalconst
from pyearth.gis.gdal.gdal_function import obtain_raster_metadata
from pyearth.gis.gdal.gdal_function import reproject_coordinates
from pyearth.gis.projection.degree_to_meter import degree_to_meter

from pystream.operation.create_mesh_op import create_mesh_op
from pystream.operation.preprocess_flowline_op import preprocess_flowline_op
from pystream.operation.intersect_flowline_with_mesh_with_postprocess_op import intersect_flowline_with_mesh_with_postprocess_op



def full_op(oModel_in):
    preprocess_flowline_op(oModel_in)
    create_mesh_op(oModel_in)
    intersect_flowline_with_mesh_with_postprocess_op(oModel_in)


    