import os
import json
from osgeo import ogr, osr, gdal, gdalconst
from pyflowline.formats.read_flowline_shapefile import read_flowline_shapefile
from pyflowline.formats.export_flowline_to_json import export_flowline_to_json

def convert_shapefile_to_json(iFlag_type_in, sFilename_shapefile_in, sFilename_geojson_in):

    #read shapefile into a list
    aFlowline_basin, pSpatial_reference = read_flowline_shapefile( sFilename_shapefile_in )    
    #convert it

    iFlag_projected = 0
    pSpatialRef_gcs = osr.SpatialReference()
    pSpatialRef_gcs.ImportFromEPSG(4326)
    pSpatialRef_gcs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    export_flowline_to_json(iFlag_projected, aFlowline_basin, pSpatialRef_gcs, sFilename_geojson_in)


    return