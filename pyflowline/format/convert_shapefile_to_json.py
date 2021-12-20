import os
import json
from osgeo import ogr, osr, gdal, gdalconst
from pyflowline.format.read_flowline_shapefile import read_flowline_shapefile
def convert_shapefile_to_json(sFilename_shapefile_in, sFilename_geojson_in):

    #read shapefile into a list
    aFlowline_basin, pSpatialRef_pcs = read_flowline_shapefile( sFilename_shapefile_in )    
    #convert it

    iFlag_projected = 0
    pSpatialRef_gcs = osr.SpatialReference()
    pSpatialRef_gcs.ImportFromEPSG(4326)
    pSpatialRef_gcs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    export_flowline_to_json(iFlag_projected, aFlowline_basin, pSpatialRef_gcs, sFilename_geojson_in)


    return