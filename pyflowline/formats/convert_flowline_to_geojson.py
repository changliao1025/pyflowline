import os
from osgeo import osr
from pyflowline.formats.read_flowline import read_flowline_shapefile, read_flowline_shapefile_swat, read_flowline_geojson
from pyflowline.formats.export_flowline import export_flowline_to_geojson

def convert_flowline_to_geojson(iFlag_type_in, sFilename_geojson_in, sFilename_geojson_out):
    """Convert a shapefile to a json file

    Args:
        iFlag_type_in (int): [0: vertex/point; 1: flowline/polyline; 2:polygon ]
        sFilename_shapefile_in ([type]): [description]
        sFilename_geojson_in ([type]): [description]
    """
    if os.path.isfile(sFilename_geojson_in):
        pass
    else:
        print('This shapefile does not exist: ', sFilename_geojson_in )
        iReturn_code = 0
        return iReturn_code


    if iFlag_type_in ==0:
        pass
    else:
        if iFlag_type_in == 1:
            
            aFlowline_basin, pSpatial_reference = read_flowline_geojson( sFilename_geojson_in )   
            #convert it

            iFlag_projected = 0
            pSpatial_reference_gcs = osr.SpatialReference()
            pSpatial_reference_gcs.ImportFromEPSG(4326)
            pSpatial_reference_gcs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
            export_flowline_to_geojson(aFlowline_basin, sFilename_geojson_out)
        else:
            if iFlag_type_in == 2:
                #polygon
                pass


    return
#this one need to be updated
def convert_shapefile_to_geojson(iFlag_type_in, sFilename_shapefile_in, sFilename_geojson_in):
    """Convert a shapefile to a json file

    Args:
        iFlag_type_in (int): [0: vertex/point; 1: flowline/polyline; 2:polygon ]
        sFilename_shapefile_in ([type]): [description]
        sFilename_geojson_in ([type]): [description]
    """
    if os.path.isfile(sFilename_shapefile_in):
        pass
    else:
        print('This shapefile does not exist: ', sFilename_shapefile_in )
        iReturn_code = 0
        return iReturn_code


    if iFlag_type_in ==0:
        pass
    else:
        if iFlag_type_in == 1:
            aFlowline_basin, pSpatial_reference = read_flowline_shapefile( sFilename_shapefile_in )    
            #aFlowline_basin, pSpatial_reference = read_flowline_geojson( sFilename_shapefile_in )   
            #convert it

            iFlag_projected = 0
            pSpatial_reference_gcs = osr.SpatialReference()
            pSpatial_reference_gcs.ImportFromEPSG(4326)
            pSpatial_reference_gcs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
            export_flowline_to_geojson(aFlowline_basin, sFilename_geojson_in)
        else:
            if iFlag_type_in == 2:
                #polygon
                pass


    return

def convert_shapefile_to_geojson_swat(iFlag_type_in, sFilename_shapefile_in, sFilename_geojson_in):
    """Convert a shapefile to a json file

    Args:
        iFlag_type_in (int): [0: vertex/point; 1: flowline/polyline; 2:polygon ]
        sFilename_shapefile_in ([type]): [description]
        sFilename_geojson_in ([type]): [description]
    """

    if iFlag_type_in ==0:
        pass
    else:
        if iFlag_type_in == 1:
            aFlowline_basin, pSpatial_reference = read_flowline_shapefile_swat( sFilename_shapefile_in )     
            iFlag_projected = 0
            pSpatial_reference_gcs = osr.SpatialReference()
            pSpatial_reference_gcs.ImportFromEPSG(4326)
            pSpatial_reference_gcs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
            export_flowline_to_geojson(aFlowline_basin, sFilename_geojson_in)
        else:
            if iFlag_type_in == 2:
                #polygon
                pass


    return