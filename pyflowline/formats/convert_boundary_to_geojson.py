import os
from shutil import copy2
from pyearth.toolbox.data.shapefile.convert_shapefile_to_geojson import convert_shapefile_to_geojson
def convert_boundary_to_geojson(sFilename_in, sFilename_geojson_out, iFlag_merge_in = 0):

    if os.path.isfile(sFilename_in):
        #check the file type of the input file
        sExtension = os.path.splitext(sFilename_in)[1]
        if sExtension == '.shp':
            iFile_type = 1 #shapefile
        else:
            if sExtension == '.geojson':
                iFile_type = 2
            else:
                iFile_type = 0
        pass
    else:
        print('This input file does not exist: ', sFilename_in )
        iReturn_code = 0
        return iReturn_code

    if iFile_type == 1: #shapefile
       convert_shapefile_to_geojson(sFilename_in, sFilename_geojson_out , sLayername_in = 'boundary')
    else:
        if iFile_type == 2:
            #should check whether this file is single polygon or multiple polygons
            #copy the file
            copy2(sFilename_in, sFilename_geojson_out)
            pass

    return