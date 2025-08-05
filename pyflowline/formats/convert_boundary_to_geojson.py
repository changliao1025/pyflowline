import os
from shutil import copy2
from osgeo import ogr, gdal, osr
from pyearth.toolbox.data.shapefile.convert_shapefile_to_geojson import convert_shapefile_to_geojson
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area
def convert_boundary_to_geojson(sFilename_in, sFilename_geojson_out, iFlag_largest_polygon_in = 0):
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

    if iFlag_largest_polygon_in == 1:
        #first convert to a temporary geojson file
        sFilename_tmp = os.path.splitext(sFilename_geojson_out)[0] + '_tmp.geojson'
        if iFile_type == 1: #shapefile
           convert_shapefile_to_geojson(sFilename_in, sFilename_tmp , sLayername_in = 'boundary')
        else:
            if iFile_type == 2:
                #now check whether it is a multi-polygon or single polygon
                pDriver_json = ogr.GetDriverByName('GeoJSON')
                pDataset_boundary = pDriver_json.Open(sFilename_in, gdal.GA_ReadOnly)
                pLayer_boundary = pDataset_boundary.GetLayer(0)
                pFeature_boundary= pLayer_boundary.GetNextFeature()
                #we also need to spatial reference
                dArea_max = 0.0
                while pFeature_boundary:
                    pGeometry_tmp = pFeature_boundary.GetGeometryRef()
                    pGeometrytype_boundary = pGeometry_tmp.GetGeometryName()
                    if(pGeometrytype_boundary == 'POLYGON'):
                        aCoords_gcs = get_geometry_coordinates(pGeometry_tmp)
                        aLon=list()
                        aLat=list()
                        for aCoord in aCoords_gcs:
                            aLon.append(aCoord[0])
                            aLat.append(aCoord[1])
                        dArea = calculate_polygon_area(aLon, aLat)
                        if dArea > dArea_max:
                            dArea_max = dArea
                            pGeometry_out = pGeometry_tmp.Clone()

                    else:
                        if(pGeometrytype_boundary == 'MULTIPOLYGON'):
                            nPolygon = pGeometry_tmp.GetGeometryCount()
                            for i in range(nPolygon):
                                pGeometry = pGeometry_tmp.GetGeometryRef(i)
                                aCoords_gcs = get_geometry_coordinates(pGeometry)
                                aLon=list()
                                aLat=list()
                                for aCoord in aCoords_gcs:
                                    aLon.append(aCoord[0])
                                    aLat.append(aCoord[1])
                                dArea = calculate_polygon_area(aLon, aLat)
                                if dArea > dArea_max:
                                    dArea_max = dArea
                                    pGeometry_out = pGeometry.Clone()
                                    pass
                        else:
                            print('This geometry type is not supported: ', pGeometrytype_boundary)
                            return 0

                    pFeature_boundary = pLayer_boundary.GetNextFeature()

                #export the geometry to geojson
                if os.path.isfile(sFilename_geojson_out):
                    os.remove(sFilename_geojson_out)
                    pass
                pDataset_geojson = pDriver_json.CreateDataSource(sFilename_geojson_out)
                pSrs = osr.SpatialReference()
                pSrs.ImportFromEPSG(4326)  # WGS84
                pLayer_geojson = pDataset_geojson.CreateLayer('boundary', geom_type=pGeometry_out.GetGeometryType(), srs=pSrs)
                pLayer_geojson.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
                pFeature_geojson = ogr.Feature(pLayer_geojson.GetLayerDefn())
                pFeature_geojson.SetGeometry(pGeometry_out)
                pFeature_geojson.SetField('id', 1)
                pLayer_geojson.CreateFeature(pFeature_geojson)
                pFeature_geojson = None
                pLayer_geojson = None
                pDataset_geojson = None
    else:
        if iFile_type == 1: #shapefile
           convert_shapefile_to_geojson(sFilename_in, sFilename_geojson_out , sLayername_in = 'boundary')
        else:
            if iFile_type == 2:
                #should check whether this file is single polygon or multiple polygons
                #copy the file
                copy2(sFilename_in, sFilename_geojson_out)
                pass

    return