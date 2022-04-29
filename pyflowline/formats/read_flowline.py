import os
import numpy as np
from osgeo import ogr, osr, gdal

from shapely.wkt import loads

from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_flowline

def read_flowline_shapefile(sFilename_shapefile_in):
    """
    convert a shpefile to json format.
    This function should be used for stream flowline only.
    """
    iReturn_code = 1
    if os.path.isfile(sFilename_shapefile_in):
        pass
    else:
        print('This shapefile does not exist: ', sFilename_shapefile_in )
        iReturn_code = 0
        return iReturn_code

    aFlowline=list()
    pDriver_shapefile = ogr.GetDriverByName('ESRI Shapefile')   
    pDataset_shapefile = pDriver_shapefile.Open(sFilename_shapefile_in, gdal.GA_ReadOnly)
    pLayer_shapefile = pDataset_shapefile.GetLayer(0)
    pSpatialRef_shapefile = pLayer_shapefile.GetSpatialRef()
    pSpatial_reference_gcs = osr.SpatialReference()
    pSpatial_reference_gcs.ImportFromEPSG(4326)
    pSpatial_reference_gcs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    comparison = pSpatialRef_shapefile.IsSame(pSpatial_reference_gcs)
    if(comparison != 1):
        iFlag_transform =1
        pTransform = osr.CoordinateTransformation(pSpatialRef_shapefile, pSpatial_reference_gcs)
    else:
        iFlag_transform =0
   
    lID = 0
    for pFeature_shapefile in pLayer_shapefile:        
        pGeometry_in = pFeature_shapefile.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()
        lNHDPlusID = int(pFeature_shapefile.GetField("NHDPlusID"))
        if (iFlag_transform ==1): #projections are different
            pGeometry_in.Transform(pTransform)
        if (pGeometry_in.IsValid()):
            pass
        else:
            print('Geometry issue')

        if(sGeometry_type == 'MULTILINESTRING'):
            aLine = ogr.ForceToLineString(pGeometry_in)
            for Line in aLine: 
                dummy = loads( Line.ExportToWkt() )
                aCoords = dummy.coords                
                dummy1= np.array(aCoords)
                pLine = convert_gcs_coordinates_to_flowline(dummy1)
                pLine.lIndex = lID
                pLine.lNHDPlusID= lNHDPlusID
                aFlowline.append(pLine)
                lID = lID + 1
               
        else:
            if sGeometry_type =='LINESTRING':
                dummy = loads( pGeometry_in.ExportToWkt() )
                aCoords = dummy.coords                
                dummy1= np.array(aCoords)
                pLine = convert_gcs_coordinates_to_flowline(dummy1)
                pLine.lIndex = lID
                pLine.lNHDPlusID= lNHDPlusID
                aFlowline.append(pLine)
                lID = lID + 1
                
            else:
                print(sGeometry_type)
                pass
        
        
    
    #we also need to spatial reference

    return aFlowline, pSpatialRef_shapefile

def read_flowline_shapefile_swat(sFilename_shapefile_in):
    """
    convert a shpefile to json format.
    This function should be used for stream flowline only.
    """
    aFlowline=list()
    pDriver_shapefile = ogr.GetDriverByName('ESRI Shapefile')   
    pDataset_shapefile = pDriver_shapefile.Open(sFilename_shapefile_in, gdal.GA_ReadOnly)
    pLayer_shapefile = pDataset_shapefile.GetLayer(0)
    pSpatialRef_shapefile = pLayer_shapefile.GetSpatialRef()
    pSpatial_reference_gcs = osr.SpatialReference()
    pSpatial_reference_gcs.ImportFromEPSG(4326)
    pSpatial_reference_gcs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    comparison = pSpatialRef_shapefile.IsSame(pSpatial_reference_gcs)
    if(comparison != 1):
        iFlag_transform =1
        pTransform = osr.CoordinateTransformation(pSpatialRef_shapefile, pSpatial_reference_gcs)
    else:
        iFlag_transform =0

    lID = 0
    for pFeature_shapefile in pLayer_shapefile:        
        pGeometry_in = pFeature_shapefile.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()
        if (iFlag_transform ==1): #projections are different
            pGeometry_in.Transform(pTransform)

        if (pGeometry_in.IsValid()):
            pass
        else:
            print('Geometry issue')

        if(sGeometry_type == 'MULTILINESTRING'):
            aLine = ogr.ForceToLineString(pGeometry_in)
            for Line in aLine: 
                dummy = loads( Line.ExportToWkt() )
                aCoords = dummy.coords      
                dummy1= np.array(aCoords)
                pLine = convert_gcs_coordinates_to_flowline(dummy1)
                pLine.lIndex = lID            
                aFlowline.append(pLine)
                lID = lID + 1
               
        else:
            if sGeometry_type =='LINESTRING':
                dummy = loads( pGeometry_in.ExportToWkt() )
                aCoords = dummy.coords               
                dummy1= np.array(aCoords)
                pLine = convert_gcs_coordinates_to_flowline(dummy1)
                pLine.lIndex = lID       
                aFlowline.append(pLine)
                lID = lID + 1
                
            else:
                print(sGeometry_type)
                pass
    

    return aFlowline, pSpatialRef_shapefile


def read_flowline_geojson(sFilename_geojson_in):
    """
    read a geojson flowline
    This function should be used for stream flowline only.
    """

    aFlowline=list()
    pDriver_geojson = ogr.GetDriverByName('GeoJSON')   
    pDataset_geojson = pDriver_geojson.Open(sFilename_geojson_in, gdal.GA_ReadOnly)
    pLayer_geojson = pDataset_geojson.GetLayer(0)
    pSpatialRef_geojson = pLayer_geojson.GetSpatialRef()
    ldefn = pLayer_geojson.GetLayerDefn()
    schema =list()
    for n in range(ldefn.GetFieldCount()):
        fdefn = ldefn.GetFieldDefn(n)
        schema.append(fdefn.name)
    if 'iseg' in schema:
        iFlag_segment = 1
    else:
        iFlag_segment = 0    
    if 'id' in schema:
        iFlag_id = 1
    else:
        iFlag_id = 0

    lID = 0
    for pFeature_geojson in pLayer_geojson:
        pGeometry_geojson = pFeature_geojson.GetGeometryRef()
        pGeometry_in = pFeature_geojson.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()
        if iFlag_segment ==1:
            iStream_segment = pFeature_geojson.GetField("iseg")
        else:
            iStream_segment = -1
        
        if iFlag_id ==1:
            lFlowlineID = pFeature_geojson.GetField("id")
        else:
            lFlowlineID = -1
        
        if sGeometry_type =='LINESTRING':
            dummy = loads( pGeometry_in.ExportToWkt() )
            aCoords = dummy.coords
            dummy1= np.array(aCoords)
            pLine = convert_gcs_coordinates_to_flowline(dummy1)
            pLine.lIndex = lID
            pLine.iStream_segment = iStream_segment
            pLine.lFlowlineID = lFlowlineID
            aFlowline.append(pLine)
            lID = lID + 1            
        else:
            print(sGeometry_type)
            pass        


    return aFlowline, pSpatialRef_geojson