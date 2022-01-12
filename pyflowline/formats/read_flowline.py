import os
import json
import numpy as np
from osgeo import ogr, osr, gdal, gdalconst
from shapely.geometry import Point, LineString, MultiLineString
from shapely.wkt import loads

from pyflowline.classes.edge import pyedge
from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.flowline import pyflowline

from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_flowline, convert_pcs_coordinates_to_flowline


def read_flowline_shapefile(sFilename_shapefile_in):
    """
    convert a shpefile to json format.
    This function should be used for stream flowline only.
    """

    aFlowline=list()

    #pDriver_json = ogr.GetDriverByName('GeoJSON')
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
                #pLine= LineString( aCoords[::-1 ] )

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
                #pLine= LineString( aCoords[::-1 ] )
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

    #pDriver_json = ogr.GetDriverByName('GeoJSON')
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
                #pLine= LineString( aCoords[::-1 ] )

                dummy1= np.array(aCoords)
                pLine = convert_gcs_coordinates_to_flowline(dummy1)
                pLine.lIndex = lID            
                aFlowline.append(pLine)
                lID = lID + 1
               
        else:
            if sGeometry_type =='LINESTRING':
                dummy = loads( pGeometry_in.ExportToWkt() )
                aCoords = dummy.coords
                #pLine= LineString( aCoords[::-1 ] )
                dummy1= np.array(aCoords)
                pLine = convert_gcs_coordinates_to_flowline(dummy1)
                pLine.lIndex = lID       
                aFlowline.append(pLine)
                lID = lID + 1
                
            else:
                print(sGeometry_type)
                pass
        
        
    
    #we also need to spatial reference

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

    lID = 0
    for pFeature_geojson in pLayer_geojson:
        pGeometry_geojson = pFeature_geojson.GetGeometryRef()
        pGeometry_in = pFeature_geojson.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()
        
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
        
        
    
    #we also need to spatial reference

    return aFlowline, pSpatialRef_geojson