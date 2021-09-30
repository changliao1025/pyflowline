import os
import json
from pyflowline.shared.edge import pyedge
from osgeo import ogr, osr, gdal, gdalconst
import numpy as np

from shapely.geometry import Point, LineString, MultiLineString
from shapely.wkt import loads

from pyflowline.shared.vertex import pyvertex
from pyflowline.shared.flowline import pyflowline

from pyflowline.format.convert_coordinates_to_flowline import convert_gcs_coordinates_to_flowline, convert_pcs_coordinates_to_flowline

def read_flowline_geojson(sFilename_geojson_in):
    """
    convert a shpefile to json format.
    This function should be used for stream flowline only.
    """

    aFlowline=list()

    pDriver_geojson =  ogr.GetDriverByName('ESRI Shapefile')# ogr.GetDriverByName('GeoJSON')
   
    pDataset_geojson = pDriver_geojson.Open(sFilename_geojson_in, gdal.GA_ReadOnly)
    pLayer_geojson = pDataset_geojson.GetLayer(0)
    pSpatialRef_geojson = pLayer_geojson.GetSpatialRef()

    lID = 0
    for pFeature_geojson in pLayer_geojson:
        pGeometry_geojson = pFeature_geojson.GetGeometryRef()
        pGeometry_in = pFeature_geojson.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()
        if(sGeometry_type == 'MULTILINESTRING'):
            aLine = ogr.ForceToLineString(pGeometry_in)
            for Line in aLine: 
                dummy = loads( Line.ExportToWkt() )
                aCoords = dummy.coords
                #pLine= LineString( aCoords[::-1 ] )

                dummy1= np.array(aCoords)
                pLine = convert_pcs_coordinates_to_flowline(dummy1)
                pLine.lIndex = lID
                aFlowline.append(pLine)
                lID = lID + 1
               
        else:
            if sGeometry_type =='LINESTRING':
                dummy = loads( pGeometry_in.ExportToWkt() )
                aCoords = dummy.coords
                #pLine= LineString( aCoords[::-1 ] )
                dummy1= np.array(aCoords)
                pLine = convert_pcs_coordinates_to_flowline(dummy1)
                pLine.lIndex = lID
                aFlowline.append(pLine)
                lID = lID + 1
                
            else:
                print(sGeometry_type)
                pass
        
        
    
    #we also need to spatial reference

    return aFlowline, pSpatialRef_geojson