import os
import json
from osgeo import ogr, osr, gdal, gdalconst

from shapely.geometry import Point, LineString, MultiLineString
from shapely.wkt import loads

from pystream.shared.flowline import pyflowline

def read_flowline_shapefile(sFilename_shapefile_in):
    """
    convert a shpefile to json format.
    This function should be used for stream flowline only.
    """

    aFlowline=list()

    pDriver = ogr.GetDriverByName('GeoJSON')
    pDriver_shapefile = ogr.GetDriverByName('ESRI Shapefile')
   
    pDataset_shapefile = pDriver_shapefile.Open(sFilename_shapefile_in, gdal.GA_ReadOnly)
    pLayer_shapefile = pDataset_shapefile.GetLayer(0)
    pSpatialRef_shapefile = pLayer_shapefile.GetSpatialRef()

    lID = 0
    for pFeature_shapefile in pLayer_shapefile:
        pGeometry_shapefile = pFeature_shapefile.GetGeometryRef()
        pGeometry_in = pFeature_shapefile.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()
        if(sGeometry_type == 'MULTILINESTRING'):
            aLine = ogr.ForceToLineString(pGeometry_in)
            for Line in aLine: 
                dummy = loads( pGeometry_in.ExportToWkt() )
                aCoords = dummy.coords
                #pLine= LineString( aCoords[::-1 ] )

                #aEdge = 
                pLine = pyflowline( aCoords)
                pLine.lIndex = lID


                aFlowline.append(pLine)
                lID = lID + 1
               
        else:
            if sGeometry_type =='LINESTRING':
                dummy = loads( pGeometry_in.ExportToWkt() )
                aCoords = dummy.coords
                #pLine= LineString( aCoords[::-1 ] )
                pLine = pyflowline( aCoords)
                pLine.lIndex = lID
                aFlowline.append(pLine)
                lID = lID + 1
                
            else:
                print(sGeometry_type)
                pass
        
        
    
    #we also need to spatial reference

    return aFlowline, pSpatialRef_shapefile