import os
import json
from pystream.shared.edge import pyedge
from osgeo import ogr, osr, gdal, gdalconst
import numpy as np

from shapely.geometry import Point, LineString, MultiLineString
from shapely.wkt import loads

from pystream.shared.vertex import pyvertex
from pystream.shared.flowline import pyflowline

def convert_coordinates_to_flowline(aCoordinates):
    npoint = len(aCoordinates)
    
    aVertex=list()
    for i in range(npoint):
        x = aCoordinates[i][0]
        y = aCoordinates[i][1]
        dummy = dict()
        dummy['x'] =x
        dummy['y'] =y
        pVertex = pyvertex(dummy)
        aVertex.append(pVertex)
        
    aEdge=list()
    for j in range(npoint-1):
        pEdge = pyedge( aVertex[j], aVertex[j+1] )
        aEdge.append(pEdge)
    
    pLine = pyflowline( aEdge)
    
    return pLine

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
                dummy = loads( Line.ExportToWkt() )
                aCoords = dummy.coords
                #pLine= LineString( aCoords[::-1 ] )

                dummy1= np.array(aCoords)
                pLine = convert_coordinates_to_flowline(dummy1)
                pLine.lIndex = lID
                aFlowline.append(pLine)
                lID = lID + 1
               
        else:
            if sGeometry_type =='LINESTRING':
                dummy = loads( pGeometry_in.ExportToWkt() )
                aCoords = dummy.coords
                #pLine= LineString( aCoords[::-1 ] )
                dummy1= np.array(aCoords)
                pLine = convert_coordinates_to_flowline(dummy1)
                pLine.lIndex = lID
                aFlowline.append(pLine)
                lID = lID + 1
                
            else:
                print(sGeometry_type)
                pass
        
        
    
    #we also need to spatial reference

    return aFlowline, pSpatialRef_shapefile