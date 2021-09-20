import os
import json
from pystream.shared.edge import pyedge
from osgeo import ogr, osr, gdal, gdalconst
import numpy as np

from shapely.geometry import Point, LineString, MultiLineString
from shapely.wkt import loads

from pystream.shared.vertex import pyvertex
from pystream.shared.flowline import pyflowline

from pystream.format.convert_coordinates_to_flowline import convert_pcs_coordinates_to_flowline
from pystream.format.convert_coordinates_to_flowline import convert_gcs_coordinates_to_flowline

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

    pSpatialRef_gcs = osr.SpatialReference()
    pSpatialRef_gcs.ImportFromEPSG(4326)
    pSpatialRef_gcs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

    comparison = pSpatialRef_shapefile.IsSame(pSpatialRef_gcs)
    if(comparison != 1):
        iFlag_transform =1
        pTransform = osr.CoordinateTransformation(pSpatialRef_shapefile, pSpatialRef_gcs)
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