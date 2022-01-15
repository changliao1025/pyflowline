import os
import json
import numpy as np
from osgeo import ogr, osr, gdal, gdalconst
from shapely.geometry import Point, LineString, MultiLineString
from shapely.wkt import loads

from pyflowline.classes.edge import pyedge
from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.flowline import pyflowline

from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_conceptual, convert_pcs_coordinates_to_conceptual



def calculate_area_of_difference_raw(sFilename_a, sFilename_b):

    


    return

def calculate_area_of_difference_simplified(sFilename_simplified_in, sFilename_conceptual_in):

    pDriver_geojson = ogr.GetDriverByName( "GeoJSON")
    aCell=list()
    aCell_intersect=list()    

    pDataset_simplified = pDriver_geojson.Open(sFilename_simplified_in, 0)
    pDataset_conceptual = pDriver_geojson.Open(sFilename_conceptual_in, 0)   
    pLayer_simplified = pDataset_simplified.GetLayer(0)
    pSpatial_reference_simplified = pLayer_simplified.GetSpatialRef()
    nfeature_simplified = pLayer_simplified.GetFeatureCount()

    pLayer_conceptual = pDataset_conceptual.GetLayer(0)
    pSpatial_reference_conceptual = pLayer_conceptual.GetSpatialRef()
    nfeature_conceptual = pLayer_conceptual.GetFeatureCount()
    pLayerDefinition = pLayer_conceptual.GetLayerDefn()
    
    #print( pSpatial_reference_conceptual)
    comparison = pSpatial_reference_simplified.IsSame
    if(comparison != 1):
        iFlag_transform = 1
        transform = osr.CoordinateTransformation(pSpatial_reference_simplified, pSpatial_reference_conceptual)
    else:
        iFlag_transform = 0
    

    for pFeature_simplified in pLayer_simplified:
       
        #pFeature_simplified= pLayer_mesh.GetFeature(i)
        pGeometry_simplified = pFeature_simplified.GetGeometryRef()        
        dummy0 = loads( pGeometry_simplified.ExportToWkt() )
        aCoords_gcs = dummy0.coords
        aCoords_gcs= np.array(aCoords_gcs)       

        
        if (iFlag_transform ==1): #projections are different
            pGeometry_simplified.Transform(transform)

        if (pGeometry_simplified.IsValid()):
            pass
        else:
            print('Geometry issue')
        pGeometrytype_simplified = pGeometry_simplified.GetGeometryName()
        if(pGeometrytype_simplified == 'LINESTRING'):            
            
                     
            aFlowline_intersect = list()
            iFlag_intersected = 0 

            for j in range (nfeature_conceptual):
        
                pFeature_conceptual = pLayer_conceptual.GetFeature(j)
                pGeometry_conceptual = pFeature_conceptual.GetGeometryRef()

             

                if (pGeometry_conceptual.IsValid()):
                    pass
                else:
                    print('Geometry issue')
                #print(pGeometry_conceptual.GetGeometryName())

                iFlag_intersect = pGeometry_conceptual.Intersects( pGeometry_simplified )
                if( iFlag_intersect == True):

                    iFlag_intersected = 1
                    pGeometry_intersect = pGeometry_conceptual.Intersection(pGeometry_simplified)                     

                    #add more process here to 
                    pGeometrytype_intersect = pGeometry_intersect.GetGeometryName()
                    if pGeometrytype_intersect == 'LINESTRING':
                        pFeatureOut.SetGeometry(pGeometry_intersect)
                        pFeatureOut.SetField("id", lID_flowline)         
                        pFeatureOut.SetField("iseg", iStream_segment)    
                        pFeatureOut.SetField("iord", iStream_order)           
                        pLayerOut.CreateFeature(pFeatureOut)    

                        dummy = loads( pGeometry_intersect.ExportToWkt() )
                        aCoords = dummy.coords                
                        dummy1= np.array(aCoords)
                        pLine = convert_gcs_coordinates_to_flowline(dummy1)
                        pLine.calculate_length()
                        pLine.lIndex = lID_flowline
                        pLine.iStream_segment = iStream_segment
                        pLine.iStream_order = iStream_order
                        aFlowline_intersect.append(pLine)
                        aFlowline_intersect_all.append(pLine)
                        lID_flowline = lID_flowline + 1
                    
                    else:
                        if(pGeometrytype_intersect == 'MULTILINESTRING'):
                            aLine = ogr.ForceToLineString(pGeometry_intersect)
                            for Line in aLine: 
                                pFeatureOut.SetGeometry(Line)
                                pFeatureOut.SetField("id", lID_flowline)         
                                pFeatureOut.SetField("iseg", iStream_segment)    
                                pFeatureOut.SetField("iord", iStream_order)           
                                pLayerOut.CreateFeature(pFeatureOut)    

                                dummy = loads( Line.ExportToWkt() )
                                aCoords = dummy.coords
                                dummy1= np.array(aCoords)
                                pLine = convert_gcs_coordinates_to_flowline(dummy1)
                                pLine.calculate_length()
                                pLine.lIndex = lID_flowline
                                pLine.iStream_segment = iStream_segment
                                pLine.iStream_order = iStream_order
                                aFlowline_intersect.append(pLine)
                                aFlowline_intersect_all.append(pLine)
                                lID_flowline = lID_flowline + 1
                            pass
                        else:
                            pass
                            

                else:
                    pass