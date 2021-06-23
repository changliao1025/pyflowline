
import imp
import os
import json
from pystream.shared.vertex import pyvertex
import numpy as np
from osgeo import ogr, osr

from shapely.geometry import Point, LineString, MultiLineString
from shapely.wkt import loads
from pystream.shared.hexagon import pyhexagon

from pystream.format.convert_coordinates_to_hexagon import convert_coordinates_to_hexagon
from pystream.format.convert_coordinates_to_flowline import convert_coordinates_to_flowline
from pystream.find_hexagon_through_edge import find_hexagon_through_edge

from pystream.shared.link import pyhexagonlink
def intersect_flowline_with_mesh(sFilename_mesh, sFilename_flowline, sFilename_output):

    if  os.path.exists(sFilename_mesh) and  os.path.exists(sFilename_flowline) : 
        pass
    else:
        print('The input file does not exist')
        return

    if os.path.exists(sFilename_output): 
        #delete it if it exists
        os.remove(sFilename_output)

    sDriverName = "GeoJSON"
    pDriver = ogr.GetDriverByName( sDriverName )

    #geojson
    aHexagon=list()
    aHexagon_intersect=list()
    
   
    pDataset_mesh = pDriver.Open(sFilename_mesh, 0)
    pDataset_flowline = pDriver.Open(sFilename_flowline, 0)   

    pLayer_mesh = pDataset_mesh.GetLayer(0)
    pSpatialRef_mesh = pLayer_mesh.GetSpatialRef()
    nfeature_mesh = pLayer_mesh.GetFeatureCount()
    print( pSpatialRef_mesh)
    pLayer_flowline = pDataset_flowline.GetLayer(0)
    pSpatialRef_flowline = pLayer_flowline.GetSpatialRef()
    nfeature_flowline = pLayer_flowline.GetFeatureCount()
    
    print( pSpatialRef_flowline)
    comparison = pSpatialRef_mesh.IsSame(pSpatialRef_flowline)
    if(comparison != 1):
        iFlag_transform =1
        transform = osr.CoordinateTransformation(pSpatialRef_mesh, pSpatialRef_flowline)
    else:
        iFlag_transform =0

    pDataset_out = pDriver.CreateDataSource(sFilename_output)

    pLayerOut = pDataset_out.CreateLayer('flowline', pSpatialRef_flowline, ogr.wkbMultiLineString)
    # Add one attribute
    pLayerOut.CreateField(ogr.FieldDefn('id', ogr.OFTInteger64)) #long type for high resolution
    
    pLayerDefn = pLayerOut.GetLayerDefn()
    pFeatureOut = ogr.Feature(pLayerDefn)

    
    lID_mesh = 0
    lID_flowline =0 
    
        


    for i in range (nfeature_mesh):
    #for pFeature_mesh in pLayer_mesh:       
        pFeature_mesh= pLayer_mesh.GetFeature(i)
        pGeometry_mesh = pFeature_mesh.GetGeometryRef()
        if (iFlag_transform ==1): #projections are different
            pGeometry_mesh.Transform(transform)

        if (pGeometry_mesh.IsValid()):
            pass
        else:
            print('Geometry issue')

        #convert geometry to edge
        pGeometrytype_mesh = pGeometry_mesh.GetGeometryName()
        if(pGeometrytype_mesh == 'POLYGON'):
            dummy = loads( pGeometry_mesh.ExportToWkt() )
            aCoords = dummy.exterior.coords
            dummy1= np.array(aCoords)
            pHexagon = convert_coordinates_to_hexagon(dummy1)
            pHexagon.lIndex = lID_mesh
                     
            #print(pGeometry_mesh.GetGeometryName())
            aFlowline_intersect = list()
            iFlag_intersected = 0 
            for j in range (nfeature_flowline):
            #for pFeature_flowline in pLayer_flowline:
                pFeature_flowline = pLayer_flowline.GetFeature(j)
                pGeometry_flowline = pFeature_flowline.GetGeometryRef()

                if (pGeometry_flowline.IsValid()):
                    pass
                else:
                    print('Geometry issue')
                #print(pGeometry_flowline.GetGeometryName())

                iFlag_intersect = pGeometry_flowline.Intersects( pGeometry_mesh )
                if( iFlag_intersect == True):

                    iFlag_intersected = 1
                    pGeometry_intersect = pGeometry_flowline.Intersection(pGeometry_mesh) 
                    pFeatureOut.SetGeometry(pGeometry_intersect)
                    pFeatureOut.SetField("id", lID_flowline)                
                    pLayerOut.CreateFeature(pFeatureOut)    

                    #add more process here to 
                    pGeometrytype_intersect = pGeometry_intersect.GetGeometryName()
                    if pGeometrytype_intersect == 'LINESTRING':

                        dummy = loads( pGeometry_intersect.ExportToWkt() )
                        aCoords = dummy.coords                
                        dummy1= np.array(aCoords)
                        pLine = convert_coordinates_to_flowline(dummy1)
                        pLine.lIndex = lID_flowline
                        aFlowline_intersect.append(pLine)
                        lID_flowline = lID_flowline + 1
                    
                    else:
                        if(pGeometrytype_intersect == 'MULTILINESTRING'):
                            aLine = ogr.ForceToLineString(pGeometry_intersect)
                            for Line in aLine: 
                                dummy = loads( Line.ExportToWkt() )
                                aCoords = dummy.coords
                                dummy1= np.array(aCoords)
                                pLine = convert_coordinates_to_flowline(dummy1)
                                pLine.lIndex = lID_flowline
                                aFlowline_intersect.append(pLine)
                                lID_flowline = lID_flowline + 1
                            pass
                        else:
                            pass
                    

                    

                else:
                    pass

            #only save the intersected hexagon to output? 
            #now add back to the cell object
            pHexagon.aFlowline = aFlowline_intersect
            pHexagon.nFlowline = len(aFlowline_intersect)
            if iFlag_intersected ==1:     
                pHexagon.iFlag_intersected = 1                       
                aHexagon_intersect.append(pHexagon)
            else:
                pHexagon.iFlag_intersected = 0   
                pass


            
            aHexagon.append(pHexagon)
            lID_mesh = lID_mesh + 1   
        else:
            pass

    
    return aHexagon, aHexagon_intersect