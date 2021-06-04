

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
def build_hexagon_grid_topology(aHexagon_in):
    
    

    nHexagon = len(aHexagon_in)

    
    lID_mesh = 0
    lID_flowline =0 
    #create mesh list first

    #use confluence information
    for i in range (nHexagon):
        
        

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

                        #link
                        #1:find start and end of flowlie
                        pVertex_start = pLine.pVertex_start
                        pVertex_end = pLine.pVertex_end
                        pHexagon_start = pHexagon
                        #the start vertex must be on the enterance edge
                        iFlag_found, aEdge_shared = aHexagon.which_edge_cross_this_vertex(pVertex_start)
                        if iFlag_found ==1:
                            #find the other hexagon
                            aHexagon_shared = find_hexagon_through_edge(aHexagon, aEdge_shared)
                            pass
                        else:
                            pass
                        
                    
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

            #now add back to the cell object
            pHexagon.aFlowline = aFlowline_intersect
            pHexagon.nFlowline = len(aFlowline_intersect)
            aHexagon.append(pHexagon)
            lID_mesh = lID_mesh + 1   
        else:
            pass

    
    return aHexagon