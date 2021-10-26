import os
import json
from pyflowline.shared.square import pysquare
from pyflowline.shared.edge import pyedge
from osgeo import ogr, osr, gdal, gdalconst
import numpy as np

from shapely.geometry import Point, LineString, MultiLineString
from shapely.wkt import loads

from pyflowline.shared.vertex import pyvertex
from pyflowline.shared.flowline import pyflowline
from pyflowline.shared.hexagon import pyhexagon
from pyflowline.shared.square import pysquare
from pyflowline.shared.latlon import pylatlon
from pyflowline.shared.mpas import pympas
from pyflowline.shared.tin import pytin

from pyearth.gis.gdal.gdal_function import reproject_coordinates



def convert_gcs_coordinates_to_cell(iMesh_type, aCoordinates_gcs, dLon, dLat):
    npoint = len(aCoordinates_gcs)    
    aVertex=list()              
    aEdge=list()    
    for i in range(npoint-1):
        x = aCoordinates_gcs[i][0]
        y = aCoordinates_gcs[i][1]
        dummy = dict()
        dummy['lon'] = x
        dummy['lat'] = y
        pVertex = pyvertex(dummy)
        aVertex.append(pVertex)

    npoint2 = len(aVertex) 
    for j in range(npoint2-1):
        pEdge = pyedge( aVertex[j], aVertex[j+1] )
        aEdge.append(pEdge)
    
    #add the last one    
    pEdge = pyedge( aVertex[npoint2-1], aVertex[0] )
    aEdge.append(pEdge)

    if iMesh_type ==1: #hexagon
        

        pHexagon = pyhexagon( aEdge, aVertex)
        return pHexagon
    else:
        if iMesh_type ==2: #sqaure
            pSquare = pysquare( aEdge, aVertex)
            return pSquare
        else:
            if iMesh_type ==3: #latlon
                pLatlon = pylatlon( aEdge, aVertex)
                return pLatlon
            else:
                if iMesh_type ==4: #mpas       
                    pMpas = pympas(  dLon, dLat, aEdge, aVertex)
                    return pMpas
                else:
                    if iMesh_type ==5: #tin
                        pTin = pytin( aEdge, aVertex)
                        return pTin
                        pass
                    else:
                        print('What mesh type are you using?')
                        return None

def convert_pcs_coordinates_to_cell(iMesh_type, aCoordinates_pcs, pSpatialRef_in):
    npoint = len(aCoordinates_pcs)    
    aVertex=list()              
    aEdge=list()    


    for i in range(npoint):
        x = aCoordinates_pcs[i][0]
        y = aCoordinates_pcs[i][1]
        dummy = dict()
        dummy['x'] = x
        dummy['y'] = y

        dummy['lon'], dummy['lat'] = reproject_coordinates(x, y , pSpatialRef_in)
        pVertex = pyvertex(dummy)
        aVertex.append(pVertex)
    for j in range(npoint-1):
        pEdge = pyedge( aVertex[j], aVertex[j+1] )
        aEdge.append(pEdge)

    if iMesh_type ==1: #hexagon     

        pHexagon = pyhexagon( aEdge, aVertex)
        return pHexagon
    else:
        if iMesh_type ==2: #sqaure
            pSquare = pysquare( aEdge, aVertex)
            return pSquare
        else:
            if iMesh_type ==3: #latlon
                pLatlon = pylatlon( aEdge, aVertex)
                return pLatlon
            else:
                if iMesh_type ==4: #mpas   

                    pMpas = pympas( aEdge, aVertex)
                    return pMpas
                else:
                    if iMesh_type ==5: #tin
                        pTin = pytin( aEdge, aVertex)
                        return pTin
                        pass
                    else:
                        print('What mesh type are you using?')
                        return None

