import os
import json
from pyflowline.classes.square import pysquare
from pyflowline.classes.edge import pyedge
from osgeo import ogr, osr, gdal, gdalconst
import numpy as np

from shapely.geometry import Point, LineString, MultiLineString
from shapely.wkt import loads

from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.flowline import pyflowline
from pyflowline.classes.hexagon import pyhexagon
from pyflowline.classes.square import pysquare
from pyflowline.classes.latlon import pylatlon
from pyflowline.classes.mpas import pympas
from pyflowline.classes.tin import pytin
    
def convert_gcs_attribute_to_cell(iMesh_type, dLon, dLat, \
        aCoordinates_gcs, \
        aVertexID, \
        aEdgeID, \
        aVertexIndexOnEdge):  
        
    npoint = len(aVertexID)     
    aVertex=list()        
    aEdge=list()
    
    if iMesh_type == 1: #hexagon
        pHexagon = pyhexagon(dLon, dLat, aEdge, aVertex)
        return pHexagon
    else:
        if iMesh_type == 2: #sqaure
            pSquare = pysquare(dLon, dLat, aEdge, aVertex)
            return pSquare
        else:
            if iMesh_type == 3: #latlon
                pLatlon = pylatlon(dLon, dLat, aEdge, aVertex)
                return pLatlon
            else:
                if iMesh_type == 4: #mpas
                    for i in range(npoint):
                        lon = aCoordinates_gcs[i][0]
                        lat = aCoordinates_gcs[i][1]
                        pVertex = dict()        
                        pVertex['lon'] =lon
                        pVertex['lat'] =lat
                        pVertex = pyvertex(pVertex)
                        pVertex.lVertexID = int(aVertexID[i])
                        aVertex.append(pVertex)

                    for j in range(npoint):
                        aVertexID_dummy = aVertexIndexOnEdge[j,:]
                        pVertex = dict()
                        dummy_index = np.where(aVertexID == aVertexID_dummy[0])     
                        pVertex['lon'] =aCoordinates_gcs[dummy_index,0]
                        pVertex['lat'] =aCoordinates_gcs[dummy_index,1]
                        pVertex_start = pyvertex(pVertex)
                        pVertex_start.lVertexID = int( aVertexID_dummy[0] )
                        dummy_index = np.where(aVertexID == aVertexID_dummy[1])     
                        pVertex['lon'] =aCoordinates_gcs[dummy_index,0]
                        pVertex['lat'] =aCoordinates_gcs[dummy_index,1]
                        pVertex_end = pyvertex(pVertex)
                        pVertex_end.lVertexID = int(aVertexID_dummy[1])
                        pEdge = pyedge( pVertex_start, pVertex_end )
                        pEdge.lEdgeID = int(aEdgeID[j])
                        aEdge.append(pEdge)

                    pMpas = pympas(dLon, dLat, aEdge, aVertex)
                    return pMpas
                else:
                    if iMesh_type ==5: #tin
                        pTin = pytin(dLon, dLat, aEdge, aVertex)
                        return pTin
                        pass
                    else:
                        print('What mesh type are you using?')
                        return None
                    
  
def convert_pcs_attribute_to_cell(iMesh_type, aVertexID, aEdgeID,aVertexIndexOnEdge, aCoordinates_pcs):  
    npoint = len(aVertexID)     
    aVertex=list()        
    aEdge=list()
    
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
                    for i in range(npoint):
                        lon = aCoordinates_pcs[i][0]
                        lat = aCoordinates_pcs[i][1]
                        pVertex = dict()        
                        pVertex['lon'] =lon
                        pVertex['lat'] =lat
                        pVertex = pyvertex(pVertex)
                        pVertex.lVertexID = int(aVertexID[i])
                        aVertex.append(pVertex)

                    for j in range(npoint):
                        aVertexID_dummy = aVertexIndexOnEdge[j,:]
                        pVertex = dict()
                        dummy_index = np.where(aVertexID == aVertexID_dummy[0])     
                        pVertex['lon'] =aCoordinates_pcs[dummy_index,0]
                        pVertex['lat'] =aCoordinates_pcs[dummy_index,1]
                        pVertex_start = pyvertex(pVertex)
                        pVertex_start.lVertexID = int( aVertexID_dummy[0] )
                        dummy_index = np.where(aVertexID == aVertexID_dummy[1])     
                        pVertex['lon'] =aCoordinates_pcs[dummy_index,0]
                        pVertex['lat'] =aCoordinates_pcs[dummy_index,1]
                        pVertex_end = pyvertex(pVertex)
                        pVertex_end.lVertexID = int(aVertexID_dummy[1])
                        pEdge = pyedge( pVertex_start, pVertex_end )
                        pEdge.lEdgeID = int(aEdgeID[j])
                        aEdge.append(pEdge)

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

