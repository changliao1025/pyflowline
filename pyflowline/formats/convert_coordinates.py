

from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.edge import pyedge
from pyflowline.classes.flowline import pyflowline

from pyflowline.classes.hexagon import pyhexagon
from pyflowline.classes.square import pysquare
from pyflowline.classes.latlon import pylatlon
from pyflowline.classes.mpas import pympas
from pyflowline.classes.tin import pytin

from pyflowline.algorithms.auxiliary.gdal_functions import reproject_coordinates

def convert_gcs_coordinates_to_cell(iMesh_type_in, \
    dLongitude_center_in, \
    dLatitude_center_in, \
        aCoordinates_gcs_in):

    npoint = len(aCoordinates_gcs_in)    
    aVertex=list()              
    aEdge=list()    
    for i in range(npoint-1):
        x = aCoordinates_gcs_in[i][0]
        y = aCoordinates_gcs_in[i][1]
        dummy = dict()
        dummy['dLongitude_degree'] = x
        dummy['dLatitude_degree'] = y
        pVertex = pyvertex(dummy)
        aVertex.append(pVertex)

    npoint2 = len(aVertex) 
    for j in range(npoint2-1):
        pEdge = pyedge( aVertex[j], aVertex[j+1] )
        aEdge.append(pEdge)
    
    #add the last one    
    pEdge = pyedge( aVertex[npoint2-1], aVertex[0] )
    aEdge.append(pEdge)

    if iMesh_type_in ==1: #hexagon
        

        pHexagon = pyhexagon( dLongitude_center_in, dLatitude_center_in, aEdge, aVertex)
        return pHexagon
    else:
        if iMesh_type_in ==2: #sqaure
            pSquare = pysquare(dLongitude_center_in, dLatitude_center_in, aEdge, aVertex)
            return pSquare
        else:
            if iMesh_type_in ==3: #latlon
                pLatlon = pylatlon(dLongitude_center_in, dLatitude_center_in, aEdge, aVertex)
                return pLatlon
            else:
                if iMesh_type_in ==4: #mpas       
                    pMpas = pympas(  dLongitude_center_in, dLatitude_center_in, aEdge, aVertex)
                    return pMpas
                else:
                    if iMesh_type_in ==5: #tin
                        pTin = pytin(dLongitude_center_in, dLatitude_center_in, aEdge, aVertex)
                        return pTin
                        pass
                    else:
                        print('What mesh type are you using?')
                        return None

def convert_pcs_coordinates_to_cell(iMesh_type_in, aCoordinates_pcs_in, pSpatial_reference_in):

    npoint = len(aCoordinates_pcs_in)    
    aVertex=list()              
    aEdge=list()    


    for i in range(npoint):
        x = aCoordinates_pcs_in[i][0]
        y = aCoordinates_pcs_in[i][1]
        dummy = dict()
        dummy['x'] = x
        dummy['y'] = y

        dummy['dLongitude_degree'], dummy['dLatitude_degree'] = reproject_coordinates(x, y , pSpatial_reference_in)
        pVertex = pyvertex(dummy)
        aVertex.append(pVertex)
    for j in range(npoint-1):
        pEdge = pyedge( aVertex[j], aVertex[j+1] )
        aEdge.append(pEdge)

    

    if iMesh_type_in ==1: #hexagon     

        pHexagon = pyhexagon(  aEdge, aVertex)
        return pHexagon
    else:
        if iMesh_type_in ==2: #sqaure
            pSquare = pysquare( aEdge, aVertex)
            return pSquare
        else:
            if iMesh_type_in ==3: #latlon
                pLatlon = pylatlon( aEdge, aVertex)
                return pLatlon
            else:
                if iMesh_type_in ==4: #mpas   

                    pMpas = pympas( aEdge, aVertex)
                    return pMpas
                else:
                    if iMesh_type_in ==5: #tin
                        pTin = pytin( aEdge, aVertex)
                        return pTin
                        pass
                    else:
                        print('What mesh type are you using?')
                        return None


def convert_gcs_coordinates_to_flowline(aCoordinates_in):
    npoint = len(aCoordinates_in)
    
    aVertex=list()
    for i in range(npoint):
        x = aCoordinates_in[i][0]
        y = aCoordinates_in[i][1]
        dummy = dict()
        dummy['dLongitude_degree'] = x
        dummy['dLatitude_degree'] = y
        pVertex = pyvertex(dummy)
        aVertex.append(pVertex)
        
    aEdge=list()
    for j in range(npoint-1):
        pEdge = pyedge( aVertex[j], aVertex[j+1] )
        aEdge.append(pEdge)
    
    pLine = pyflowline( aEdge)    
    return pLine

def convert_pcs_coordinates_to_flowline(aCoordinates_in, pSpatial_reference_in):
    npoint = len(aCoordinates_in)
    
    aVertex=list()
    for i in range(npoint):
        x = aCoordinates_in[i][0]
        y = aCoordinates_in[i][1]
        dummy = dict()
        dummy['x'] =x
        dummy['y'] =y
        lon, lat = reproject_coordinates(x, y, pSpatial_reference_in)
        dummy['lon'] = lon
        dummy['lat'] = lat
        pVertex = pyvertex(dummy)
        aVertex.append(pVertex)
        
    aEdge=list()
    for j in range(npoint-1):
        pEdge = pyedge( aVertex[j], aVertex[j+1] )
        aEdge.append(pEdge)
    
    pLine = pyflowline( aEdge)
    
    return pLine