
import numpy as np
from pyflowline.classes.edge import pyedge
from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.hexagon import pyhexagon
from pyflowline.classes.square import pysquare
from pyflowline.classes.latlon import pylatlon
from pyflowline.classes.mpas import pympas
from pyflowline.classes.tin import pytin

def convert_gcs_attributes_to_cell(iMesh_type_in,
                                   dLongitude_center_in,
                                   dLatitude_center_in,
                                   aCoordinates_gcs_in,
                                   aVertexID_in,
                                   aEdgeID_in,
                                   aVertexIndexOnEdge_in):
    """_summary_

    Args:
        iMesh_type_in (_type_): _description_
        dLongitude_center_in (_type_): _description_
        dLatitude_center_in (_type_): _description_
        aCoordinates_gcs_in (_type_): _description_
        aVertexID_in (_type_): _description_
        aEdgeID_in (_type_): _description_
        aVertexIndexOnEdge_in (_type_): _description_

    Returns:
        _type_: _description_
    """

    npoint = len(aVertexID_in)
    aVertex=list()
    aEdge=list()

    if iMesh_type_in == 1: #hexagon
        pHexagon = pyhexagon(dLongitude_center_in, dLatitude_center_in, aEdge, aVertex)
        return pHexagon
    else:
        if iMesh_type_in == 2: #sqaure
            pSquare = pysquare(dLongitude_center_in, dLatitude_center_in, aEdge, aVertex)
            return pSquare
        else:
            if iMesh_type_in == 3: #latlon
                pLatlon = pylatlon(dLongitude_center_in, dLatitude_center_in, aEdge, aVertex)
                return pLatlon
            else:
                if iMesh_type_in == 4: #mpas
                    for i in range(npoint):
                        lon = float(aCoordinates_gcs_in[i][0])
                        lat = float(aCoordinates_gcs_in[i][1])
                        pVertex = dict()
                        pVertex['dLongitude_degree'] =lon
                        pVertex['dLatitude_degree'] =lat
                        pVertex = pyvertex(pVertex)
                        pVertex.lVertexID = int(aVertexID_in[i])
                        aVertex.append(pVertex)

                    for j in range(npoint):
                        aVertexID_dummy = aVertexIndexOnEdge_in[j,:]
                        pVertex1 = dict()
                        dummy_index = np.where(aVertexID_in == aVertexID_dummy[0])
                        if len(dummy_index[0]) == 0 :
                            print('Vertex ID not found')
                            return None
                        pVertex1['dLongitude_degree'] = float(aCoordinates_gcs_in[dummy_index,0])
                        pVertex1['dLatitude_degree']  = float(aCoordinates_gcs_in[dummy_index,1])
                        pVertex_start = pyvertex(pVertex1)
                        pVertex_start.lVertexID = int( aVertexID_dummy[0] )

                        pVertex2 = dict()
                        dummy_index = np.where(aVertexID_in == aVertexID_dummy[1])
                        if len(dummy_index[0]) == 0 :
                            print('Vertex ID not found')
                            return None
                        pVertex2['dLongitude_degree'] = float(aCoordinates_gcs_in[dummy_index,0])
                        pVertex2['dLatitude_degree']  = float(aCoordinates_gcs_in[dummy_index,1])
                        pVertex_end = pyvertex(pVertex2)
                        pVertex_end.lVertexID = int(aVertexID_dummy[1])
                        pEdge = pyedge( pVertex_start, pVertex_end )
                        pEdge.lEdgeID = int(aEdgeID_in[j])
                        aEdge.append(pEdge)

                    pMpas = pympas(dLongitude_center_in, dLatitude_center_in, aEdge, aVertex)
                    return pMpas
                else:
                    if iMesh_type_in ==5: #tin
                        pTin = pytin(dLongitude_center_in, dLatitude_center_in, aEdge, aVertex)
                        return pTin

                    else:
                        print('What mesh type are you using?')
                        return None


def convert_pcs_attributes_to_cell(iMesh_type_in,
    aCoordinates_pcs_in,
    aVertexID_in,
    aEdgeID_in,
    aVertexIndexOnEdge_in):

    npoint = len(aVertexID_in)
    aVertex=list()
    aEdge=list()

    if iMesh_type_in ==1: #hexagon
        pHexagon = pyhexagon( aEdge, aVertex)
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
                    for i in range(npoint):
                        lon = aCoordinates_pcs_in[i][0]
                        lat = aCoordinates_pcs_in[i][1]
                        pVertex = dict()
                        pVertex['dLongitude_degree'] =lon
                        pVertex['dLatitude_degree'] =lat
                        pVertex = pyvertex(pVertex)
                        pVertex.lVertexID = int(aVertexID_in[i])
                        aVertex.append(pVertex)

                    for j in range(npoint):
                        aVertexID_dummy = aVertexIndexOnEdge_in[j,:]
                        pVertex = dict()
                        dummy_index = np.where(aVertexID_in == aVertexID_dummy[0])
                        pVertex['dLongitude_degree'] =aCoordinates_pcs_in[dummy_index,0]
                        pVertex['dLatitude_degree'] =aCoordinates_pcs_in[dummy_index,1]
                        pVertex_start = pyvertex(pVertex)
                        pVertex_start.lVertexID = int( aVertexID_dummy[0] )
                        dummy_index = np.where(aVertexID_in == aVertexID_dummy[1])
                        pVertex['dLongitude_degree'] =aCoordinates_pcs_in[dummy_index,0]
                        pVertex['dLatitude_degree'] =aCoordinates_pcs_in[dummy_index,1]
                        pVertex_end = pyvertex(pVertex)
                        pVertex_end.lVertexID = int(aVertexID_dummy[1])
                        pEdge = pyedge( pVertex_start, pVertex_end )
                        pEdge.lEdgeID = int(aEdgeID_in[j])
                        aEdge.append(pEdge)

                    pMpas = pympas( aEdge, aVertex)
                    return pMpas
                else:
                    if iMesh_type_in ==5: #tin
                        pTin = pytin( aEdge, aVertex)
                        return pTin

                    else:
                        print('What mesh type are you using?')
                        return None


