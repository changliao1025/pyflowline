import importlib
from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.edge import pyedge
from pyflowline.classes.flowline import pyflowline

from pyflowline.classes.hexagon import pyhexagon
from pyflowline.classes.square import pysquare
from pyflowline.classes.latlon import pylatlon
from pyflowline.classes.mpas import pympas
from pyflowline.classes.tin import pytin
from pyflowline.classes.dggrid import pydggrid

from pyearth.gis.spatialref.reproject_coodinates import reproject_coordinates
iFlag_cython = importlib.util.find_spec("cython") 
if iFlag_cython is not None:
    from pyflowline.algorithms.cython.kernel import calculate_angle_betwen_vertex
    from pyflowline.algorithms.cython.kernel import calculate_distance_to_plane
else:   
    from pyearth.gis.geometry.calculate_angle_betwen_vertex import calculate_angle_betwen_vertex
    from pyearth.gis.geometry.calculate_distance_to_plane import calculate_distance_to_plane

def convert_gcs_coordinates_to_cell(iMesh_type_in,     
                                    dLongitude_center_in,     
                                    dLatitude_center_in,         
                                    aCoordinates_gcs_in,        
                                    iFlag_simplify_in=None):
    
    
   

    if iFlag_simplify_in is None:
        iFlag_simplify_in = 0
    else:
        iFlag_simplify_in = iFlag_simplify_in

    #the closed polygon has a duplicate point (start and end are the same)
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

    if iFlag_simplify_in ==1:
        aVertex_simple=list()  

        for i in range(npoint-1):
            if i == 0 : #the first one
                pv_start = aVertex[npoint-2]
                pv_middle = aVertex[i]
                pv_end = aVertex[i+1]                
            else:
                if i == npoint-2: #the last one
                    pv_start = aVertex[i-1]
                    pv_middle = aVertex[i]
                    pv_end = aVertex[0]
                else:
                    pv_start = aVertex[i-1]
                    pv_middle = aVertex[i]
                    pv_end = aVertex[i+1]
                    
            #calculate the angle between the three points
            x1 = pv_start.dLongitude_degree
            y1 = pv_start.dLatitude_degree
            x2 = pv_middle.dLongitude_degree
            y2 = pv_middle.dLatitude_degree
            x3 = pv_end.dLongitude_degree
            y3 = pv_end.dLatitude_degree
            angle3deg = calculate_angle_betwen_vertex(x1,y1, x2,y2, x3,y3)
            if  angle3deg > 175: #care
                pass
            else:
                #the center is a actual vertex
                aVertex_simple.append(pv_middle)
                
        #replace the original vertex
        aVertex = aVertex_simple
        pass
    else:
        pass

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
                    if iMesh_type_in ==5: #dggrid
                        pdggrid = pydggrid(dLongitude_center_in, dLatitude_center_in, aEdge, aVertex)
                        return pdggrid
                       
                    else:
                        if iMesh_type_in ==6 : #tin
                            pTin = pytin(dLongitude_center_in, dLatitude_center_in, aEdge, aVertex)
                            return pTin
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
                    if iMesh_type_in ==5: #dggrid
                        pdggrid = pydggrid( aEdge, aVertex)
                        return pdggrid
                 
                    else:
                        if iMesh_type_in ==6: #tin
                            pTin = pytin( aEdge, aVertex)
                            return pTin
                        else:
                            print('What mesh type are you using?')
                            return None


def convert_gcs_coordinates_to_flowline(aCoordinates_in):
    """convert coordinates to flowline, but we cannot setup index and id yet

    Args:
        aCoordinates_in (_type_): _description_

    Returns:
        _type_: _description_
    """
    
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
        if aVertex[j] == aVertex[j+1]:
            print('Two vertices are the same')
            pass
        else:
            pEdge = pyedge( aVertex[j], aVertex[j+1] )
            aEdge.append(pEdge)
    
    if len(aEdge) == 0:
        print('No edge is created')
        pFlowline = None
        pass
    else:
        pFlowline = pyflowline( aEdge)    
    return pFlowline

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
    
    pFlowline = pyflowline( aEdge)
    
    return pFlowline