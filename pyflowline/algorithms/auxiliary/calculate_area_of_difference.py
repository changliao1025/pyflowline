import os
import json
from telnetlib import X3PAD
import numpy as np
from osgeo import ogr, osr, gdal, gdalconst
from shapely.geometry import Point, LineString, MultiLineString
from shapely.wkt import loads

from shapely.ops import polygonize, polygonize_full
from pyflowline.algorithms.auxiliary.gdal_functions import calculate_angle_betwen_vertex,  calculate_angle_betwen_vertex_normal
from pyflowline.classes.edge import pyedge
from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.flowline import pyflowline

from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_flowline 

from pyflowline.algorithms.intersect.intersect_flowline_with_flowline import intersect_flowline_with_flowline
from pyflowline.algorithms.auxiliary.gdal_functions import calculate_polygon_area



def calculate_area_of_difference_raw(sFilename_a, sFilename_b):

    #not yet supported


    return

def calculate_area_of_difference_simplified(aFlowline_in, aVertex_all_in, \
     sFilename_output_in):

    if os.path.exists(sFilename_output_in): 
        #delete it if it exists
        os.remove(sFilename_output_in)
        pass

    pDriver_geojson = ogr.GetDriverByName( "GeoJSON")

    nFlowline = len(aFlowline_in)
    nVeretx = len(aVertex_all_in)
    aFlowline = list()


    #rebuild
    for i  in range(nFlowline):
        aFlowline_in[i].lFlowlineID = i
    for i  in range(nVeretx):
        aVertex_all_in[i].lVertexID = i

    for i  in range(nVeretx):
    #for n, v in enumerate(aVertex_all_in):
    
        for j in range(nFlowline):
            if aVertex_all_in[i] == aFlowline_in[j].pVertex_start:
                aFlowline_in[j].pVertex_start.lVertexID = i

            elif aVertex_all_in[i] == aFlowline_in[j].pVertex_end:
                aFlowline_in[j].pVertex_end.lVertexID = i



    for pFlowline in aFlowline_in:      

        for j  in range(len(aFlowline_in)):
            pFlowline2 = aFlowline_in[j]
            if pFlowline.lFlowlineID != pFlowline2.lFlowlineID:
                if pFlowline.pVertex_start in (pFlowline2.pVertex_start, pFlowline2.pVertex_end):
                    pFlowline.aFlowlineID_start.append(pFlowline2.lFlowlineID)                    
                elif pFlowline.pVertex_end in (pFlowline2.pVertex_start, pFlowline2.pVertex_end):                                    
                    pFlowline.aFlowlineID_end.append(pFlowline2.lFlowlineID)

  


    def get_right_branch(iFlag_reverse, pFlowline_in):
        """
        This function starts at a contact (called `edge_0`) and calculates
        the angles between the adjacent contacts at a side (either `stop` or
        `start`) and then chooses the right-branching contact, and the opposite
        side for the next contact.

        For example, if we are at `edge_0` and `edge_2` is the right-branching contact,
        and the `start` side of `edge_2` meets `edge_0`, then the next side will be `stop`.

        """

        nEdge = pFlowline_in.nEdge
        if iFlag_reverse ==0 :
            #normal direction
            pVertex_stop = pFlowline_in.pVertex_end
            pEdge = pFlowline_in.aEdge[nEdge-1]

            x1 = pEdge.pVertex_start.dLongitude_degree
            y1 = pEdge.pVertex_start.dLatitude_degree
            x2 = pEdge.pVertex_end.dLongitude_degree
            y2 = pEdge.pVertex_end.dLatitude_degree
            
            angle_min = 360
            iIndex_right = -1
            for i in pFlowline_in.aFlowlineID_end:
                pFlowline_dummy = aFlowline_in[i]
                if (pFlowline_dummy.pVertex_start == pVertex_stop):
                    nEdge2 = pFlowline_dummy.nEdge
                    pEdge2 = pFlowline_dummy.aEdge[0]
                    pVertex_dummy = pEdge2.pVertex_end
                else:
                    #revered
                    nEdge2 = pFlowline_dummy.nEdge
                    pEdge2 = pFlowline_dummy.aEdge[nEdge2-1]
                    pVertex_dummy = pEdge2.pVertex_start
                
                #calculate angle              
                
                x3 = pVertex_dummy.dLongitude_degree
                y3 = pVertex_dummy.dLatitude_degree

                #angle_dummy0 = calculate_angle_betwen_vertex( x1, y1, x2, y2, x3, y3  )

                angle_dummy = calculate_angle_betwen_vertex_normal( x1, y1, x2, y2, x3, y3  )
                print(angle_dummy)
                if angle_dummy < angle_min:
                    angle_min = angle_dummy
                    iIndex_right = i
            
            #mini

            pFlowline_out = aFlowline_in[iIndex_right]

            if pFlowline_out.pVertex_start == pVertex_stop:
                #aFlowline_in[iIndex_right].iFlag_right =1
                iFlag_reverse_new=0
                pVertex_stop_out =   pFlowline_out.pVertex_end
            else:
                #aFlowline_in[iIndex_right].iFlag_left =1
                iFlag_reverse_new=1
                pVertex_stop_out=pFlowline_out.pVertex_start
            return iFlag_reverse_new, pFlowline_out, pVertex_stop_out


        else:
            #rever direction
            pVertex_stop = pFlowline_in.pVertex_start

        
    def get_left_branch(iFlag_reverse, pFlowline_in):
        """
        This function starts at a contact (called `edge_0`) and calculates
        the angles between the adjacent contacts at a side (either `stop` or
        `start`) and then chooses the right-branching contact, and the opposite
        side for the next contact.

        For example, if we are at `edge_0` and `edge_2` is the right-branching contact,
        and the `start` side of `edge_2` meets `edge_0`, then the next side will be `stop`.

        """

        nEdge = pFlowline_in.nEdge
        if iFlag_reverse ==0 :
            #normal direction
            pVertex_stop = pFlowline_in.pVertex_end
            pEdge = pFlowline_in.aEdge[nEdge-1]

            x1 = pEdge.pVertex_start.dLongitude_degree
            y1 = pEdge.pVertex_start.dLatitude_degree
            x2 = pEdge.pVertex_end.dLongitude_degree
            y2 = pEdge.pVertex_end.dLatitude_degree
            
            angle_min = 360
            iIndex_right = -1
            for i in pFlowline_in.aFlowlineID_end:
                pFlowline_dummy = aFlowline_in[i]
                if (pFlowline_dummy.pVertex_start == pVertex_stop):
                    pVertex_start_dummy = pFlowline_dummy.pVertex_start
                else:
                    pVertex_start_dummy = pFlowline_dummy.pVertex_end
                
                #calculate angle
                
                
                x3 = pVertex_start_dummy.dLongitude_degree
                y3 = pVertex_start_dummy.dLatitude_degree

                angle_dummy = calculate_angle_betwen_vertex( x1, y1, x2, y2, x3,y3  )
                print(angle_dummy)
                if angle_dummy < angle_min:
                    angle_min = angle_dummy
                    iIndex_right = i
            
            #take the mini

            


        else:
            #rever direction
            pVertex_stop = pFlowline_in.pVertex_start

    def walk_cycle(iFlag_reverse, pFlowline_in ):

        if iFlag_reverse ==0:

            iFlag_right = pFlowline_in.iFlag_right
            if iFlag_right == 1:
                return
            else:
                pVertex_start = pFlowline_in.pVertex_start
                pVertex_end = pFlowline_in.pVertex_end
                iFlag_reverse_new, pFlowline_next, pVertex_stop  = get_right_branch(0, pFlowline_in)
                aFlowline_in[pFlowline_in.lFlowlineID].iFlag_right=1

                if (pVertex_stop == pVertex_start):
                    #a loop is finished
                    pass
                else:
                    pFlowline_in = pFlowline_next
                    walk_cycle( iFlag_reverse_new, pFlowline_in)

        else:
            iFlag_left = pFlowline_in.iFlag_left

    
        return

    

    for i in range(nFlowline):
        pFlowline_in = aFlowline_in[i]
        walk_cycle(0, pFlowline_in)



        
    
    #shapely method, but not complete
    #for i in range(nFlowline):
    #   
    #    pFlowline = aFlowline_in[i]
    #    nVertex= pFlowline.nVertex
    #    aCoords_gcs = np.full((nVertex,2), -9999. ,dtype=float)
    #    for k in range(nVertex):
    #        aCoords_gcs[k,0] = pFlowline.aVertex[k].dLongitude_degree
    #        aCoords_gcs[k,1] = pFlowline.aVertex[k].dLatitude_degree
    #    
    #    aCoords_gcs= tuple(j for j in aCoords_gcs) #np.array(aCoords_gcs)    
    #    aFlowline.append(aCoords_gcs)   
    #dummy = polygonize(aFlowline)
    #aPolygon_out = list(dummy)


    pDataset = pDriver_geojson.CreateDataSource(sFilename_output_in) 
    pSpatial_reference_gcs = osr.SpatialReference()  
    pSpatial_reference_gcs.ImportFromEPSG(4326)    # WGS84 lat/lon
    pSpatial_reference_gcs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

    pLayer = pDataset.CreateLayer('intersect', pSpatial_reference_gcs, ogr.wkbPolygon)
    # Add one attribute
    pLayer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger64)) #long type for high resolution

    pLayerDefn = pLayer.GetLayerDefn()
    pFeature = ogr.Feature(pLayerDefn)    
    lCellID =0
    dArea =0.0

    
    for po in aPolygon_out:
        ring = ogr.Geometry(ogr.wkbLinearRing)
        aCoords_gcs = po.exterior.coords
        aCoords_gcs= np.array(aCoords_gcs)  
        nPoint  = (aCoords_gcs.shape)[0]
        lons=list()
        lats=list()
        for i in range(nPoint):
            ring.AddPoint(aCoords_gcs[i,0], aCoords_gcs[i,1])
            lons.append( aCoords_gcs[i,0] )
            lats.append( aCoords_gcs[i,1] )
        
        pPolygon = ogr.Geometry(ogr.wkbPolygon)
        pPolygon.AddGeometry(ring)
        pFeature.SetGeometry(pPolygon)
        pFeature.SetField("id", lCellID)
        pLayer.CreateFeature(pFeature)
        lCellID= lCellID+1

        dArea0 = calculate_polygon_area(lons, lats)
        dArea = dArea + dArea0

    return aPolygon_out, dArea
   
    

    

      

             

    
