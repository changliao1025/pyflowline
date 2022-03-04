import os
import numpy as np
from osgeo import ogr, osr
from shapely.ops import polygonize
from pyflowline.algorithms.auxiliary.gdal_functions import  calculate_angle_betwen_vertex_normal
from pyflowline.algorithms.auxiliary.gdal_functions import calculate_polygon_area
from pyflowline.algorithms.auxiliary.find_index_in_list import find_list_in_list ,find_vertex_in_list

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
                if (pFlowline.pVertex_start == pFlowline2.pVertex_start ):
                    pFlowline.aFlowlineID_start_start.append(pFlowline2.lFlowlineID)

                if (pFlowline.pVertex_start==pFlowline2.pVertex_end):
                    pFlowline.aFlowlineID_start_end.append(pFlowline2.lFlowlineID)   
                
                if (pFlowline.pVertex_end==pFlowline2.pVertex_start ):
                    pFlowline.aFlowlineID_end_start.append(pFlowline2.lFlowlineID)
                if (pFlowline.pVertex_end==pFlowline2.pVertex_end):
                    pFlowline.aFlowlineID_end_end.append(pFlowline2.lFlowlineID)

    def get_next_branch(iFlag_rightleft, iFlag_reverse, pFlowline_in):
        """
        This function starts at a contact (called `edge_0`) and calculates
        the angles between the adjacent contacts at a side (either `stop` or
        `start`) and then chooses the right-branching contact, and the opposite
        side for the next contact.

        For example, if we are at `edge_0` and `edge_2` is the right-branching contact,
        and the `start` side of `edge_2` meets `edge_0`, then the next side will be `stop`.

        """

        nEdge = pFlowline_in.nEdge
        iFlag_exist = 0 
        if iFlag_reverse ==0 :
            #normal direction
            pVertex_stop = pFlowline_in.pVertex_end
            pEdge = pFlowline_in.aEdge[nEdge-1]

            x1 = pEdge.pVertex_start.dLongitude_degree
            y1 = pEdge.pVertex_start.dLatitude_degree
            x2 = pEdge.pVertex_end.dLongitude_degree
            y2 = pEdge.pVertex_end.dLatitude_degree
            
           
            iIndex_right = -1

            n = len(pFlowline_in.aFlowlineID_end_start) + len(pFlowline_in.aFlowlineID_end_end)
            if n == 0:
                iFlag_exist = 0 
                iFlag_reverse_new= 0
                pFlowline_out = None
                pVertex_stop_out = None
                return iFlag_exist, iFlag_reverse_new, pFlowline_out, pVertex_stop_out
            else:
                aAngle=list()
                for i in pFlowline_in.aFlowlineID_end_start:
                    pFlowline_dummy = aFlowline_in[i]
                    if (pFlowline_dummy.pVertex_start == pVertex_stop):
                        nEdge2 = pFlowline_dummy.nEdge
                        pEdge2 = pFlowline_dummy.aEdge[0]
                        pVertex_dummy = pEdge2.pVertex_end
                    else:                        #revered
                        nEdge2 = pFlowline_dummy.nEdge
                        pEdge2 = pFlowline_dummy.aEdge[nEdge2-1]
                        pVertex_dummy = pEdge2.pVertex_start

                    #calculate angle              
                    x3 = pVertex_dummy.dLongitude_degree
                    y3 = pVertex_dummy.dLatitude_degree
                    angle_dummy = calculate_angle_betwen_vertex_normal( x1, y1, x2, y2, x3, y3  )                    
                    aAngle.append(angle_dummy)
                for i in pFlowline_in.aFlowlineID_end_end:
                    pFlowline_dummy = aFlowline_in[i]
                    if (pFlowline_dummy.pVertex_start == pVertex_stop):
                        nEdge2 = pFlowline_dummy.nEdge
                        pEdge2 = pFlowline_dummy.aEdge[0]
                        pVertex_dummy = pEdge2.pVertex_end
                    else:                        #revered
                        nEdge2 = pFlowline_dummy.nEdge
                        pEdge2 = pFlowline_dummy.aEdge[nEdge2-1]
                        pVertex_dummy = pEdge2.pVertex_start

                    #calculate angle              
                    x3 = pVertex_dummy.dLongitude_degree
                    y3 = pVertex_dummy.dLatitude_degree
                    angle_dummy = calculate_angle_betwen_vertex_normal( x1, y1, x2, y2, x3, y3  )                    
                    aAngle.append(angle_dummy)
                    

                #mini
                dummy = pFlowline_in.aFlowlineID_end_start + pFlowline_in.aFlowlineID_end_end
                if iFlag_rightleft ==0 :
                    min_angle = min(aAngle)
                    iIndex_right = aAngle.index(min_angle)                   
                    pFlowline_out = aFlowline_in[ dummy[iIndex_right] ]
                    pass
                else:
                    max_angle = max(aAngle)
                    iIndex_left = aAngle.index(max_angle)
                    pFlowline_out = aFlowline_in[ dummy[iIndex_left]]
                    pass

                if pFlowline_out.pVertex_start == pVertex_stop:                 
                    iFlag_reverse_new=0
                    pVertex_stop_out =   pFlowline_out.pVertex_end
                else:                  
                    iFlag_reverse_new=1
                    pVertex_stop_out=pFlowline_out.pVertex_start

                iFlag_exist = 1
                return iFlag_exist, iFlag_reverse_new, pFlowline_out, pVertex_stop_out


        else:
            #rever direction
            pVertex_stop = pFlowline_in.pVertex_start
            pEdge = pFlowline_in.aEdge[0]            
            x1 = pEdge.pVertex_end.dLongitude_degree
            y1 = pEdge.pVertex_end.dLatitude_degree
            x2 = pEdge.pVertex_start.dLongitude_degree
            y2 = pEdge.pVertex_start.dLatitude_degree
            angle_min = 360
            iIndex_right = -1
            n = len(pFlowline_in.aFlowlineID_start_start) + len(pFlowline_in.aFlowlineID_start_end)
            if n == 0:
                iFlag_exist = 0 
                iFlag_reverse_new= 0
                pFlowline_out = None
                pVertex_stop_out = None
                return iFlag_exist, iFlag_reverse_new, pFlowline_out, pVertex_stop_out
            else:
                aAngle=list()
                for i in pFlowline_in.aFlowlineID_start_start:
                    pFlowline_dummy = aFlowline_in[i]
                    if (pFlowline_dummy.pVertex_start == pVertex_stop):
                        nEdge2 = pFlowline_dummy.nEdge
                        pEdge2 = pFlowline_dummy.aEdge[0]
                        pVertex_dummy = pEdge2.pVertex_end
                    else:
                        nEdge2 = pFlowline_dummy.nEdge
                        pEdge2 = pFlowline_dummy.aEdge[nEdge2-1]
                        pVertex_dummy = pEdge2.pVertex_start

                    x3 = pVertex_dummy.dLongitude_degree
                    y3 = pVertex_dummy.dLatitude_degree
                    angle_dummy = calculate_angle_betwen_vertex_normal( x1, y1, x2, y2, x3, y3  )
                    aAngle.append(angle_dummy)
                for i in pFlowline_in.aFlowlineID_start_end:
                    pFlowline_dummy = aFlowline_in[i]
                    if (pFlowline_dummy.pVertex_start == pVertex_stop):
                        nEdge2 = pFlowline_dummy.nEdge
                        pEdge2 = pFlowline_dummy.aEdge[0]
                        pVertex_dummy = pEdge2.pVertex_end
                    else:
                        nEdge2 = pFlowline_dummy.nEdge
                        pEdge2 = pFlowline_dummy.aEdge[nEdge2-1]
                        pVertex_dummy = pEdge2.pVertex_start

                    x3 = pVertex_dummy.dLongitude_degree
                    y3 = pVertex_dummy.dLatitude_degree
                    angle_dummy = calculate_angle_betwen_vertex_normal( x1, y1, x2, y2, x3, y3  )
                    aAngle.append(angle_dummy)
                    
                #mini
                dummy = pFlowline_in.aFlowlineID_start_start + pFlowline_in.aFlowlineID_start_end
                if iFlag_rightleft ==0 :
                    min_angle = min(aAngle)
                    iIndex_right = aAngle.index(min_angle)
                    pFlowline_out = aFlowline_in[dummy[iIndex_right]]
                    pass
                else:
                    max_angle = max(aAngle)
                    iIndex_left = aAngle.index(max_angle)
                    pFlowline_out = aFlowline_in[dummy[iIndex_left]]
                    pass

                if pFlowline_out.pVertex_start == pVertex_stop:            
                    iFlag_reverse_new=0
                    pVertex_stop_out =   pFlowline_out.pVertex_end
                else:           
                    iFlag_reverse_new=1
                    pVertex_stop_out=pFlowline_out.pVertex_start

                iFlag_exist = 1
                return iFlag_exist, iFlag_reverse_new, pFlowline_out, pVertex_stop_out
  

    def walk_cycle(iFlag_rightleft, iFlag_reverse, pFlowline_in, pVertex_origin, aFlowline_list, aVertex_list ):
        iFlag_loop = 1
        if iFlag_reverse ==0:

            iFlag_right = pFlowline_in.iFlag_right
            if iFlag_right == 1:
                iFlag_loop = 0
                return iFlag_loop, aFlowline_list
            else:
                
                iFlag_exist, iFlag_reverse_new, pFlowline_next, pVertex_stop  = get_next_branch(iFlag_rightleft, 0, pFlowline_in)
                if iFlag_exist ==1:
                    aFlowline_list.append(pFlowline_next.lFlowlineID)
                    if (pVertex_stop == pVertex_origin):
                        #a loop is finished
                        iFlag_loop = 1
                        return iFlag_loop, aFlowline_list
                        pass
                    else:
                        iFlag_vertex_exist, index_dummy = find_vertex_in_list(aVertex_list, pVertex_stop)
                        if iFlag_vertex_exist ==1:
                            #this one show before
                            iFlag_loop = 0
                            return iFlag_loop, aFlowline_list
                        else:
                            aVertex_list.append(pVertex_stop)

                        pFlowline_in = pFlowline_next
                        iFlag_loop, aFlowline_list = walk_cycle(iFlag_rightleft, iFlag_reverse_new, pFlowline_in,\
                            pVertex_origin, aFlowline_list, aVertex_list)
                else:
                    iFlag_loop = 0
                    return iFlag_loop, aFlowline_list

        else:
            iFlag_right = pFlowline_in.iFlag_right
            if iFlag_right == 1:
                iFlag_loop = 0
                return iFlag_loop, aFlowline_list
            else:
                
                iFlag_exist, iFlag_reverse_new, pFlowline_next, pVertex_stop = get_next_branch(iFlag_rightleft,1, pFlowline_in)
                
                if iFlag_exist ==1:
                    #aFlowline_in[pFlowline_in.lFlowlineID].iFlag_right=1

                    aFlowline_list.append(pFlowline_next.lFlowlineID)

                    if (pVertex_stop == pVertex_origin):
                        #a loop is finished
                        iFlag_loop = 1
                        return iFlag_loop, aFlowline_list                        
                    else:
                        iFlag_vertex_exist, index_dummy = find_vertex_in_list(aVertex_list, pVertex_stop)
                        if iFlag_vertex_exist ==1:
                            #this one show before
                            iFlag_loop = 0
                            return iFlag_loop, aFlowline_list
                        else:
                            aVertex_list.append(pVertex_stop)

                        pFlowline_in = pFlowline_next
                        iFlag_loop, aFlowline_list = walk_cycle(iFlag_rightleft, iFlag_reverse_new, pFlowline_in, \
                        pVertex_origin, aFlowline_list,aVertex_list)
                else:
                    iFlag_loop = 0
                    return iFlag_loop, aFlowline_list

                
        return iFlag_loop, aFlowline_list

    aList_all = list()

    for i in range(0,nFlowline):
        pFlowline_in = aFlowline_in[i]
        aFlowline_list=list()
        aVertex_list=list()
        iFlag_loop , aFlowline_list = walk_cycle(0, 0, pFlowline_in, pFlowline_in.pVertex_start, aFlowline_list, aVertex_list)
        if iFlag_loop == 1:            
            aFlowline_list.insert(0, i)
            iFlag = find_list_in_list(aList_all,aFlowline_list )
            if iFlag == 1:
                pass
            else:
                aList_all.append(aFlowline_list)
                

    
    pDataset = pDriver_geojson.CreateDataSource(sFilename_output_in) 
    pSpatial_reference_gcs = osr.SpatialReference()  
    pSpatial_reference_gcs.ImportFromEPSG(4326)    # WGS84 lat/lon
    pSpatial_reference_gcs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

    pLayer = pDataset.CreateLayer('aod', pSpatial_reference_gcs, ogr.wkbPolygon)
    # Add one attribute
    pLayer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger64)) #long type for high resolution

    pLayerDefn = pLayer.GetLayerDefn()
    pFeature = ogr.Feature(pLayerDefn)    
    lCellID = 0
    dArea = 0.0

    #shapely method, but not complete
    for n in range(0,len(aList_all)):
        aFlowline_list = aList_all[n]
        nFlowline = len(aFlowline_list)
        aFlowline_polygon=list()
        for m in range(nFlowline):
            dummy_index = aFlowline_list[m]
            pFlowline = aFlowline_in[dummy_index]
            nVertex= pFlowline.nVertex
            aCoords_gcs = np.full((nVertex,2), -9999. ,dtype=float)
            for k in range(nVertex):
                aCoords_gcs[k,0] = pFlowline.aVertex[k].dLongitude_degree
                aCoords_gcs[k,1] = pFlowline.aVertex[k].dLatitude_degree
            
            aCoords_gcs= tuple(j for j in aCoords_gcs) #np.array(aCoords_gcs)    
            aFlowline_polygon.append(aCoords_gcs)   

        dummy = polygonize(aFlowline_polygon)
        aPolygon_out =  list(dummy) 
  
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
   
    

    

      

             

    
