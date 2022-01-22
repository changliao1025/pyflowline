import os
import json
import numpy as np
from osgeo import ogr, osr, gdal, gdalconst
from shapely.geometry import Point, LineString, MultiLineString
from shapely.wkt import loads

from shapely.ops import polygonize, polygonize_full
from pyflowline.algorithms.auxiliary.gdal_functions import  calculate_angle_betwen_vertex
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
    aFlowline = list()


    #rebuild
    for n, v in enumerate(aFlowline_in):
        v.lFlowlineID = n
    for n, v in enumerate(aVertex_all_in):
        v.lVertexID = n
    for n, v in enumerate(aVertex_all_in):
    
        for i, feat in enumerate(aFlowline_in):
            if v == feat.pVertex_start:
                feat.pVertex_start.lVertexID = n

            elif v == feat.pVertex_end:
                feat.pVertex_end.lVertexID = n

    contact_adj_dict = {}

    for i, c in enumerate(aFlowline_in):
        #contact_adj_dict[i] = {'start': [],
        #                       'stop': []}

        for j, cc in enumerate(aFlowline_in):
            if j != i:
                if c.pVertex_start in (cc.pVertex_start, cc.pVertex_end):
                    c.aFlowlineID_start.append(cc.lFlowlineID)
                    #contact_adj_dict[i]['start'].append(j)
                elif c.pVertex_end in (cc.pVertex_start, cc.pVertex_end):                    
                 #c['properties']['stop'] in (
                 #   cc['properties']['start'], cc['properties']['stop']):
                    #contact_adj_dict[i]['stop'].append(j)
                    c.aFlowlineID_end.append(cc.lFlowlineID)

    def angle_between(c1, c2, c3):
    
        angle = (np.arctan2(c3[1] - c1[1], c3[0] - c1[0]) -
                 np.arctan2(c2[1] - c1[1], c2[0] - c1[0]))

        if angle < 0:
            angle +=  2 * np.pi

        return angle


    


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
            for i in pFlowline_in.aFlowlineID_end:
                pFlowline_dummy = aFlowline_in[i-1]
                if (pFlowline_dummy.pVertex_start == pVertex_stop):
                    pVertex_start_dummy = pFlowline_dummy.pVertex_start
                else:
                    pVertex_start_dummy = pFlowline_dummy.pVertex_end
                
                #calculate angle
                angle_dummy = calculate_angle_betwen_vertex()


        else:
            #rever direction
            pVertex_stop = pFlowline_in.pVertex_start

        


    def walk_cycle_right( pFlowline_in):
        pVertex_start = pFlowline_in.pVertex_start
        pVertex_end = pFlowline_in.pVertex_end
        pFlowline_next, pVertex_stop  = get_right_branch(0, pFlowline_in)

        if (pVertex_stop == pVertex_start):
            #a loop is finished

            pass
        else:
            pFlowline_in = pFlowline_next
            walk_cycle_right( pFlowline_in)

    
        return

    for i in range(nFlowline):
        pFlowline_in = aFlowline_in[i]
        walk_cycle_right(pFlowline_in)



        
    
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
   
    

    

      

             

    
