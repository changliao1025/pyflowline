import os, sys
import numpy as np
from osgeo import ogr, osr, gdal
sys.setrecursionlimit(100000)
from tinyr import RTree
from pyearth.gis.spatialref.convert_between_degree_and_meter import meter_to_degree
from pyflowline.classes.vertex import pyvertex
from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_flowline
from pyflowline.formats.export_flowline import export_flowline_to_geojson
from pyflowline.algorithms.cython.kernel import calculate_distance_based_on_longitude_latitude

def calculate_distance_between_flowlines(pFlowline_a, pFlowline_b, dDistance_buffer):
    #this function calculates the distance between two flowleine using all the points, not the midpoint
    dDistance_minimal = 1.0E10
    #start with the first flowline
    aVertex_a = pFlowline_a.aVertex
    #use r tree
    index_a = RTree(max_cap=5, min_cap=2)
    for i in range(len(aVertex_a)):
        pVertex = aVertex_a[i]
        dLon_min = pVertex.dLongitude_degree
        dLat_min = pVertex.dLatitude_degree
        dLon_max = pVertex.dLongitude_degree
        dLat_max = pVertex.dLatitude_degree
        #add some buffer using the threshold
        dBuffer_degree = meter_to_degree(dDistance_buffer, (dLat_min + dLat_max) / 2.0)
        pBound = ( dLon_min - dBuffer_degree, dLat_min - dBuffer_degree,
                   dLon_max + dBuffer_degree, dLat_max + dBuffer_degree)
        index_a.insert(i, pBound)

    #now start with the second flowline
    aVertex_b = pFlowline_b.aVertex
    for i in range(len(aVertex_b)):
        pVertex = aVertex_b[i]
        #define a buffer using the radius
        dummy=dict()
        dummy['dLongitude_degree'] =pVertex.dLongitude_degree
        dummy['dLatitude_degree'] = pVertex.dLatitude_degree

        #create the bound
        dLon_min = pVertex.dLongitude_degree
        dLat_min = pVertex.dLatitude_degree
        dBuffer_degree = meter_to_degree(dDistance_buffer, (dLat_min + dLat_max) / 2.0)
        pBound = ( dLon_min - dBuffer_degree, dLat_min - dBuffer_degree,
                   dLon_max + dBuffer_degree, dLat_max + dBuffer_degree)

        aIntersect = list(index_a.search(pBound))
        if len(aIntersect) > 0:
            #check the distance
            for j in range(len(aIntersect)):
                #get the vertex
                pVertex_a = aVertex_a[aIntersect[j]]
                #calculate the distance
                dDistance = calculate_distance_based_on_longitude_latitude(pVertex.dLongitude_degree, pVertex.dLatitude_degree,
                                                                            pVertex_a.dLongitude_degree, pVertex_a.dLatitude_degree)
                if dDistance < dDistance_minimal:
                    dDistance_minimal = dDistance
        else:
            pass

    return dDistance_minimal


def simplify_hydrosheds_river_network(sFilename_flowline_hydroshed_in,
                       sFilename_flowline_hydroshed_out,
                       dDistance_tolerance_in,
                        dDrainage_area_threshold_in):
    ### Simplify hydroshed flowlines
    #check file exists
    if not os.path.isfile(sFilename_flowline_hydroshed_in):
        print('This input file does not exist: ', sFilename_flowline_hydroshed_in )
        return 0
    pDriver_geojson = ogr.GetDriverByName("GeoJSON")
    pDriver_shapefile = ogr.GetDriverByName("ESRI Shapefile")

    #check the file type using extension
    sFile_extension = os.path.splitext(sFilename_flowline_hydroshed_in)[1]
    if sFile_extension == '.geojson':
        pDataset_in = pDriver_geojson.Open(sFilename_flowline_hydroshed_in, gdal.GA_ReadOnly)
    else:
        pDataset_in = pDriver_shapefile.Open(sFilename_flowline_hydroshed_in, gdal.GA_ReadOnly)

    pLayer_shapefile = pDataset_in.GetLayer(0)
    pSpatialRef_shapefile = pLayer_shapefile.GetSpatialRef()

    pProjection_geojson = pSpatialRef_shapefile.ExportToWkt()

    pSpatial_reference_gcs = osr.SpatialReference()
    pSpatial_reference_gcs.ImportFromEPSG(4326)
    pSpatial_reference_gcs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

    #we dont care about whether it goes into the ocean or it is a endorheic basin
    #make a list for the flowlines that will be


    if os.path.exists(sFilename_flowline_hydroshed_out):
        os.remove(sFilename_flowline_hydroshed_out)

    aFlowline_hydroshed_outlet = []
    aFlowline_hydroshed_upstream_all = []
    aFlowlineID_outlet = []


    #get the number of features
    iNumber_of_features = pLayer_shapefile.GetFeatureCount()

    lFlowlineIndex = 0
    def convert_geometry_flowline(pGeometry_in, lFlowlineIndex, lID, lOutletID, lStream_order):
        aCoords = list()
        nPoint = pGeometry_in.GetPointCount()
        if nPoint < 2:
            print('This is an empty flowline')
            return None
        else:
            for i in range(0, nPoint):
                pt = pGeometry_in.GetPoint(i)
                aCoords.append( [ pt[0], pt[1]])

        dummy1= np.array(aCoords)
        pFlowline = convert_gcs_coordinates_to_flowline(dummy1)
        pFlowline.lFlowlineIndex = lFlowlineIndex
        pFlowline.lFlowlineID= lID
        pFlowline.iStream_segment = lOutletID
        pFlowline.iStream_order = lStream_order
        return pFlowline


    for i in range(0, iNumber_of_features):
        pFeature_shapefile = pLayer_shapefile.GetFeature(i)
        pGeometry_shapefile = pFeature_shapefile.GetGeometryRef()
        sGeometry_type = pGeometry_shapefile.GetGeometryName()
        lID = pFeature_shapefile.GetFieldAsInteger("HYRIV_ID")
        lOutletID = pFeature_shapefile.GetFieldAsInteger("MAIN_RIV")
        lStream_order = pFeature_shapefile.GetFieldAsInteger("ORD_STRA")
        #get the drainage area
        dDrainage_area = pFeature_shapefile.GetFieldAsDouble("UPLAND_SKM") * 1.0E6 #unit: km^2 to m^2
        #get the flag whether it flows into the ocean
        iFlag_edge = pFeature_shapefile.GetFieldAsInteger("NEXT_DOWN") #0 is ocean, non-0 is the id of next downstream flowline
        #get the flag whether it is an endorheic basin
        iFlag_endorheic = pFeature_shapefile.GetFieldAsInteger("ENDORHEIC") #0 = not part of an endorheic basin; 1 = part of an endorheic basin.

        if iFlag_edge == 0 and dDrainage_area > dDrainage_area_threshold_in: #river flow to ocean or inland sink
            if sGeometry_type == 'LINESTRING':
                pFlowline = convert_geometry_flowline(pGeometry_shapefile, lFlowlineIndex, lID, lOutletID, lStream_order)
                if pFlowline is not None:
                    #aFlowlineID_outlet.append(lOutletID)
                    pFlowline.dDrainage_area = dDrainage_area
                    aFlowline_hydroshed_outlet.append(pFlowline)
                    lFlowlineIndex = lFlowlineIndex + 1
                    pass
            else:
                if sGeometry_type == 'MULTILINESTRING':
                    #loop through all the lines
                    nLine = pGeometry_shapefile.GetGeometryCount()
                    for j in range(0, nLine):
                        pGeometry_line = pGeometry_shapefile.GetGeometryRef(j)
                        pFlowline = convert_geometry_flowline(pGeometry_line, lFlowlineIndex, lID, lOutletID, lStream_order)
                        if pFlowline is not None:
                            #aFlowlineID_outlet.append(lOutletID)
                            pFlowline.dDrainage_area = dDrainage_area
                            aFlowline_hydroshed_outlet.append(pFlowline)
                            lFlowlineIndex = lFlowlineIndex + 1
            pass
        else:
            #not endorheic basin, but it has a large drainage area
            if dDrainage_area > dDrainage_area_threshold_in: #not next to ocean but has a large drainage area
                if sGeometry_type == 'LINESTRING':
                    pFlowline = convert_geometry_flowline(pGeometry_shapefile, lFlowlineIndex, lID, lOutletID, lStream_order)
                    if pFlowline is not None:
                        pFlowline.lFlowlineID_downstream = iFlag_edge #this is the downstream flowline
                        pFlowline.dDrainage_area = dDrainage_area
                        pFlowline.iStream_segment = lOutletID
                        pFlowline.iStream_order = lStream_order
                        aFlowline_hydroshed_upstream_all.append(pFlowline)
                        aFlowlineID_outlet.append(lOutletID)
                        lFlowlineIndex = lFlowlineIndex + 1
                else:
                    if sGeometry_type == 'MULTILINESTRING':
                        #loop through all the lines
                        nLine = pGeometry_shapefile.GetGeometryCount()
                        for j in range(0, nLine):
                            pGeometry_line = pGeometry_shapefile.GetGeometryRef(j)
                            pFlowline = convert_geometry_flowline(pGeometry_line, lFlowlineIndex, lID, lOutletID, lStream_order)
                            if pFlowline is not None:
                                pFlowline.lFlowlineID_downstream = iFlag_edge
                                pFlowline.dDrainage_area = dDrainage_area
                                pFlowline.iStream_segment = lOutletID
                                pFlowline.iStream_order = lStream_order
                                aFlowline_hydroshed_upstream_all.append(pFlowline)
                                aFlowlineID_outlet.append(lOutletID)
                                lFlowlineIndex = lFlowlineIndex + 1
                pass
            else:
                #small drainage area, we dont need them
                pass

    #step 1, filter outlet, if two outlets are too close, we need to remove smaller ones
    aFlowline_hydroshed_outlet_simplified=list()
    nFlowline_outlet = len(aFlowline_hydroshed_outlet)
    index_outlet = RTree(max_cap=5, min_cap=2)
    for i in range(nFlowline_outlet):
        pBound= aFlowline_hydroshed_outlet[i].pBound
        #add some buffer using the threshold
        dLon_min, dLat_min, dLon_max, dLat_max = pBound
        dBuffer_degree = meter_to_degree(dDistance_tolerance_in, (dLat_min + dLat_max) / 2.0)
        pBound = ( dLon_min - dBuffer_degree, dLat_min - dBuffer_degree,
                   dLon_max + dBuffer_degree, dLat_max + dBuffer_degree)
        index_outlet.insert(i, pBound)  #
        pass

    #use intersect to find the outlet flowlines that are too close
    for i in range(nFlowline_outlet):
        pBound= aFlowline_hydroshed_outlet[i].pBound
        dLon_min, dLat_min, dLon_max, dLat_max = pBound
        dBuffer_degree = meter_to_degree(dDistance_tolerance_in, (dLat_min + dLat_max) / 2.0)
        pBound = ( dLon_min - dBuffer_degree, dLat_min - dBuffer_degree,
                   dLon_max + dBuffer_degree, dLat_max + dBuffer_degree)
        aIntersect = list(index_outlet.search(pBound))
        lon1 = 0.5 * (aFlowline_hydroshed_outlet[i].pVertex_start.dLongitude_degree + aFlowline_hydroshed_outlet[i].pVertex_end.dLongitude_degree)
        lat1 = 0.5 * (aFlowline_hydroshed_outlet[i].pVertex_start.dLatitude_degree + aFlowline_hydroshed_outlet[i].pVertex_end.dLatitude_degree)
        for j in range(len(aIntersect)):
            if i != aIntersect[j] :
                lon2 = 0.5 * (aFlowline_hydroshed_outlet[aIntersect[j]].pVertex_start.dLongitude_degree + aFlowline_hydroshed_outlet[aIntersect[j]].pVertex_end.dLongitude_degree)
                lat2 = 0.5 * (aFlowline_hydroshed_outlet[aIntersect[j]].pVertex_start.dLatitude_degree + aFlowline_hydroshed_outlet[aIntersect[j]].pVertex_end.dLatitude_degree)
                dDistance = calculate_distance_based_on_longitude_latitude(lon1, lat1, lon2, lat2)
                dDrainage_area_a = aFlowline_hydroshed_outlet[i].dDrainage_area
                dDrainage_area_b = aFlowline_hydroshed_outlet[aIntersect[j]].dDrainage_area
                if dDistance < dDistance_tolerance_in:
                    #we need remove one of them
                    if dDrainage_area_a > dDrainage_area_b :
                        aFlowline_hydroshed_outlet[i].iFlag_keep = 1
                        aFlowline_hydroshed_outlet[aIntersect[j]].iFlag_keep = 0
                    else:
                        aFlowline_hydroshed_outlet[i].iFlag_keep = 0
                        aFlowline_hydroshed_outlet[aIntersect[j]].iFlag_keep = 1

    #collect the outlet flowlines
    for i in range(nFlowline_outlet):
        if aFlowline_hydroshed_outlet[i].iFlag_keep == 1:
            aFlowline_hydroshed_outlet_simplified.append(aFlowline_hydroshed_outlet[i])

    print('Number of valid outlet flowlines in the hydroshed: ', len(aFlowline_hydroshed_outlet_simplified))
    sFilename_flowline_hydroshed_tmp = sFilename_flowline_hydroshed_out.replace('.geojson', '_outlet.geojson')
    export_flowline_to_geojson(aFlowline_hydroshed_outlet_simplified, sFilename_flowline_hydroshed_tmp)

    aFlowlineID_downslope = []
    aFlowline_upstream = list()
    lFlowlineIndex = 0
    #flatten the array
    aFlowlineID_outlet = np.array(aFlowlineID_outlet).flatten()

    #step 2, upstream flowlines, if its most downstream flowline is in the outlet flowlines, we need to consider it
    for pFlowline in aFlowline_hydroshed_outlet_simplified:
        lFlowlineID = pFlowline.lFlowlineID
        dummy_index = np.where(aFlowlineID_outlet == lFlowlineID)
        #add all of them
        for i in range(len(dummy_index[0])):
            pFlowline_up = aFlowline_hydroshed_upstream_all[dummy_index[0][i]]
            pFlowline_up.lFlowlineIndex = lFlowlineIndex
            aFlowline_upstream.append(pFlowline_up)
            lFlowlineIndex = lFlowlineIndex + 1
            aFlowlineID_downslope.append(pFlowline_up.lFlowlineID_downstream)
            pass

    aFlowlineID_downslope = np.array(aFlowlineID_downslope)
    #we do not smplify the outlet flowlines, only the upstream flowlines
    #build the rtree
    index_reach = RTree(max_cap=5, min_cap=2)
    nFlowline = len(aFlowline_upstream)
    aFlag_process=np.full(nFlowline, 0, dtype =int)

    for i in range(nFlowline):
        pBound= aFlowline_upstream[i].pBound
        #add a little buffer zone around the bound

        dLon_min, dLat_min, dLon_max, dLat_max = pBound
        #convert the distance to degree
        dBuffer_degree = meter_to_degree(dDistance_tolerance_in, (dLat_min + dLat_max) / 2.0)
        pBound = ( dLon_min - dBuffer_degree, dLat_min - dBuffer_degree,
                   dLon_max + dBuffer_degree, dLat_max + dBuffer_degree)
        index_reach.insert(i, pBound)  #
        pass

    #save the flowlines
    aFlowline_before_distance_operation = aFlowline_hydroshed_outlet + aFlowline_upstream
    sFilename_flowline_hydroshed_tmp = sFilename_flowline_hydroshed_out.replace('.geojson', '_all.geojson')
    aAttribute_field = ['lineid', 'downstream_id', 'drainage_area']
    aAttribute_dtype= ['int', 'int', 'float']
    aAttribute_data=list()
    aAttribute_data.append( np.array([pFlowline.lFlowlineID for pFlowline in aFlowline_before_distance_operation]) )
    aAttribute_data.append( np.array([pFlowline.lFlowlineID_downstream for pFlowline in aFlowline_before_distance_operation]) )
    aAttribute_data.append( np.array([pFlowline.dDrainage_area for pFlowline in aFlowline_before_distance_operation]) )
    export_flowline_to_geojson(aFlowline_before_distance_operation, sFilename_flowline_hydroshed_tmp,
                                aAttribute_field=aAttribute_field,
    aAttribute_data=aAttribute_data,
    aAttribute_dtype=aAttribute_dtype)
    print('Number of all flowlines in the hydroshed: ', len(aFlowline_upstream))

    def find_upstream_flowline(lFlowlineIndex_in):
        lFlowlineID = aFlowline_upstream[lFlowlineIndex_in].lFlowlineID
        dummy_index = np.where(aFlowlineID_downslope == lFlowlineID)
        aUpstream = dummy_index[0]
        nUpstream = len(aUpstream)
        return nUpstream, aUpstream


    def mark_upstream(lFlowlineIndex_in):
        print("mark_upstream", aFlowline_upstream[lFlowlineIndex_in].lFlowlineID)
        if(aFlowline_upstream[lFlowlineIndex_in].iStream_order==1):
            pass
        else:
            nUpstream, aUpstream = find_upstream_flowline(lFlowlineIndex_in)
            if nUpstream > 0:
                for j in range(nUpstream):
                    pFlowline = aFlowline_upstream[ aUpstream[j] ]
                    pFlowline.iFlag_keep = 0
                    mark_upstream(  pFlowline.lFlowlineIndex )

    def is_downstream(lFlowlineIndex__down, lFlowlineIndex__up):
        lFlowlineID_a = aFlowline_upstream[lFlowlineIndex__down].lFlowlineID
        lFlowlineID_b = aFlowline_upstream[lFlowlineIndex__up].lFlowlineID
        if lFlowlineID_a == lFlowlineID_b:
            return 1
        else:
            nUpstream, aUpstream = find_upstream_flowline(lFlowlineIndex__down)
            if nUpstream > 0:
                for j in range(nUpstream):
                    if is_downstream(aUpstream[j], lFlowlineIndex__up) == 1:
                        return 1
                    pass
                pass
            return 0


    #use intersect to find the upstream flowlines
    for i in range(nFlowline):
        if aFlag_process[i] == 1:
            continue

        print('Processing flowline: ', i, aFlowline_upstream[i].lFlowlineID)
        pBound= aFlowline_upstream[i].pBound
        dLon_min, dLat_min, dLon_max, dLat_max = pBound
        #convert the distance to degree
        dBuffer_degree = meter_to_degree(dDistance_tolerance_in, (dLat_min + dLat_max) / 2.0)
        pBound = ( dLon_min - dBuffer_degree, dLat_min - dBuffer_degree,
                   dLon_max + dBuffer_degree, dLat_max + dBuffer_degree)
        aIntersect = list(index_reach.search(pBound))
        lon1 = 0.5 * (aFlowline_upstream[i].pVertex_start.dLongitude_degree + aFlowline_upstream[i].pVertex_end.dLongitude_degree)
        lat1 = 0.5 * (aFlowline_upstream[i].pVertex_start.dLatitude_degree + aFlowline_upstream[i].pVertex_end.dLatitude_degree)
        for j in range(len(aIntersect)):
            if i != aIntersect[j] :
                if aFlowline_upstream[aIntersect[j]].iFlag_keep == 0:
                    continue
                #calculate the distance between two flowlines
                #check whether they share a end point
                lon2 = 0.5 * (aFlowline_upstream[aIntersect[j]].pVertex_start.dLongitude_degree + aFlowline_upstream[aIntersect[j]].pVertex_end.dLongitude_degree)
                lat2 = 0.5 * (aFlowline_upstream[aIntersect[j]].pVertex_start.dLatitude_degree + aFlowline_upstream[aIntersect[j]].pVertex_end.dLatitude_degree)
                if aFlowline_upstream[i].pVertex_end == aFlowline_upstream[aIntersect[j]].pVertex_end:
                    print(lon1, lat1, lon2, lat2)
                    pass
                else:
                    dDistance = calculate_distance_between_flowlines(aFlowline_upstream[i], aFlowline_upstream[aIntersect[j]], dDistance_tolerance_in * 2)
                    dDistance2 = calculate_distance_based_on_longitude_latitude(lon1, lat1, lon2, lat2)
                    iStream_segment_a= aFlowline_upstream[i].iStream_segment
                    iStream_segment_b= aFlowline_upstream[aIntersect[j]].iStream_segment
                    iStream_order_a = aFlowline_upstream[i].iStream_order
                    iStream_order_b = aFlowline_upstream[aIntersect[j]].iStream_order
                    dDrainage_area_a = aFlowline_upstream[i].dDrainage_area
                    dDrainage_area_b = aFlowline_upstream[aIntersect[j]].dDrainage_area
                    if iStream_segment_a == iStream_segment_b: #same basin
                        if is_downstream(aIntersect[j], i) == 1 or is_downstream(i, aIntersect[j]) == 1:
                            #they are on the same channel
                            continue
                        else: #same river basin but different subbasin
                            if dDistance < dDistance_tolerance_in:
                                if dDrainage_area_a > dDrainage_area_b:
                                    aFlowline_upstream[aIntersect[j]].iFlag_keep = 0
                                    #mark all upstream flowlines as 0
                                    mark_upstream( aIntersect[j] )
                                else:
                                    aFlowline_upstream[i].iFlag_keep = 0
                                    mark_upstream( i )
                                pass

                    else: #different river basin
                        if dDistance < dDistance_tolerance_in:
                            #we need remove one of them
                            if dDrainage_area_a > dDrainage_area_b :
                                aFlowline_upstream[aIntersect[j]].iFlag_keep = 0
                                #mark all upstream flowlines as 0
                                mark_upstream( aIntersect[j] )
                            else:
                                aFlowline_upstream[i].iFlag_keep = 0
                                mark_upstream( i )

                pass
            pass

    aFlowline_upstream_simplified = list()
    for i in range(nFlowline):
        if aFlowline_upstream[i].iFlag_keep == 1:
            aFlowline_upstream_simplified.append(aFlowline_upstream[i])

    #save the flowlines
    aFlowline_after_distance_operation = aFlowline_hydroshed_outlet_simplified + aFlowline_upstream_simplified
    export_flowline_to_geojson(aFlowline_after_distance_operation, sFilename_flowline_hydroshed_out)
    #close the file
    pDataset_in = pLayer_shapefile = pFeature_shapefile = None
    print('Number of flowlines in the hydroshed: ', len(aFlowline_upstream_simplified))
    return

