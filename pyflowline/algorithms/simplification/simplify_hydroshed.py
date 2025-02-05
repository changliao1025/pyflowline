import os, sys
import numpy as np
from osgeo import ogr, osr, gdal
sys.setrecursionlimit(100000)
from tinyr import RTree
from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_flowline
from pyflowline.formats.export_flowline import export_flowline_to_geojson
from pyflowline.algorithms.cython.kernel import calculate_distance_based_on_longitude_latitude

def simplify_hydroshed(sFilename_flowline_hydroshed_in, sFilename_flowline_hydroshed_out, dDistance_tolerance_in, dDrainage_area_threshold_in):
    ### Simplify hydroshed flowlines
    #check file exists
    if not os.path.isfile(sFilename_flowline_hydroshed_in):
        print('This input file does not exist: ', sFilename_flowline_hydroshed_in )
        return 0

    pDriver_shapefile = ogr.GetDriverByName("ESRI Shapefile")
    pDataset_shapefile = pDriver_shapefile.Open(sFilename_flowline_hydroshed_in, gdal.GA_ReadOnly)
    pLayer_shapefile = pDataset_shapefile.GetLayer(0)
    pSpatialRef_shapefile = pLayer_shapefile.GetSpatialRef()

    pProjection_geojson = pSpatialRef_shapefile.ExportToWkt()

    pSpatial_reference_gcs = osr.SpatialReference()
    pSpatial_reference_gcs.ImportFromEPSG(4326)
    pSpatial_reference_gcs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

    #we dont care about whether it goes into the ocean or it is a endorheic basin
    #make a list for the flowlines that will be

    pDriver_geojson = ogr.GetDriverByName("GeoJSON")
    if os.path.exists(sFilename_flowline_hydroshed_out):
        os.remove(sFilename_flowline_hydroshed_out)

    aFlowline_hydroshed_outlet = []
    aFlowline_hydroshed_upstream_all = []
    aFlowlineID_outlet = []


    #get the number of features
    iNumber_of_features = pLayer_shapefile.GetFeatureCount()

    lFlowlineIndex = 0
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
        iFlag_ocean = pFeature_shapefile.GetFieldAsInteger("NEXT_DOWN") #0 is ocean, non-0 is the id of next downstream flowline
        #get the flag whether it is an endorheic basin
        iFlag_endorheic = pFeature_shapefile.GetFieldAsInteger("ENDORHEIC") #0 = not part of an endorheic basin; 1 = part of an endorheic basin.
        if iFlag_ocean == 0 and dDrainage_area > dDrainage_area_threshold_in:
            #keep this
            aCoords = list()
            for i in range(0, pGeometry_shapefile.GetPointCount()):
                pt = pGeometry_shapefile.GetPoint(i)
                aCoords.append( [ pt[0], pt[1]])
            dummy1= np.array(aCoords)
            pFlowline = convert_gcs_coordinates_to_flowline(dummy1)
            pFlowline.lFlowlineIndex = lFlowlineIndex
            pFlowline.lFlowlineID= lID
            pFlowline.iStream_segment = lOutletID
            pFlowline.iStream_order = lStream_order
            aFlowline_hydroshed_outlet.append(pFlowline)
            lFlowlineIndex = lFlowlineIndex + 1

            pass

        else:
            if iFlag_endorheic == 1 and dDrainage_area > dDrainage_area_threshold_in:
                #keep this
                #print('This is an endorheic basin')
                #keep this
                aCoords = list()
                for i in range(0, pGeometry_shapefile.GetPointCount()):
                    pt = pGeometry_shapefile.GetPoint(i)
                    aCoords.append( [ pt[0], pt[1]])
                dummy1= np.array(aCoords)
                pFlowline = convert_gcs_coordinates_to_flowline(dummy1)
                pFlowline.lFlowlineIndex = lFlowlineIndex
                pFlowline.lFlowlineID= lID
                pFlowline.iStream_segment = lOutletID
                pFlowline.iStream_order = lStream_order
                aFlowline_hydroshed_outlet.append(pFlowline)
                lFlowlineIndex = lFlowlineIndex + 1

                pass
            else:
                if dDrainage_area > dDrainage_area_threshold_in: #not next to ocean but has a large drainage area
                    aCoords = list()
                    for i in range(0, pGeometry_shapefile.GetPointCount()):
                        pt = pGeometry_shapefile.GetPoint(i)
                        aCoords.append( [ pt[0], pt[1]])
                    dummy1= np.array(aCoords)
                    pFlowline = convert_gcs_coordinates_to_flowline(dummy1)
                    pFlowline.lFlowlineID= lID
                    pFlowline.lFlowlineID_downstream = iFlag_ocean #this is the downstream flowline
                    pFlowline.dDrainage_area = dDrainage_area
                    pFlowline.iStream_segment = lOutletID
                    pFlowline.iStream_order = lStream_order
                    aFlowline_hydroshed_upstream_all.append(pFlowline)

                    aFlowlineID_outlet.append(lOutletID)
                    pass
                else:
                    #small drainage area, we dont need them
                    pass

    #filter outlet first
    aFlowline_hydroshed_outlet_simplified=list()
    nFlowline_outlet = len(aFlowline_hydroshed_outlet)
    index_outlet = RTree(max_cap=5, min_cap=2)
    for i in range(nFlowline_outlet):
        pBound= aFlowline_hydroshed_outlet[i].pBound
        index_outlet.insert(i, pBound)  #
        pass
    #use intersect to find the upstream flowlines
    for i in range(nFlowline_outlet):
        pBound= aFlowline_hydroshed_outlet[i].pBound
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

    #collect the flowlines
    for i in range(nFlowline_outlet):
        if aFlowline_hydroshed_outlet[i].iFlag_keep == 1:
            aFlowline_hydroshed_outlet_simplified.append(aFlowline_hydroshed_outlet[i])

    print('Number of flowlines in the hydroshed: ', len(aFlowline_hydroshed_outlet_simplified))

    aFlowlineID_downslope = []
    aFlowline_upstream = list()
    lFlowlineIndex=0
    aFlowlineID_outlet = np.array(aFlowlineID_outlet)
    #for pFlowline in aFlowline_hydroshed_outlet:
    for pFlowline in aFlowline_hydroshed_outlet_simplified:
        lFlowlineID = pFlowline.lFlowlineID
        dummy_index = np.where(aFlowlineID_outlet == lFlowlineID)
        #add all of them
        for i in range(len(dummy_index[0])):
            pFlowline_up = aFlowline_hydroshed_upstream_all[dummy_index[0][i]]
            pFlowline_up = pFlowline_up.split_by_length(dDistance_tolerance_in)
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
    print('Number of flowlines in the hydroshed: ', len(aFlowline_upstream))

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
        aIntersect = list(index_reach.search(pBound))
        lon1 = 0.5 * (aFlowline_upstream[i].pVertex_start.dLongitude_degree + aFlowline_upstream[i].pVertex_end.dLongitude_degree)
        lat1 = 0.5 * (aFlowline_upstream[i].pVertex_start.dLatitude_degree + aFlowline_upstream[i].pVertex_end.dLatitude_degree)
        for j in range(len(aIntersect)):
            if i != aIntersect[j] :
                if aFlowline_upstream[aIntersect[j]].iFlag_keep == 0:
                    continue
                #calculate the distance
                #check whether they share a end point

                if aFlowline_upstream[i].pVertex_end == aFlowline_upstream[aIntersect[j]].pVertex_end:
                    pass
                else:
                    lon2 = 0.5 * (aFlowline_upstream[aIntersect[j]].pVertex_start.dLongitude_degree + aFlowline_upstream[aIntersect[j]].pVertex_end.dLongitude_degree)
                    lat2 = 0.5 * (aFlowline_upstream[aIntersect[j]].pVertex_start.dLatitude_degree + aFlowline_upstream[aIntersect[j]].pVertex_end.dLatitude_degree)
                    dDistance = calculate_distance_based_on_longitude_latitude(lon1, lat1, lon2, lat2)
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
    aFlowline_after_distance_operation = aFlowline_hydroshed_outlet + aFlowline_upstream_simplified
    export_flowline_to_geojson(aFlowline_after_distance_operation, sFilename_flowline_hydroshed_out)
    #close the file
    pDataset_shapefile = pLayer_shapefile = pFeature_shapefile = None
    print('Number of flowlines in the hydroshed: ', len(aFlowline_upstream_simplified))
    return

