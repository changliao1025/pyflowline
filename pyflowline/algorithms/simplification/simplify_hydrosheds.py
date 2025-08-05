import os, sys
sys.setrecursionlimit(100000)
import numpy as np
from osgeo import ogr, osr, gdal
sys.setrecursionlimit(100000)
from tinyr import RTree
from pyearth.gis.spatialref.convert_between_degree_and_meter import meter_to_degree
from pyearth.toolbox.geometry.create_gcs_buffer_zone import create_polyline_buffer_zone, create_point_buffer_zone
from pyflowline.classes.vertex import pyvertex
from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_flowline
from pyflowline.formats.export_flowline import export_flowline_to_geojson
from pyflowline.algorithms.cython.kernel import calculate_distance_based_on_longitude_latitude, calculate_distance_based_on_longitude_latitude_numpy
from pyflowline.algorithms.index.define_stream_order import define_stream_order, update_head_water_stream_order
from pyflowline.algorithms.merge.merge_flowline import merge_flowline
from pyflowline.algorithms.index.define_stream_segment_index import define_stream_segment_index
from pyflowline.algorithms.split.find_flowline_confluence import find_flowline_confluence
from pyflowline.algorithms.index.define_stream_topology import define_stream_topology
from pyflowline.classes.confluence import pyconfluence
from pyflowline.algorithms.split.find_flowline_vertex import find_flowline_vertex
from pyflowline.algorithms.split.split_flowline import split_flowline

from pyflowline.configuration.config_manager import create_template_configuration_file
from pyflowline.configuration.change_json_key_value import change_json_key_value


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

def precompute_flowline_geometries(aFlowlines, dDistance_tolerance):
    """
    Precomputes and caches flowline geometries to avoid repeated calculations.
    Args:
        aFlowlines: List of flowline objects
        dDistance_tolerance: Distance tolerance for buffer creation
    Returns:
        Tuple of (bounds_cache, buffer_cache) dictionaries keyed by flowline ID
    """
    print(f"Precomputing geometries for {len(aFlowlines)} flowlines...")
    bounds_cache = {}
    #buffer_cache = {}
    for i, pFlowline in enumerate(aFlowlines):
        wkt = pFlowline.wkt
        wkt_buffer = create_polyline_buffer_zone(wkt, dDistance_tolerance)
        pBuffer = ogr.CreateGeometryFromWkt(wkt_buffer)
        pBound0 = pBuffer.GetEnvelope()
        dLon_min, dLon_max, dLat_min, dLat_max = pBound0
        pBound = (dLon_min, dLat_min, dLon_max, dLat_max)
        bounds_cache[pFlowline.lFlowlineID] = pBound
        #buffer_cache[pFlowline.lFlowlineID] = wkt_buffer
    return bounds_cache #, buffer_cache

def precompute_flowline_geometries_by_segment(aFlowlines, dDistance_tolerance):
    """
    Precomputes and caches flowline geometries to avoid repeated calculations.
    Args:
        aFlowlines: List of flowline objects
        dDistance_tolerance: Distance tolerance for buffer creation
    Returns:
        Tuple of (bounds_cache, buffer_cache) dictionaries keyed by flowline ID
    """
    print(f"Precomputing geometries for {len(aFlowlines)} flowlines...")
    bounds_cache = {}
    #buffer_cache = {}
    for i, pFlowline in enumerate(aFlowlines):
        wkt = pFlowline.wkt
        #if pFlowline.iStream_segment == 580:
        #    print('debugging flowline: ', pFlowline.lFlowlineID, pFlowline.iStream_segment)
        wkt_buffer = create_polyline_buffer_zone(wkt, dDistance_tolerance)
        pBuffer = ogr.CreateGeometryFromWkt(wkt_buffer)
        pBound0 = pBuffer.GetEnvelope()
        dLon_min, dLon_max, dLat_min, dLat_max = pBound0
        pBound = (dLon_min, dLat_min, dLon_max, dLat_max)
        bounds_cache[pFlowline.iStream_segment] = pBound
        #buffer_cache[pFlowline.lFlowlineID] = wkt_buffer
    return bounds_cache #, buffer_cache

def basin_build_confluence( aFlowline_basin_in, aVertex_confluence_in):
    #this can only be calculated for confluence
    # Create a dictionary to map each vertex to its upstream and downstream flowlines
    vertex_to_flowlines = {}
    for pFlowline in aFlowline_basin_in:
        pVertex_start = pFlowline.pVertex_start
        pVertex_end = pFlowline.pVertex_end
        if pVertex_end not in vertex_to_flowlines:
            vertex_to_flowlines[pVertex_end] = {'upstream': [], 'downstream': None}
        vertex_to_flowlines[pVertex_end]['upstream'].append(pFlowline)
        if pVertex_start not in vertex_to_flowlines:
            vertex_to_flowlines[pVertex_start] = {'upstream': [], 'downstream': None}
        vertex_to_flowlines[pVertex_start]['downstream'] = pFlowline
    # Build the confluence for each vertex
    aConfluence_basin = []
    for pVertex in aVertex_confluence_in:
        aFlowline_upstream = vertex_to_flowlines[pVertex]['upstream']
        pFlowline_downstream = vertex_to_flowlines[pVertex]['downstream']
        pConfluence = pyconfluence(pVertex, aFlowline_upstream, pFlowline_downstream)
        aConfluence_basin.append(pConfluence)
    return aConfluence_basin

def simplify_hydrosheds_river_network(sFilename_flowline_hydroshed_in,
                       sFilename_flowline_hydroshed_out,
                       dDistance_tolerance_in,
                        dDrainage_area_threshold_in,
                        iFlag_pyflowline_configuration_in=1,
                        nOutlet_largest= 10):

    dDrainage_area_threshold_ratio = 0.05
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


    if os.path.exists(sFilename_flowline_hydroshed_out):
        os.remove(sFilename_flowline_hydroshed_out)

    aFlowline_hydroshed_outlet = []
    aFlowline_hydroshed_upstream_all = []
    aFlowlineID_outlet = []

    lFlowlineIndex = 0

    #first we will find all the flowlines that flow into the ocean or inland sink
    for i, pFeature_shapefile in enumerate(pLayer_shapefile):
        fid = pFeature_shapefile.GetFID()
        #pFeature_shapefile = pLayer_shapefile.GetFeature(i)
        pGeometry_shapefile = pFeature_shapefile.GetGeometryRef()
        sGeometry_type = pGeometry_shapefile.GetGeometryName()
        lID = pFeature_shapefile.GetFieldAsInteger("HYRIV_ID")
        #if lID == 70784289:
        #    print('debugging flowline: ', lID)

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
                    pFlowline.iFlag_endorheic = iFlag_endorheic
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
                            pFlowline.iFlag_endorheic = iFlag_endorheic
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
                        pFlowline.iFlag_endorheic= 0
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
                                pFlowline.iFlag_endorheic= 0
                                pFlowline.iStream_segment = lOutletID
                                pFlowline.iStream_order = lStream_order
                                aFlowline_hydroshed_upstream_all.append(pFlowline)
                                aFlowlineID_outlet.append(lOutletID)
                                lFlowlineIndex = lFlowlineIndex + 1
                pass
            else:
                #small drainage area, we dont need them
                pass


    print("Precomputing outlet flowline geometries...")
    outlet_bounds_cache = precompute_flowline_geometries(aFlowline_hydroshed_outlet, dDistance_tolerance_in)
    # step 1, filter outlet, if two outlets are too close, we need to remove smaller ones

    nFlowline_outlet = len(aFlowline_hydroshed_outlet)
    index_outlet = RTree(max_cap=5, min_cap=2)

    # Use the precomputed bounds to build the RTree
    for i in range(nFlowline_outlet):
        pflowline = aFlowline_hydroshed_outlet[i]
        pflowline.iFlag_keep = 1
        pBound = outlet_bounds_cache[pflowline.lFlowlineID]
        index_outlet.insert(i, pBound)

    #use intersect to find the outlet flowlines that are too close
    for i in range(nFlowline_outlet):
        # Use cached bounds directly
        pBound = outlet_bounds_cache[aFlowline_hydroshed_outlet[i].lFlowlineID]
        aIntersect = list(index_outlet.search(pBound))
        #lon1 = 0.5 * (aFlowline_hydroshed_outlet[i].pVertex_start.dLongitude_degree + aFlowline_hydroshed_outlet[i].pVertex_end.dLongitude_degree)
        #lat1 = 0.5 * (aFlowline_hydroshed_outlet[i].pVertex_start.dLatitude_degree + aFlowline_hydroshed_outlet[i].pVertex_end.dLatitude_degree)
        lon1 = aFlowline_hydroshed_outlet[i].pVertex_end.dLongitude_degree
        lat1 = aFlowline_hydroshed_outlet[i].pVertex_end.dLatitude_degree
        for j in range(len(aIntersect)):
            if i != aIntersect[j] :
                #lon2 = 0.5 * (aFlowline_hydroshed_outlet[aIntersect[j]].pVertex_start.dLongitude_degree + aFlowline_hydroshed_outlet[aIntersect[j]].pVertex_end.dLongitude_degree)
                #lat2 = 0.5 * (aFlowline_hydroshed_outlet[aIntersect[j]].pVertex_start.dLatitude_degree + aFlowline_hydroshed_outlet[aIntersect[j]].pVertex_end.dLatitude_degree)
                lon2 = aFlowline_hydroshed_outlet[aIntersect[j]].pVertex_end.dLongitude_degree
                lat2 = aFlowline_hydroshed_outlet[aIntersect[j]].pVertex_end.dLatitude_degree
                dDistance = calculate_distance_based_on_longitude_latitude(lon1, lat1, lon2, lat2)
                dDrainage_area_a = aFlowline_hydroshed_outlet[i].dDrainage_area
                dDrainage_area_b = aFlowline_hydroshed_outlet[aIntersect[j]].dDrainage_area
                if aFlowline_hydroshed_outlet[i].iFlag_endorheic == 1 and aFlowline_hydroshed_outlet[aIntersect[j]].iFlag_endorheic == 1:
                    #both are endorheic, we need to keep both
                    aFlowline_hydroshed_outlet[i].iFlag_keep = 1
                    aFlowline_hydroshed_outlet[aIntersect[j]].iFlag_keep = 1
                else:
                    if dDistance < dDistance_tolerance_in:
                        #we need remove one of them
                        if dDrainage_area_a > dDrainage_area_b :
                            aFlowline_hydroshed_outlet[i].iFlag_keep = 1
                            aFlowline_hydroshed_outlet[aIntersect[j]].iFlag_keep = 0
                        else:
                            aFlowline_hydroshed_outlet[i].iFlag_keep = 0
                            aFlowline_hydroshed_outlet[aIntersect[j]].iFlag_keep = 1

    #collect the outlet flowlines
    aFlowline_hydroshed_outlet_simplified = list()
    for i in range(nFlowline_outlet):
        if aFlowline_hydroshed_outlet[i].iFlag_keep == 1:
            aFlowline_hydroshed_outlet_simplified.append(aFlowline_hydroshed_outlet[i])

    print('Number of valid outlet flowlines in the hydroshed: ', len(aFlowline_hydroshed_outlet_simplified))
    sFilename_flowline_hydroshed_tmp = sFilename_flowline_hydroshed_out.replace('.geojson', '_outlet.geojson')
    export_flowline_to_geojson(aFlowline_hydroshed_outlet_simplified, sFilename_flowline_hydroshed_tmp)

    # Sort the outlet flowlines by drainage area (largest to smallest)
    aFlowline_hydroshed_outlet_simplified = sorted(
        aFlowline_hydroshed_outlet_simplified,
        key=lambda x: x.dDrainage_area,
        reverse=True  # Descending order (largest first)
        )

    aFlowlineID_outlet = np.array(aFlowlineID_outlet)

    #collect the flowlines that flow in to the outlet flowlines
    aFlowline_all = list()
    for pFlowline in aFlowline_hydroshed_outlet_simplified:
        lFlowlineID = pFlowline.lFlowlineID
        aFlowline_all.append(pFlowline)
        dummy_index = np.where(aFlowlineID_outlet == lFlowlineID)
        for i in range(len(dummy_index[0])):
            pFlowline_current = aFlowline_hydroshed_upstream_all[dummy_index[0][i]]
            pFlowline_current.lFlowlineIndex = lFlowlineIndex
            aFlowline_all.append(pFlowline_current)
            pass

    #save the flowlines
    sFilename_flowline_hydroshed_tmp = sFilename_flowline_hydroshed_out.replace('.geojson', '_all.geojson')
    aAttribute_field = ['lineid', 'downstream_id', 'drainage_area']
    aAttribute_dtype= ['int', 'int', 'float']
    aAttribute_data=list()
    aAttribute_data.append( np.array([pFlowline.lFlowlineID for pFlowline in aFlowline_all]) )
    aAttribute_data.append( np.array([pFlowline.lFlowlineID_downstream for pFlowline in aFlowline_all]) )
    aAttribute_data.append( np.array([pFlowline.dDrainage_area for pFlowline in aFlowline_all]) )
    export_flowline_to_geojson(aFlowline_all, sFilename_flowline_hydroshed_tmp,
                                aAttribute_field=aAttribute_field,
    aAttribute_data=aAttribute_data,
    aAttribute_dtype=aAttribute_dtype)
    print('Number of all flowlines in the hydroshed: ', len(aFlowline_all))

    def is_downstream(iStream_segment_a, iStream_segment_b):
        """
        Check if flowline a is downstream of flowline b.
        Args:
            iStream_segment_a: Stream segment ID of flowline a
            iStream_segment_b: Stream segment ID of flowline b
        Returns:
            1 if a is downstream of b, 0 otherwise
        """
        if iStream_segment_a == iStream_segment_b:
            return 1
        else:
            index_current = np.where(aStream_segment == iStream_segment_a)
            lStream_segment_index_next = aFlowline_basin_simplified[index_current[0][0]].lFlowlineIndex_downstream
            while lStream_segment_index_next !=-1 and lStream_segment_index_next is not None:
                lStream_segment_next= aFlowline_basin_simplified[lStream_segment_index_next].iStream_segment
                if lStream_segment_next == iStream_segment_b:
                    return 1
                else:
                    lStream_segment_index_next = aFlowline_basin_simplified[lStream_segment_index_next].lFlowlineIndex_downstream
            return 0

    #now we will use the pyflowline package algorithm to simplify the flowlines
    nFlowline_outlet = len(aFlowline_hydroshed_outlet_simplified)
    def find_index_flowline_list(aFlowline_list, target_stream_segment):
        """
        Find the index of a flowline in a list by its ID

        Args:
            aFlowline_list: List of flowline objects
            target_stream_segment: The flowline ID to search for

        Returns:
            int: Index of the flowline if found, -1 if not found
        """
        for i, pFlowline in enumerate(aFlowline_list):
            if hasattr(pFlowline, 'iStream_segment') and pFlowline.iStream_segment == target_stream_segment:
                return i
        return -1

    def remove_flowline_by_id(aFlowline_rtree, iStream_segment_in):
        """
        Remove a flowline from the list by its ID

        Args:
            aFlowline_rtree: List of flowline objects
            iStream_segment_in: The ID of the flowline to remove

        Returns:
            bool: True if removed successfully, False if not found
        """
        for i, pFlowline in enumerate(aFlowline_rtree):
            if hasattr(pFlowline, 'iStream_segment') and pFlowline.iStream_segment == iStream_segment_in:
                removed_flowline = aFlowline_rtree.pop(i)
                print(f"Removed flowline with ID {iStream_segment_in}")
                return True

        print(f"Flowline with ID {iStream_segment_in} not found in aFlowline_rtree")
        return False

    def remove_from_rtree_by_id_and_bounds(index_reach, item_id, bounding_box):
        """
        Remove an item from R-tree using ID and bounding box

        Args:
            index_reach: The R-tree index object
            item_id: The ID that was used when inserting the item
            bounding_box: The bounding box (minx, miny, maxx, maxy) used during insertion

        Returns:
            bool: True if removed successfully, False otherwise
        """
        try:
            # R-tree delete requires both the ID and the exact bounding box
            index_reach.remove(item_id, bounding_box)
            print(f"Successfully removed item {item_id} from R-tree")
            return True
        except Exception as e:
            print(f"Failed to remove item {item_id} from R-tree: {e}")
            return False

    def tag_upstream(iStream_segment_in, dDrainage_area_threshold):
        lIndex_dummy= np.where(aStream_segment == iStream_segment_in)
        #check whther it is empty here
        if len(lIndex_dummy[0]) == 0:
            print('This flowline id does not exist: ', iStream_segment_in)
            return
        lFlowlineIndex = lIndex_dummy[0][0]
        pFlowline_curent = aFlowline_basin_simplified[lFlowlineIndex]
        aUpstream_segment = pFlowline_curent.aFlowline_upstream
        if aUpstream_segment is None:
            return
        nUpstream = len(aUpstream_segment)
        aUpstream_segment = np.array(aUpstream_segment)
        if nUpstream > 0:
            if nUpstream == 1:
                #check whether it intersects with existing any existing flowlines in the rtree
                iStream_segment_a = aUpstream_segment[0]
                index_current= np.where(aStream_segment == iStream_segment_a)
                pFlowline_a = aFlowline_basin_simplified[index_current[0][0]]
                #get its bound
                pBound = all_bounds_cache[iStream_segment_a]
                aIntersect = list(index_reach.search(pBound))
                nIntersect = len(aIntersect)
                if nIntersect > 0: #there are some flowlines that intersect with the current flowline
                    iFlag_keep = 1
                    for j in range(nIntersect):
                        iStream_segment_b = aIntersect[j]
                        idx = find_index_flowline_list(aFlowline_rtree, iStream_segment_b)
                        pFlowline_b = aFlowline_rtree[idx] #this flowline is already in the rtree
                        #iStream_segment_b = pFlowline_b.iStream_segment
                        if is_downstream(iStream_segment_a, iStream_segment_b) == 1:
                            #this is the flowline that we are looking for
                            continue
                        else:
                            dDistance = pFlowline_a.calculate_distance_to_flowline(pFlowline_b )
                            if dDistance < dDistance_tolerance_in:
                                iFlag_keep = 0
                                print('Flowline ', iStream_segment_a, ' intersects with flowline ', iStream_segment_b, ' with distance: ', dDistance)
                            else:
                                pass
                    #add it into the index tree
                    if iFlag_keep == 1:
                        aFlowline_rtree.append(pFlowline_a)
                        index_reach.insert(iStream_segment_a, pBound)
                        if iStream_segment_a ==624:
                            print('debug')
                        tag_upstream(iStream_segment_a, dDrainage_area_threshold)
                    else:
                        pass
                else: #no intersecting flowlines
                    aFlowline_rtree.append(pFlowline_a)
                    index_reach.insert(iStream_segment_a, pBound)
                    if iStream_segment_a ==624:
                        print('debug')
                    tag_upstream(iStream_segment_a, dDrainage_area_threshold)
                    pass
            else:
                if nUpstream >= 2:
                    #this is a confluence, we need to check the distance
                    #always start with the first upstream that has the largest drainage area
                    aDrainage_area = list()
                    aStream_order = list()
                    aStream_segment_confluence = list()
                    for k in range(nUpstream):
                        iSegment_upstream = aUpstream_segment[k] #segment
                        index_current = np.where(aStream_segment == iSegment_upstream)
                        pFlowline_a = aFlowline_basin_simplified[index_current[0][0]]
                        aDrainage_area.append(pFlowline_a.dDrainage_area)
                        aStream_order.append(pFlowline_a.iStream_order)
                        aStream_segment_confluence.append(pFlowline_a.iStream_segment)

                    aDrainage_area = np.array(aDrainage_area)
                    aStream_order = np.array(aStream_order)
                    aIndex_sorted = np.argsort(aDrainage_area)[::-1]
                    # To get the sorted original indices:
                    sorted_indices = aUpstream_segment[aIndex_sorted]
                    # now we can start the loop with the sorted indices
                    aUpstream_segment = sorted_indices
                    aStream_order = aStream_order[aIndex_sorted]
                    for k in range(nUpstream):
                        #repeat the process for each upstream flowline
                        iSegment_upstream = aUpstream_segment[k]
                        dDrainage_area_upstream = aDrainage_area[k]
                        index_current = np.where(aStream_segment == iSegment_upstream)
                        pFlowline_a = aFlowline_basin_simplified[index_current[0][0]]
                        iStream_order_a = pFlowline_a.iStream_order
                        dDrainage_area_a = pFlowline_a.dDrainage_area
                        pBound_a = all_bounds_cache[iSegment_upstream]
                        if iStream_order_a == 1:
                            aIntersect = list(index_reach.search(pBound_a))
                            nIntersect = len(aIntersect)
                            if nIntersect>0:
                                iFlag_keep = 1
                                for j in range(len(aIntersect)):
                                    iStream_segment_b = aIntersect[j]
                                    idx = find_index_flowline_list(aFlowline_rtree, iStream_segment_b)
                                    pFlowline_b = aFlowline_rtree[idx]
                                    if iStream_segment_b in aStream_segment_confluence or is_downstream(iSegment_upstream, iStream_segment_b) == 1:
                                        #this is the flowline that we are looking for
                                        continue
                                    else:
                                        dDistance = pFlowline_a.calculate_distance_to_flowline( pFlowline_b )
                                        if dDistance < dDistance_tolerance_in:
                                            iFlag_keep = 0
                                            print('Flowline ', iSegment_upstream, ' intersects with flowline ', iStream_segment_b, ' with distance: ', dDistance)
                                        else:
                                            pass
                                #add it into the index tree
                                if iFlag_keep == 1:
                                    aFlowline_rtree.append(pFlowline_a)
                                    index_reach.insert(iSegment_upstream, pBound_a)
                                    if iSegment_upstream ==624:
                                        print('debug')
                                    tag_upstream(iSegment_upstream, dDrainage_area_threshold)
                                else:
                                    #aFlowline_all[index_current].iFlag_keep = 0
                                    print('Flowline ', iSegment_upstream, ' is not kept due to intersection with other flowlines.')
                                    pass

                            else: #no intersecting flowlines
                                aFlowline_rtree.append(pFlowline_a)
                                index_reach.insert(iSegment_upstream, pBound_a)
                                if iSegment_upstream ==624:
                                    print('debug')
                                tag_upstream(iSegment_upstream, dDrainage_area_threshold)
                                pass
                            pass
                        else:
                            aIntersect = list(index_reach.search(pBound_a))
                            nIntersect = len(aIntersect)
                            if nIntersect>0:
                                iFlag_keep = 1
                                for j in range(len(aIntersect)):
                                    iStream_segment_b = aIntersect[j]
                                    idx = find_index_flowline_list(aFlowline_rtree, iStream_segment_b)
                                    pFlowline_b = aFlowline_rtree[idx]
                                    iStream_order_b = pFlowline_b.iStream_order
                                    dDrainage_area_b = pFlowline_b.dDrainage_area
                                    if iStream_segment_b in aStream_segment_confluence or is_downstream(iSegment_upstream, iStream_segment_b) == 1:
                                        #this is the flowline that we are looking for
                                        continue
                                    else:
                                        if dDrainage_area_a >=dDrainage_area_threshold:
                                            #we can keep this flowline
                                            iFlag_keep = 1
                                            #how about the other flowline? we dont need to check the distance
                                            #and we dont need to remove others for now
                                        else:
                                            dDistance = pFlowline_a.calculate_distance_to_flowline( pFlowline_b )
                                            if dDistance < dDistance_tolerance_in:
                                                if dDrainage_area_b < dDrainage_area_a: #the other flowline is a smaller one, we can remove it?
                                                    success = remove_flowline_by_id(aFlowline_rtree, iStream_segment_b)
                                                    pBound_b = all_bounds_cache[iStream_segment_b]
                                                    success2 = remove_from_rtree_by_id_and_bounds(index_reach, iStream_segment_b, pBound_b)
                                                    print(f"Removed flowline {iStream_segment_b} from aFlowline_rtree and index_reach: {success} {success2}")
                                                    pass
                                                else:
                                                    iFlag_keep = 0
                                                    print('Flowline ', iSegment_upstream, ' intersects with flowline ', iStream_segment_b, ' with distance: ', dDistance)

                                #add it into the index tree
                                if iFlag_keep == 1:
                                    aFlowline_rtree.append(pFlowline_a)
                                    if iSegment_upstream ==624:
                                        print('debug')
                                    index_reach.insert(iSegment_upstream, pBound_a)
                                    tag_upstream(iSegment_upstream, dDrainage_area_threshold)
                                else:
                                    print('Flowline ', iSegment_upstream, ' is not kept due to intersection with other flowlines.')
                                    pass
                            else:
                                aFlowline_rtree.append(pFlowline_a)
                                if iSegment_upstream ==624:
                                    print('debug')
                                index_reach.insert(iSegment_upstream, pBound_a)
                                tag_upstream(iSegment_upstream, dDrainage_area_threshold)
                    pass
                else:
                    # more than 2?
                    print('This is a confluence with more than 2 upstream flowlines, which is not supported yet.')
                    pass
        else:
            #if a flowline has no upstream, then it is a headwater
            pass

    aFlowline_upstream_simplified = list()
    nBasin = nFlowline_outlet
    aOulet_coordate= np.full( (nBasin, 2), -9999, dtype=float)

    #create a configuration file
    #this configuration will be used for the pyflowline standalone simulation,
    #only the largest basin will be used for the simulation
    sWorkspace_output = os.path.dirname(sFilename_flowline_hydroshed_out)
    if iFlag_pyflowline_configuration_in ==1:
        sFilename_configuration_json = os.path.join(sWorkspace_output, 'pyflowline_configuration.json')
        create_template_configuration_file(sFilename_configuration_json,
            sWorkspace_output = sWorkspace_output,
            iFlag_standalone_in=1,
            nOutlet = nOutlet_largest,
            sMesh_type_in='mpas',
            sModel_in='pyflowline')
        sFilename_configuration_basin_json = os.path.join(sWorkspace_output, 'pyflowline_configuration_basins.json')

    aFlowline_rtree_all = list()
    aFlowline_rtree = list()
    for i in range(0, nFlowline_outlet, 1):
        sFilename_flowline_hydroshed_outlet = sFilename_flowline_hydroshed_out.replace('.geojson', f'_{i+1:04d}_outlet_simplified.geojson')
        aFlowline_rtree.clear()

        pFlowline_current = aFlowline_hydroshed_outlet_simplified[i]
        dLongitude_outlet = pFlowline_current.pVertex_end.dLongitude_degree
        dLatitude_outlet = pFlowline_current.pVertex_end.dLatitude_degree
        dDrainage_area_threshold = pFlowline_current.dDrainage_area * dDrainage_area_threshold_ratio #this is a key threshold to keep flowline with large drainage area even they are close
        aOulet_coordate[i, 0] = dLongitude_outlet
        aOulet_coordate[i, 1] = dLatitude_outlet
        pVertex_outlet = pFlowline_current.pVertex_end
        aVertex = find_flowline_vertex(aFlowline_all)
        #aFlowline_basin_simplified = split_flowline(aFlowline_all, aVertex)
        aFlowline_basin_simplified = update_head_water_stream_order(aFlowline_all)
        aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity, pVertex_outlet = find_flowline_confluence(aFlowline_basin_simplified, pVertex_outlet)
        aFlowline_basin_simplified = merge_flowline( aFlowline_basin_simplified, aVertex, pVertex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence )
        aFlowline_basin_simplified = update_head_water_stream_order(aFlowline_basin_simplified )
        aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity, pVertex_outlet = find_flowline_confluence(aFlowline_basin_simplified, pVertex_outlet)
        aFlowline_basin_simplified, aStream_segment = define_stream_segment_index(aFlowline_basin_simplified)
        if len(aFlowline_basin_simplified) == 1:
            aFlowline_rtree.append(aFlowline_basin_simplified[0])
            pass
        else:
            aVertex = np.array(aVertex)
            aIndex_confluence = np.array(aIndex_confluence)
            aVertex_confluence = aVertex[aIndex_confluence]
            aConfluence_basin_simplified = basin_build_confluence(aFlowline_basin_simplified, aVertex_confluence)
            aFlowline_basin_simplified = define_stream_topology(aFlowline_basin_simplified, aConfluence_basin_simplified)
            aFlowline_basin_simplified, aStream_order = define_stream_order(aFlowline_basin_simplified, aConfluence_basin_simplified, iFlag_so_method_in=1)

            #Add index to the filename for multiple basin support
            aStream_segment = np.array(aStream_segment)
            aStream_order = np.array(aStream_order)
            index_reach = RTree(max_cap=5, min_cap=2)
            all_bounds_cache = precompute_flowline_geometries_by_segment(aFlowline_basin_simplified, dDistance_tolerance_in)

            pFlowline_outlet = aFlowline_basin_simplified[0]
            pBound = all_bounds_cache[pFlowline_outlet.iStream_segment]
            index_reach.insert(pFlowline_outlet.iStream_segment, pBound)
            aFlowline_rtree.append(pFlowline_outlet)
            tag_upstream(pFlowline_outlet.iStream_segment, dDrainage_area_threshold)

        #now save the flowlines

        if i < nOutlet_largest:
            #produce a basin configuration file
            #update the configuration file with the basin information
            aStream_segment=list()
            aStream_order = list()
            for pFlowline in aFlowline_rtree:
                aStream_segment.append(pFlowline.iStream_segment)
                aStream_order.append(pFlowline.iStream_order)

            aStream_segment = np.array(aStream_segment)
            aStream_order = np.array(aStream_order)
            export_flowline_to_geojson(aFlowline_rtree,
                                       sFilename_flowline_hydroshed_outlet,
                    aAttribute_data=[aStream_segment, aStream_order],
                    aAttribute_field=['stream_segment','stream_order'],
                    aAttribute_dtype=['int','int'])
            if pFlowline_current.iFlag_endorheic != 1:
                if iFlag_pyflowline_configuration_in ==1:
                    change_json_key_value(sFilename_configuration_basin_json, 'dAccumulation_threshold', dDrainage_area_threshold, iFlag_basin_in=1, iBasin_index_in=i)
                    change_json_key_value(sFilename_configuration_basin_json, 'dLatitude_outlet_degree', dLatitude_outlet, iFlag_basin_in=1, iBasin_index_in=i)
                    change_json_key_value(sFilename_configuration_basin_json, 'dLongitude_outlet_degree', dLongitude_outlet, iFlag_basin_in=1, iBasin_index_in=i)
                    change_json_key_value(sFilename_configuration_basin_json, 'sFilename_flowline_filter', sFilename_flowline_hydroshed_outlet, iFlag_basin_in=1, iBasin_index_in=i)
            else:
                print('This is an endorheic basin, we do not need to save the basin configuration file for it.')
            pass


        for pFlowline in aFlowline_rtree:
            aFlowline_rtree_all.append(pFlowline)


    #save the flowlines
    #aFlowline_after_distance_operation = aFlowline_hydroshed_outlet_simplified + aFlowline_upstream_simplified
    export_flowline_to_geojson(aFlowline_rtree_all, sFilename_flowline_hydroshed_out)
    #close the file
    pDataset_in = pLayer_shapefile = pFeature_shapefile = None
    print('Number of flowlines in the hydroshed: ', len(aFlowline_upstream_simplified))
    return

