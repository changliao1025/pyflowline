
import os
import importlib.util
import numpy as np
from osgeo import ogr, osr
from rtree.index import Index as RTreeindex
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.geometry.convert_idl_polygon_to_valid_polygon import convert_idl_polygon_to_valid_polygon

from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_cell
from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_flowline

iFlag_cython = importlib.util.find_spec("cython")
if iFlag_cython is not None:
    from pyflowline.algorithms.cython.kernel import find_vertex_in_list
else:
    from pyflowline.algorithms.auxiliary.find_vertex_in_list import find_vertex_in_list


def intersect_flowline_with_mesh(iMesh_type_in, sFilename_mesh_in, sFilename_flowline_in, sFilename_output_in):
    if  os.path.exists(sFilename_mesh_in) and  os.path.exists(sFilename_flowline_in) :
        pass
    else:
        print('The input file does not exist')
        return

    if os.path.exists(sFilename_output_in):
        os.remove(sFilename_output_in)


    pDriver_geojson = ogr.GetDriverByName( "GeoJSON")
    #aCell=list()
    aCell_intersect=list()

    pDataset_mesh = pDriver_geojson.Open(sFilename_mesh_in, 0)
    pLayer_mesh = pDataset_mesh.GetLayer(0)
    pSpatial_reference_mesh = pLayer_mesh.GetSpatialRef()
    nfeature_mesh = pLayer_mesh.GetFeatureCount()

    pDataset_flowline = pDriver_geojson.Open(sFilename_flowline_in, 0)
    pLayer_flowline = pDataset_flowline.GetLayer(0)
    pSpatial_reference_flowline = pLayer_flowline.GetSpatialRef()
    nfeature_flowline = pLayer_flowline.GetFeatureCount()
    pLayerDefinition = pLayer_flowline.GetLayerDefn()

    comparison = pSpatial_reference_mesh.IsSame(pSpatial_reference_flowline)
    if(comparison != 1):
        iFlag_transform = 1
        transform = osr.CoordinateTransformation(pSpatial_reference_mesh, pSpatial_reference_flowline)
    else:
        iFlag_transform = 0

    pDataset_out = pDriver_geojson.CreateDataSource(sFilename_output_in)

    pLayerOut = pDataset_out.CreateLayer('flowline', pSpatial_reference_flowline, ogr.wkbMultiLineString)
    # Add one attribute
    pLayerOut.CreateField(ogr.FieldDefn('lineid', ogr.OFTInteger64)) #long type for high resolution
    pLayerOut.CreateField(ogr.FieldDefn('stream_segment', ogr.OFTInteger)) #long type for high resolution
    pLayerOut.CreateField(ogr.FieldDefn('stream_order', ogr.OFTInteger)) #long type for high resolution
    pLayerDefn = pLayerOut.GetLayerDefn()
    pFeatureOut = ogr.Feature(pLayerDefn)

    lFlowlineID = 1
    aFlowline_intersect_all=list()

    #interleaved = True
    #index_flowline = RTree(interleaved=interleaved, max_cap=5, min_cap=2)
    index_flowline = RTreeindex()
    for i in range(nfeature_flowline):
        lID = i
        pFeature_flowline = pLayer_flowline.GetFeature(i)
        pGeometry_flowline = pFeature_flowline.GetGeometryRef()
        left, right, bottom, top= pGeometry_flowline.GetEnvelope()
        pBound= (left, bottom, right, top)
        index_flowline.insert(lID, pBound)  #
    #now intersect using rtree
    for j in range (nfeature_mesh):
        pFeature_mesh = pLayer_mesh.GetFeature(j)
        pGeometry_mesh = pFeature_mesh.GetGeometryRef()
        aCoords_gcs = get_geometry_coordinates(pGeometry_mesh)
        lCellID = pFeature_mesh.GetField("cellid")
        dLongitude_center = pFeature_mesh.GetField("longitude")
        dLatitude_center = pFeature_mesh.GetField("latitude")
        dArea = pFeature_mesh.GetField("area")
        dLongitude_min= np.min(aCoords_gcs[:,0])
        dLongitude_max = np.max(aCoords_gcs[:,0])
        dArea = pFeature_mesh.GetField("area")
        left_orig, right_orig, bottom, top = pGeometry_mesh.GetEnvelope()
        pGeometrytype_mesh = pGeometry_mesh.GetGeometryName()
        if(pGeometrytype_mesh == 'POLYGON'):
            pCell = convert_gcs_coordinates_to_cell(iMesh_type_in, dLongitude_center, dLatitude_center, aCoords_gcs)
            pCell.lCellID = lCellID
            pCell.dArea = dArea
            pCell.dLength = pCell.calculate_edge_length()
            pCell.dLength_flowline = pCell.dLength
            aFlowline_intersect = list()
            iFlag_intersected = 0
            #left, right, bottom, top= pGeometry_mesh.GetEnvelope()
            pBound= (left_orig, bottom, right_orig, top)
            aIntersect = list(index_flowline.intersection(pBound))
            for k in aIntersect:
                pFeature_flowline = pLayer_flowline.GetFeature(k)
                pGeometry_flowline = pFeature_flowline.GetGeometryRef()
                iFlag_intersect = pGeometry_flowline.Intersects( pGeometry_mesh )
                if( iFlag_intersect == True):
                    iFlag_intersected = 1
                    pGeometry_intersect = pGeometry_flowline.Intersection(pGeometry_mesh)
                    pGeometrytype_intersect = pGeometry_intersect.GetGeometryName()
                    iStream_segment = pFeature_flowline.GetField("stream_segment")
                    iStream_order = pFeature_flowline.GetField("stream_order")
                    if pGeometrytype_intersect == 'LINESTRING':
                        pFeatureOut.SetGeometry(pGeometry_intersect)
                        pFeatureOut.SetField("lineid", lFlowlineID)
                        pFeatureOut.SetField("stream_segment", iStream_segment)
                        pFeatureOut.SetField("stream_order", iStream_order)
                        pLayerOut.CreateFeature(pFeatureOut)
                        aCoords = list()
                        for i in range(0, pGeometry_intersect.GetPointCount()):
                            pt = pGeometry_intersect.GetPoint(i)
                            aCoords.append( [ pt[0], pt[1]])
                        dummy1= np.array(aCoords)
                        pFlowline = convert_gcs_coordinates_to_flowline(dummy1)
                        pFlowline.calculate_length()
                        pFlowline.lFlowlineIndex = lFlowlineID
                        pFlowline.iStream_segment = iStream_segment
                        pFlowline.iStream_order = iStream_order
                        aFlowline_intersect.append(pFlowline)
                        aFlowline_intersect_all.append(pFlowline)
                        lFlowlineID = lFlowlineID + 1
                    else:
                        if(pGeometrytype_intersect == 'MULTILINESTRING'):
                            nLine = pGeometry_intersect.GetGeometryCount()
                            for i in range(nLine):
                                Line = pGeometry_intersect.GetGeometryRef(i)
                                pFeatureOut.SetGeometry(Line)
                                pFeatureOut.SetField("lineid", lFlowlineID)
                                pFeatureOut.SetField("stream_segment", iStream_segment)
                                pFeatureOut.SetField("stream_order", iStream_order)
                                pLayerOut.CreateFeature(pFeatureOut)
                                aCoords = list()
                                for i in range(0, Line.GetPointCount()):
                                    pt = Line.GetPoint(i)
                                    aCoords.append( [ pt[0], pt[1]])
                                dummy1= np.array(aCoords)
                                pFlowline = convert_gcs_coordinates_to_flowline(dummy1)
                                pFlowline.calculate_length()
                                pFlowline.lFlowlineIndex = lFlowlineID
                                pFlowline.iStream_segment = iStream_segment
                                pFlowline.iStream_order = iStream_order
                                aFlowline_intersect.append(pFlowline)
                                aFlowline_intersect_all.append(pFlowline)
                                lFlowlineID = lFlowlineID + 1
                            pass
                        else:
                            pass
                    pass
            #now add back to the cell object
            pCell.aFlowline = aFlowline_intersect
            pCell.nFlowline = len(aFlowline_intersect)
            if iFlag_intersected ==1:
                pCell.iFlag_intersected = 1
                pCell.dLength_flowline = 0.0 #reset the flowline length
                for i in range (pCell.nFlowline):
                    pFlowline = pCell.aFlowline[i]
                    dLength_flowline = pFlowline.dLength
                    if ( dLength_flowline > pCell.dLength_flowline ):
                        pCell.dLength_flowline = dLength_flowline
                #replace flowline length if there is an actual flowline
                aCell_intersect.append(pCell)
            else:
                pCell.iFlag_intersected = 0
                pass

    #quality control
    #find the flowline that has the largest stream segment index
    pVertex_outlet = None
    iSegment_max = 0
    iIndex_flowline = -1
    for i in range(nfeature_flowline):
        pFeature_flowline = pLayer_flowline.GetFeature(i)
        iStream_segment = pFeature_flowline.GetField("stream_segment")
        if (iStream_segment > iSegment_max):
            iSegment_max = iStream_segment
            iIndex_flowline = i

    if (iIndex_flowline != -1):
        #find all the flowlines and cells that have the same stream segment

        #get the first vertex of the flowline
        pFeature_flowline = pLayer_flowline.GetFeature(iIndex_flowline)
        pGeometry_flowline = pFeature_flowline.GetGeometryRef()
        aCoords = list()
        for i in range(0, pGeometry_flowline.GetPointCount()):
            pt = pGeometry_flowline.GetPoint(i)
            aCoords.append( [ pt[0], pt[1]])
        dummy1= np.array(aCoords)
        pFlowline = convert_gcs_coordinates_to_flowline(dummy1)
        pVertex_start = pFlowline.pVertex_start
        pVertex_outlet = pFlowline.pVertex_end
        aFlowline_last = list()
        nFlowline_intersect = len(aFlowline_intersect_all)
        for i in range(nFlowline_intersect):
            pFlowline = aFlowline_intersect_all[i]
            if ( pFlowline.iStream_segment == iSegment_max):
                aFlowline_last.append(pFlowline)


        aFlowline_stay = list()
        #now check the connectivity of the flowlines because the flowlines are not necessarily connected
        iFlag_done = 0
        while (iFlag_done ==0):
            iFlag_found = 0
            nFlowline_last = len(aFlowline_last)
            for i in range(nFlowline_last):
                pFlowline = aFlowline_last[i]
                pVertex_start_dummy = pFlowline.pVertex_start
                pVertex_end_dummy = pFlowline.pVertex_end
                if pVertex_start_dummy == pVertex_start:
                    #found the head
                    iFlag_found = 1
                    pVertex_start = pVertex_end_dummy
                    aFlowline_stay.append(i)
                    pVertex_outlet = pVertex_end_dummy
                    break

            if iFlag_found == 1:
                iFlag_done = 0
            else:
                iFlag_done = 1

    return   aCell_intersect, aFlowline_intersect_all, pVertex_outlet


def intersect_flowline_with_mesh_optimized(iMesh_type_in, sFilename_mesh_in, sFilename_flowline_in, sFilename_output_in, optimization_method='auto'):
    """
    Optimized version of intersect_flowline_with_mesh with automatic method selection

    When mesh cells >> flowlines, builds R-tree for mesh and loops through flowlines (much faster)
    When flowlines >= mesh/2, uses original flowline R-tree approach

    Args:
        iMesh_type_in: Mesh type identifier
        sFilename_mesh_in: Path to mesh GeoJSON file
        sFilename_flowline_in: Path to flowline GeoJSON file
        sFilename_output_in: Path to output GeoJSON file
        optimization_method: 'auto', 'mesh_rtree', or 'flowline_rtree'

    Returns:
        tuple: (aCell_intersect, aFlowline_intersect_all, pVertex_outlet)
    """
    import time
    start_time = time.time()

    if not (os.path.exists(sFilename_mesh_in) and os.path.exists(sFilename_flowline_in)):
        print('The input file does not exist')
        return

    if os.path.exists(sFilename_output_in):
        os.remove(sFilename_output_in)

    pDriver_geojson = ogr.GetDriverByName("GeoJSON")
    aCell_intersect = list()

    # Open datasets
    pDataset_mesh = pDriver_geojson.Open(sFilename_mesh_in, 0)
    pLayer_mesh = pDataset_mesh.GetLayer(0)
    pSpatial_reference_mesh = pLayer_mesh.GetSpatialRef()
    nfeature_mesh = pLayer_mesh.GetFeatureCount()

    pDataset_flowline = pDriver_geojson.Open(sFilename_flowline_in, 0)
    pLayer_flowline = pDataset_flowline.GetLayer(0)
    pSpatial_reference_flowline = pLayer_flowline.GetSpatialRef()
    nfeature_flowline = pLayer_flowline.GetFeatureCount()

    print(f"Optimization analysis: {nfeature_mesh} mesh cells, {nfeature_flowline} flowlines")

    # Automatic method selection
    if optimization_method == 'auto':
        if nfeature_mesh > nfeature_flowline * 2:
            optimization_method = 'mesh_rtree'
            print("Using MESH R-tree optimization (mesh >> flowlines)")
        else:
            optimization_method = 'flowline_rtree'
            print("Using FLOWLINE R-tree optimization (original method)")
    else:
        print(f"Using specified method: {optimization_method}")

    # Handle coordinate transformations
    comparison = pSpatial_reference_mesh.IsSame(pSpatial_reference_flowline)
    if comparison != 1:
        iFlag_transform = 1
        transform = osr.CoordinateTransformation(pSpatial_reference_mesh, pSpatial_reference_flowline)
    else:
        iFlag_transform = 0

    # Setup output
    pDataset_out = pDriver_geojson.CreateDataSource(sFilename_output_in)
    pLayerOut = pDataset_out.CreateLayer('flowline', pSpatial_reference_flowline, ogr.wkbMultiLineString)
    pLayerOut.CreateField(ogr.FieldDefn('lineid', ogr.OFTInteger64))
    pLayerOut.CreateField(ogr.FieldDefn('stream_segment', ogr.OFTInteger))
    pLayerOut.CreateField(ogr.FieldDefn('stream_order', ogr.OFTInteger))
    pLayerDefn = pLayerOut.GetLayerDefn()
    pFeatureOut = ogr.Feature(pLayerDefn)

    lFlowlineID = 1
    aFlowline_intersect_all = list()

    if optimization_method == 'mesh_rtree':
        # OPTIMIZED: Build R-tree for mesh, loop through flowlines
        aCell_intersect, aFlowline_intersect_all, lFlowlineID = _intersect_mesh_rtree_method(
            iMesh_type_in, pLayer_mesh, pLayer_flowline, pLayerOut, pFeatureOut,
            nfeature_mesh, nfeature_flowline, lFlowlineID
        )
    else:
        # ORIGINAL: Build R-tree for flowlines, loop through mesh
        aCell_intersect, aFlowline_intersect_all, lFlowlineID = _intersect_flowline_rtree_method(
            iMesh_type_in, pLayer_mesh, pLayer_flowline, pLayerOut, pFeatureOut,
            nfeature_mesh, nfeature_flowline, lFlowlineID
        )

    # Quality control (same as original)
    pVertex_outlet = _perform_quality_control(pLayer_flowline, aFlowline_intersect_all, nfeature_flowline)

    end_time = time.time()
    print(f"Intersection completed in {end_time - start_time:.2f} seconds")
    print(f"Found {len(aFlowline_intersect_all)} total intersections")

    return aCell_intersect, aFlowline_intersect_all, pVertex_outlet


def _intersect_mesh_rtree_method(iMesh_type_in, pLayer_mesh, pLayer_flowline, pLayerOut, pFeatureOut,
                                nfeature_mesh, nfeature_flowline, lFlowlineID):
    """
    OPTIMIZED METHOD: Build R-tree for mesh cells, loop through flowlines
    Much faster when mesh >> flowlines
    """
    aCell_intersect = list()
    aFlowline_intersect_all = list()

    print("Building mesh R-tree index...")

    # Build R-tree for mesh cells
    index_mesh = RTreeindex()
    mesh_data = {}  # Cache mesh data

    for j in range(nfeature_mesh):
        pFeature_mesh = pLayer_mesh.GetFeature(j)
        pGeometry_mesh = pFeature_mesh.GetGeometryRef()
        left_orig, right_orig, bottom, top = pGeometry_mesh.GetEnvelope()
        pBound = (left_orig, bottom, right_orig, top)
        index_mesh.insert(j, pBound)

        # Cache mesh data to avoid repeated feature access
        mesh_data[j] = {
            'feature': pFeature_mesh,
            'geometry': pGeometry_mesh,
            'coords_gcs': get_geometry_coordinates(pGeometry_mesh),
            'cell_id': pFeature_mesh.GetField("cellid"),
            'longitude': pFeature_mesh.GetField("longitude"),
            'latitude': pFeature_mesh.GetField("latitude"),
            'area': pFeature_mesh.GetField("area")
        }

    print("Processing flowline intersections...")
    processed_cells = {}  # Track processed cells to avoid duplicates

    # Loop through flowlines (fewer iterations!)
    for i in range(nfeature_flowline):
        pFeature_flowline = pLayer_flowline.GetFeature(i)
        pGeometry_flowline = pFeature_flowline.GetGeometryRef()
        left, right, bottom, top = pGeometry_flowline.GetEnvelope()
        pBound_flowline = (left, bottom, right, top)

        iStream_segment = pFeature_flowline.GetField("stream_segment")
        iStream_order = pFeature_flowline.GetField("stream_order")

        # Find potentially intersecting mesh cells
        aIntersect_mesh = list(index_mesh.intersection(pBound_flowline))

        for mesh_idx in aIntersect_mesh:
            data = mesh_data[mesh_idx]
            pGeometry_mesh = data['geometry']

            # Check actual intersection
            iFlag_intersect = pGeometry_flowline.Intersects(pGeometry_mesh)
            if iFlag_intersect:
                # Get or create cell object
                lCellID = data['cell_id']
                if lCellID not in processed_cells:
                    pGeometrytype_mesh = pGeometry_mesh.GetGeometryName()
                    if pGeometrytype_mesh == 'POLYGON':
                        pCell = convert_gcs_coordinates_to_cell(iMesh_type_in, data['longitude'],
                                                              data['latitude'], data['coords_gcs'])
                        pCell.lCellID = lCellID
                        pCell.dArea = data['area']
                        pCell.dLength = pCell.calculate_edge_length()
                        pCell.dLength_flowline = pCell.dLength
                        pCell.aFlowline = list()
                        pCell.nFlowline = 0
                        pCell.iFlag_intersected = 0
                        processed_cells[lCellID] = pCell
                        aCell_intersect.append(pCell)
                else:
                    pCell = processed_cells[lCellID]

                # Process intersection geometry
                pGeometry_intersect = pGeometry_flowline.Intersection(pGeometry_mesh)
                pGeometrytype_intersect = pGeometry_intersect.GetGeometryName()

                if pGeometrytype_intersect == 'LINESTRING':
                    lFlowlineID = _process_linestring_intersection(
                        pGeometry_intersect, pFeatureOut, pLayerOut, lFlowlineID,
                        iStream_segment, iStream_order, pCell, aFlowline_intersect_all
                    )
                elif pGeometrytype_intersect == 'MULTILINESTRING':
                    lFlowlineID = _process_multilinestring_intersection(
                        pGeometry_intersect, pFeatureOut, pLayerOut, lFlowlineID,
                        iStream_segment, iStream_order, pCell, aFlowline_intersect_all
                    )

    # Update cell properties
    for pCell in aCell_intersect:
        pCell.nFlowline = len(pCell.aFlowline)
        if pCell.nFlowline > 0:
            pCell.iFlag_intersected = 1
            pCell.dLength_flowline = 0.0
            for pFlowline in pCell.aFlowline:
                if pFlowline.dLength > pCell.dLength_flowline:
                    pCell.dLength_flowline = pFlowline.dLength

    return aCell_intersect, aFlowline_intersect_all, lFlowlineID


def _intersect_flowline_rtree_method(iMesh_type_in, pLayer_mesh, pLayer_flowline, pLayerOut, pFeatureOut,
                                    nfeature_mesh, nfeature_flowline, lFlowlineID):
    """
    ORIGINAL METHOD: Build R-tree for flowlines, loop through mesh cells
    """
    aCell_intersect = list()
    aFlowline_intersect_all = list()

    print("Building flowline R-tree index...")

    # Build R-tree for flowlines (original approach)
    index_flowline = RTreeindex()
    for i in range(nfeature_flowline):
        pFeature_flowline = pLayer_flowline.GetFeature(i)
        pGeometry_flowline = pFeature_flowline.GetGeometryRef()
        left, right, bottom, top = pGeometry_flowline.GetEnvelope()
        pBound = (left, bottom, right, top)
        index_flowline.insert(i, pBound)

    print("Processing mesh intersections...")

    # Loop through mesh cells (original approach)
    for j in range(nfeature_mesh):
        pFeature_mesh = pLayer_mesh.GetFeature(j)
        pGeometry_mesh = pFeature_mesh.GetGeometryRef()
        aCoords_gcs = get_geometry_coordinates(pGeometry_mesh)
        lCellID = pFeature_mesh.GetField("cellid")
        dLongitude_center = pFeature_mesh.GetField("longitude")
        dLatitude_center = pFeature_mesh.GetField("latitude")
        dArea = pFeature_mesh.GetField("area")
        left_orig, right_orig, bottom, top = pGeometry_mesh.GetEnvelope()
        pGeometrytype_mesh = pGeometry_mesh.GetGeometryName()

        if pGeometrytype_mesh == 'POLYGON':
            pCell = convert_gcs_coordinates_to_cell(iMesh_type_in, dLongitude_center, dLatitude_center, aCoords_gcs)
            pCell.lCellID = lCellID
            pCell.dArea = dArea
            pCell.dLength = pCell.calculate_edge_length()
            pCell.dLength_flowline = pCell.dLength
            aFlowline_intersect = list()
            iFlag_intersected = 0

            pBound = (left_orig, bottom, right_orig, top)
            aIntersect = list(index_flowline.intersection(pBound))

            for k in aIntersect:
                pFeature_flowline = pLayer_flowline.GetFeature(k)
                pGeometry_flowline = pFeature_flowline.GetGeometryRef()
                iFlag_intersect = pGeometry_flowline.Intersects(pGeometry_mesh)

                if iFlag_intersect:
                    iFlag_intersected = 1
                    pGeometry_intersect = pGeometry_flowline.Intersection(pGeometry_mesh)
                    pGeometrytype_intersect = pGeometry_intersect.GetGeometryName()
                    iStream_segment = pFeature_flowline.GetField("stream_segment")
                    iStream_order = pFeature_flowline.GetField("stream_order")

                    if pGeometrytype_intersect == 'LINESTRING':
                        lFlowlineID = _process_linestring_intersection(
                            pGeometry_intersect, pFeatureOut, pLayerOut, lFlowlineID,
                            iStream_segment, iStream_order, pCell, aFlowline_intersect_all
                        )
                        aFlowline_intersect.append(pCell.aFlowline[-1])  # Add to local list
                    elif pGeometrytype_intersect == 'MULTILINESTRING':
                        start_count = len(pCell.aFlowline)
                        lFlowlineID = _process_multilinestring_intersection(
                            pGeometry_intersect, pFeatureOut, pLayerOut, lFlowlineID,
                            iStream_segment, iStream_order, pCell, aFlowline_intersect_all
                        )
                        # Add new flowlines to local list
                        aFlowline_intersect.extend(pCell.aFlowline[start_count:])

            # Update cell properties
            pCell.aFlowline = aFlowline_intersect
            pCell.nFlowline = len(aFlowline_intersect)
            if iFlag_intersected == 1:
                pCell.iFlag_intersected = 1
                pCell.dLength_flowline = 0.0
                for pFlowline in pCell.aFlowline:
                    if pFlowline.dLength > pCell.dLength_flowline:
                        pCell.dLength_flowline = pFlowline.dLength
                aCell_intersect.append(pCell)
            else:
                pCell.iFlag_intersected = 0

    return aCell_intersect, aFlowline_intersect_all, lFlowlineID


def _process_linestring_intersection(pGeometry_intersect, pFeatureOut, pLayerOut, lFlowlineID,
                                   iStream_segment, iStream_order, pCell, aFlowline_intersect_all):
    """Helper function to process LINESTRING intersections"""
    pFeatureOut.SetGeometry(pGeometry_intersect)
    pFeatureOut.SetField("lineid", lFlowlineID)
    pFeatureOut.SetField("stream_segment", iStream_segment)
    pFeatureOut.SetField("stream_order", iStream_order)
    pLayerOut.CreateFeature(pFeatureOut)

    aCoords = []
    for i in range(0, pGeometry_intersect.GetPointCount()):
        pt = pGeometry_intersect.GetPoint(i)
        aCoords.append([pt[0], pt[1]])

    dummy1 = np.array(aCoords)
    pFlowline = convert_gcs_coordinates_to_flowline(dummy1)
    pFlowline.calculate_length()
    pFlowline.lFlowlineIndex = lFlowlineID
    pFlowline.iStream_segment = iStream_segment
    pFlowline.iStream_order = iStream_order

    pCell.aFlowline.append(pFlowline)
    aFlowline_intersect_all.append(pFlowline)

    return lFlowlineID + 1


def _process_multilinestring_intersection(pGeometry_intersect, pFeatureOut, pLayerOut, lFlowlineID,
                                        iStream_segment, iStream_order, pCell, aFlowline_intersect_all):
    """Helper function to process MULTILINESTRING intersections"""
    nLine = pGeometry_intersect.GetGeometryCount()
    for i in range(nLine):
        Line = pGeometry_intersect.GetGeometryRef(i)
        pFeatureOut.SetGeometry(Line)
        pFeatureOut.SetField("lineid", lFlowlineID)
        pFeatureOut.SetField("stream_segment", iStream_segment)
        pFeatureOut.SetField("stream_order", iStream_order)
        pLayerOut.CreateFeature(pFeatureOut)

        aCoords = []
        for j in range(0, Line.GetPointCount()):
            pt = Line.GetPoint(j)
            aCoords.append([pt[0], pt[1]])

        dummy1 = np.array(aCoords)
        pFlowline = convert_gcs_coordinates_to_flowline(dummy1)
        pFlowline.calculate_length()
        pFlowline.lFlowlineIndex = lFlowlineID
        pFlowline.iStream_segment = iStream_segment
        pFlowline.iStream_order = iStream_order

        pCell.aFlowline.append(pFlowline)
        aFlowline_intersect_all.append(pFlowline)
        lFlowlineID += 1

    return lFlowlineID


def _perform_quality_control(pLayer_flowline, aFlowline_intersect_all, nfeature_flowline):
    """Perform quality control to find outlet vertex (same as original)"""
    pVertex_outlet = None
    iSegment_max = 0
    iIndex_flowline = -1

    for i in range(nfeature_flowline):
        pFeature_flowline = pLayer_flowline.GetFeature(i)
        iStream_segment = pFeature_flowline.GetField("stream_segment")
        if iStream_segment > iSegment_max:
            iSegment_max = iStream_segment
            iIndex_flowline = i

    if iIndex_flowline != -1:
        pFeature_flowline = pLayer_flowline.GetFeature(iIndex_flowline)
        pGeometry_flowline = pFeature_flowline.GetGeometryRef()
        aCoords = []
        for i in range(0, pGeometry_flowline.GetPointCount()):
            pt = pGeometry_flowline.GetPoint(i)
            aCoords.append([pt[0], pt[1]])
        dummy1 = np.array(aCoords)
        pFlowline = convert_gcs_coordinates_to_flowline(dummy1)
        pVertex_start = pFlowline.pVertex_start
        pVertex_outlet = pFlowline.pVertex_end

        aFlowline_last = []
        nFlowline_intersect = len(aFlowline_intersect_all)
        for i in range(nFlowline_intersect):
            pFlowline = aFlowline_intersect_all[i]
            if pFlowline.iStream_segment == iSegment_max:
                aFlowline_last.append(pFlowline)

        aFlowline_stay = []
        iFlag_done = 0
        while iFlag_done == 0:
            iFlag_found = 0
            nFlowline_last = len(aFlowline_last)
            for i in range(nFlowline_last):
                pFlowline = aFlowline_last[i]
                pVertex_start_dummy = pFlowline.pVertex_start
                pVertex_end_dummy = pFlowline.pVertex_end
                if pVertex_start_dummy == pVertex_start:
                    iFlag_found = 1
                    pVertex_start = pVertex_end_dummy
                    aFlowline_stay.append(i)
                    pVertex_outlet = pVertex_end_dummy
                    break

            if iFlag_found == 1:
                iFlag_done = 0
            else:
                iFlag_done = 1

    return pVertex_outlet


class MeshRTreeCache:
    """
    Cache class for storing mesh R-tree indices and associated data by filename

    This eliminates the need to rebuild mesh R-trees when processing multiple
    basins that use the same mesh file, providing significant performance improvements.
    """

    def __init__(self):
        self.cache = {}  # filename -> cache_data

    def get_mesh_data(self, sFilename_mesh):
        """
        Get cached mesh data if available

        Args:
            sFilename_mesh: Path to mesh file

        Returns:
            dict or None: Cached mesh data or None if not cached
        """
        return self.cache.get(sFilename_mesh)

    def store_mesh_data(self, sFilename_mesh, mesh_data):
        """
        Store mesh data in cache

        Args:
            sFilename_mesh: Path to mesh file
            mesh_data: Dictionary containing mesh R-tree index and metadata
        """
        self.cache[sFilename_mesh] = mesh_data
        print(f"Cached mesh data for {sFilename_mesh} ({len(mesh_data['mesh_features'])} cells)")

    def clear_cache(self):
        """Clear all cached mesh data"""
        self.cache.clear()
        print("Mesh R-tree cache cleared")

    def get_cache_info(self):
        """Get information about cached files"""
        return {filename: len(data['mesh_features']) for filename, data in self.cache.items()}


# Global cache instance
_mesh_rtree_cache = MeshRTreeCache()


def intersect_flowline_with_mesh_cached(iMesh_type_in, sFilename_mesh_in, sFilename_flowline_in,
                                       sFilename_output_in, mesh_cache=None, optimization_method='auto'):
    """
    Cached version of intersect_flowline_with_mesh with mesh R-tree caching

    Reuses mesh R-tree indices across multiple calls with the same mesh file,
    providing significant performance improvements for multi-basin processing.

    Args:
        iMesh_type_in: Mesh type identifier
        sFilename_mesh_in: Path to mesh GeoJSON file
        sFilename_flowline_in: Path to flowline GeoJSON file
        sFilename_output_in: Path to output GeoJSON file
        mesh_cache: Optional MeshRTreeCache instance (uses global cache if None)
        optimization_method: 'auto', 'mesh_rtree', or 'flowline_rtree'

    Returns:
        tuple: (aCell_intersect, aFlowline_intersect_all, pVertex_outlet)
    """

    if not (os.path.exists(sFilename_mesh_in) and os.path.exists(sFilename_flowline_in)):
        print('The input file does not exist')
        return

    if os.path.exists(sFilename_output_in):
        os.remove(sFilename_output_in)

    # Use provided cache or global cache
    if mesh_cache is None:
        mesh_cache = _mesh_rtree_cache

    # Check if mesh data is cached
    cached_mesh_data = mesh_cache.get_mesh_data(sFilename_mesh_in)

    if cached_mesh_data is not None:
        print(f"Using cached mesh data for {sFilename_mesh_in}")
        # Use cached mesh data with more aggressive mesh R-tree preference
        return _intersect_cached_mesh_rtree_method(
            iMesh_type_in, cached_mesh_data, sFilename_flowline_in,
            sFilename_output_in, optimization_method='auto_cached'
        )
    else:
        # First time with this mesh - cache it for future use
        print(f"Caching mesh data for {sFilename_mesh_in}")
        return _cache_and_intersect(
            iMesh_type_in, sFilename_mesh_in, sFilename_flowline_in,
            sFilename_output_in, mesh_cache, optimization_method
        )


def _intersect_cached_mesh_rtree_method(iMesh_type_in, cached_mesh_data, sFilename_flowline_in,
                                       sFilename_output_in, optimization_method='auto_cached'):
    """
    Intersection method using pre-cached mesh R-tree data

    Args:
        iMesh_type_in: Mesh type identifier
        cached_mesh_data: Pre-built mesh cache data
        sFilename_flowline_in: Path to flowline GeoJSON file
        sFilename_output_in: Path to output GeoJSON file
        optimization_method: Method selection (auto_cached uses lower threshold)

    Returns:
        tuple: (aCell_intersect, aFlowline_intersect_all, pVertex_outlet)
    """
    import time
    start_time = time.time()

    pDriver_geojson = ogr.GetDriverByName("GeoJSON")

    # Open flowline dataset
    pDataset_flowline = pDriver_geojson.Open(sFilename_flowline_in, 0)
    pLayer_flowline = pDataset_flowline.GetLayer(0)
    pSpatial_reference_flowline = pLayer_flowline.GetSpatialRef()
    nfeature_flowline = pLayer_flowline.GetFeatureCount()

    # Get cached mesh data
    index_mesh = cached_mesh_data['rtree_index']
    mesh_features = cached_mesh_data['mesh_features']
    pSpatial_reference_mesh = cached_mesh_data['spatial_reference']
    nfeature_mesh = len(mesh_features)

    print(f"Cached optimization analysis: {nfeature_mesh} mesh cells, {nfeature_flowline} flowlines")

    # Automatic method selection with lower threshold for cached mesh
    use_mesh_rtree = True
    if optimization_method == 'auto_cached':
        # More aggressive mesh R-tree usage since mesh is pre-cached (1.2x vs 2x threshold)
        if nfeature_mesh > nfeature_flowline * 1.2:
            use_mesh_rtree = True
            print("Using CACHED MESH R-tree optimization (mesh data pre-cached)")
        else:
            use_mesh_rtree = False
            print("Using FLOWLINE R-tree optimization (even with cached mesh)")
    elif optimization_method == 'mesh_rtree':
        use_mesh_rtree = True
        print("Using CACHED MESH R-tree optimization (forced)")
    else:
        use_mesh_rtree = False
        print("Using FLOWLINE R-tree optimization (forced)")

    # Handle coordinate transformations
    comparison = pSpatial_reference_mesh.IsSame(pSpatial_reference_flowline)
    if comparison != 1:
        iFlag_transform = 1
        transform = osr.CoordinateTransformation(pSpatial_reference_mesh, pSpatial_reference_flowline)
    else:
        iFlag_transform = 0

    # Setup output
    pDataset_out = pDriver_geojson.CreateDataSource(sFilename_output_in)
    pLayerOut = pDataset_out.CreateLayer('flowline', pSpatial_reference_flowline, ogr.wkbMultiLineString)
    pLayerOut.CreateField(ogr.FieldDefn('lineid', ogr.OFTInteger64))
    pLayerOut.CreateField(ogr.FieldDefn('stream_segment', ogr.OFTInteger))
    pLayerOut.CreateField(ogr.FieldDefn('stream_order', ogr.OFTInteger))
    pLayerDefn = pLayerOut.GetLayerDefn()
    pFeatureOut = ogr.Feature(pLayerDefn)

    lFlowlineID = 1
    aFlowline_intersect_all = list()
    aCell_intersect = list()

    if use_mesh_rtree:
        # Use cached mesh R-tree method
        processed_cells = {}

        # Loop through flowlines
        for i in range(nfeature_flowline):
            pFeature_flowline = pLayer_flowline.GetFeature(i)
            pGeometry_flowline = pFeature_flowline.GetGeometryRef()
            left, right, bottom, top = pGeometry_flowline.GetEnvelope()
            pBound_flowline = (left, bottom, right, top)

            iStream_segment = pFeature_flowline.GetField("stream_segment")
            iStream_order = pFeature_flowline.GetField("stream_order")

            # Find potentially intersecting mesh cells using cached R-tree
            aIntersect_mesh = list(index_mesh.intersection(pBound_flowline))

            for mesh_idx in aIntersect_mesh:
                data = mesh_features[mesh_idx]
                pGeometry_mesh = data['geometry']

                # Check actual intersection
                iFlag_intersect = pGeometry_flowline.Intersects(pGeometry_mesh)
                if iFlag_intersect:
                    # Get or create cell object
                    lCellID = data['cell_id']
                    if lCellID not in processed_cells:
                        pGeometrytype_mesh = pGeometry_mesh.GetGeometryName()
                        if pGeometrytype_mesh == 'POLYGON':
                            pCell = convert_gcs_coordinates_to_cell(iMesh_type_in, data['longitude'],
                                                                  data['latitude'], data['coords_gcs'])
                            pCell.lCellID = lCellID
                            pCell.dArea = data['area']
                            pCell.dLength = pCell.calculate_edge_length()
                            pCell.dLength_flowline = pCell.dLength
                            pCell.aFlowline = list()
                            pCell.nFlowline = 0
                            pCell.iFlag_intersected = 0
                            processed_cells[lCellID] = pCell
                            aCell_intersect.append(pCell)
                    else:
                        pCell = processed_cells[lCellID]

                    # Process intersection geometry
                    pGeometry_intersect = pGeometry_flowline.Intersection(pGeometry_mesh)
                    pGeometrytype_intersect = pGeometry_intersect.GetGeometryName()

                    if pGeometrytype_intersect == 'LINESTRING':
                        lFlowlineID = _process_linestring_intersection(
                            pGeometry_intersect, pFeatureOut, pLayerOut, lFlowlineID,
                            iStream_segment, iStream_order, pCell, aFlowline_intersect_all
                        )
                    elif pGeometrytype_intersect == 'MULTILINESTRING':
                        lFlowlineID = _process_multilinestring_intersection(
                            pGeometry_intersect, pFeatureOut, pLayerOut, lFlowlineID,
                            iStream_segment, iStream_order, pCell, aFlowline_intersect_all
                        )

        # Update cell properties
        for pCell in aCell_intersect:
            pCell.nFlowline = len(pCell.aFlowline)
            if pCell.nFlowline > 0:
                pCell.iFlag_intersected = 1
                pCell.dLength_flowline = 0.0
                for pFlowline in pCell.aFlowline:
                    if pFlowline.dLength > pCell.dLength_flowline:
                        pCell.dLength_flowline = pFlowline.dLength
    else:
        # Fallback to flowline R-tree method
        print("Note: Using non-cached flowline R-tree method")
        aCell_intersect, aFlowline_intersect_all, lFlowlineID = _intersect_flowline_rtree_cached(
            iMesh_type_in, mesh_features, pLayer_flowline, pLayerOut, pFeatureOut,
            nfeature_mesh, nfeature_flowline, lFlowlineID
        )

    # Quality control
    pVertex_outlet = _perform_quality_control(pLayer_flowline, aFlowline_intersect_all, nfeature_flowline)

    end_time = time.time()
    print(f"Cached intersection completed in {end_time - start_time:.2f} seconds")
    print(f"Found {len(aFlowline_intersect_all)} total intersections")

    return aCell_intersect, aFlowline_intersect_all, pVertex_outlet


def get_mesh_cache():
    """Get the global mesh cache instance"""
    return _mesh_rtree_cache


def clear_mesh_cache():
    """Clear the global mesh cache"""
    _mesh_rtree_cache.clear_cache()


def _cache_and_intersect(iMesh_type_in, sFilename_mesh_in, sFilename_flowline_in,
                        sFilename_output_in, mesh_cache, optimization_method):
    """
    Cache mesh data and perform intersection
    """
    pDriver_geojson = ogr.GetDriverByName("GeoJSON")

    # Open and cache mesh dataset
    pDataset_mesh = pDriver_geojson.Open(sFilename_mesh_in, 0)
    pLayer_mesh = pDataset_mesh.GetLayer(0)
    pSpatial_reference_mesh = pLayer_mesh.GetSpatialRef()
    nfeature_mesh = pLayer_mesh.GetFeatureCount()

    print("Building and caching mesh R-tree index...")

    # Build mesh R-tree and cache data
    index_mesh = RTreeindex()
    mesh_features = {}

    for j in range(nfeature_mesh):
        pFeature_mesh = pLayer_mesh.GetFeature(j)
        pGeometry_mesh = pFeature_mesh.GetGeometryRef()
        left_orig, right_orig, bottom, top = pGeometry_mesh.GetEnvelope()
        pBound = (left_orig, bottom, right_orig, top)
        index_mesh.insert(j, pBound)

        # Cache mesh feature data
        mesh_features[j] = {
            'geometry': pGeometry_mesh.Clone(),  # Clone to avoid reference issues
            'coords_gcs': get_geometry_coordinates(pGeometry_mesh),
            'cell_id': pFeature_mesh.GetField("cellid"),
            'longitude': pFeature_mesh.GetField("longitude"),
            'latitude': pFeature_mesh.GetField("latitude"),
            'area': pFeature_mesh.GetField("area")
        }

    # Store in cache
    cached_data = {
        'rtree_index': index_mesh,
        'mesh_features': mesh_features,
        'spatial_reference': pSpatial_reference_mesh,
        'feature_count': nfeature_mesh
    }
    mesh_cache.store_mesh_data(sFilename_mesh_in, cached_data)

    # Now perform intersection using cached data
    return _intersect_cached_mesh_rtree_method(
        iMesh_type_in, cached_data, sFilename_flowline_in,
        sFilename_output_in, optimization_method='auto_cached'
    )


def _intersect_flowline_rtree_cached(iMesh_type_in, mesh_features, pLayer_flowline, pLayerOut, pFeatureOut,
                                   nfeature_mesh, nfeature_flowline, lFlowlineID):
    """
    Flowline R-tree method using cached mesh features (fallback method)
    """
    aCell_intersect = list()
    aFlowline_intersect_all = list()

    print("Building flowline R-tree index...")

    # Build R-tree for flowlines
    index_flowline = RTreeindex()
    for i in range(nfeature_flowline):
        pFeature_flowline = pLayer_flowline.GetFeature(i)
        pGeometry_flowline = pFeature_flowline.GetGeometryRef()
        left, right, bottom, top = pGeometry_flowline.GetEnvelope()
        pBound = (left, bottom, right, top)
        index_flowline.insert(i, pBound)

    print("Processing mesh intersections...")

    # Loop through cached mesh features
    for j, data in mesh_features.items():
        pGeometry_mesh = data['geometry']
        pGeometrytype_mesh = pGeometry_mesh.GetGeometryName()

        if pGeometrytype_mesh == 'POLYGON':
            pCell = convert_gcs_coordinates_to_cell(iMesh_type_in, data['longitude'], data['latitude'], data['coords_gcs'])
            pCell.lCellID = data['cell_id']
            pCell.dArea = data['area']
            pCell.dLength = pCell.calculate_edge_length()
            pCell.dLength_flowline = pCell.dLength
            aFlowline_intersect = list()
            iFlag_intersected = 0

            left_orig, right_orig, bottom, top = pGeometry_mesh.GetEnvelope()
            pBound = (left_orig, bottom, right_orig, top)
            aIntersect = list(index_flowline.intersection(pBound))

            for k in aIntersect:
                pFeature_flowline = pLayer_flowline.GetFeature(k)
                pGeometry_flowline = pFeature_flowline.GetGeometryRef()
                iFlag_intersect = pGeometry_flowline.Intersects(pGeometry_mesh)

                if iFlag_intersect:
                    iFlag_intersected = 1
                    pGeometry_intersect = pGeometry_flowline.Intersection(pGeometry_mesh)
                    pGeometrytype_intersect = pGeometry_intersect.GetGeometryName()
                    iStream_segment = pFeature_flowline.GetField("stream_segment")
                    iStream_order = pFeature_flowline.GetField("stream_order")

                    if pGeometrytype_intersect == 'LINESTRING':
                        lFlowlineID = _process_linestring_intersection(
                            pGeometry_intersect, pFeatureOut, pLayerOut, lFlowlineID,
                            iStream_segment, iStream_order, pCell, aFlowline_intersect_all
                        )
                        aFlowline_intersect.append(pCell.aFlowline[-1])
                    elif pGeometrytype_intersect == 'MULTILINESTRING':
                        start_count = len(pCell.aFlowline)
                        lFlowlineID = _process_multilinestring_intersection(
                            pGeometry_intersect, pFeatureOut, pLayerOut, lFlowlineID,
                            iStream_segment, iStream_order, pCell, aFlowline_intersect_all
                        )
                        aFlowline_intersect.extend(pCell.aFlowline[start_count:])

            # Update cell properties
            pCell.aFlowline = aFlowline_intersect
            pCell.nFlowline = len(aFlowline_intersect)
            if iFlag_intersected == 1:
                pCell.iFlag_intersected = 1
                pCell.dLength_flowline = 0.0
                for pFlowline in pCell.aFlowline:
                    if pFlowline.dLength > pCell.dLength_flowline:
                        pCell.dLength_flowline = pFlowline.dLength
                aCell_intersect.append(pCell)
            else:
                pCell.iFlag_intersected = 0

    return aCell_intersect, aFlowline_intersect_all, lFlowlineID