import os
import numpy as np
from osgeo import ogr, osr
from rtree.index import Index as RTreeindex
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_cell
from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_flowline
from pyearth.gis.geometry.international_date_line_utility import check_cross_international_date_line_geometry


def intersect_flowline_with_mesh(iMesh_type_in, sFilename_mesh_in, sFilename_flowline_in, sFilename_output_in):
    """
    Intersection function using cached mesh R-tree method for IDL-aware processing of MULTIPOLYGON geometries

    Args:
        iMesh_type_in: Mesh type identifier (1=hexagon, 2=square, 3=latlon, 4=mpas, 5=dggrid, 6=triangular, 7=tin, 8=cubic_sphere)
        sFilename_mesh_in: Path to mesh GeoJSON file
        sFilename_flowline_in: Path to flowline GeoJSON file
        sFilename_output_in: Path to output GeoJSON file

    Returns:
        tuple: (aCell_intersect, aFlowline_intersect_all, pVertex_outlet) or (None, None, None) on error

    Raises:
        ValueError: If mesh type is invalid
        RuntimeError: If files cannot be opened
        FileNotFoundError: If input files don't exist
    """
    # Input validation
    if not isinstance(iMesh_type_in, int) or iMesh_type_in < 1 or iMesh_type_in > 8:
        raise ValueError(f"Invalid mesh type: {iMesh_type_in}. Must be integer between 1-8")

    if not os.path.exists(sFilename_mesh_in):
        raise FileNotFoundError(f"Mesh file does not exist: {sFilename_mesh_in}")

    if not os.path.exists(sFilename_flowline_in):
        raise FileNotFoundError(f"Flowline file does not exist: {sFilename_flowline_in}")

    # Validate file extensions
    if not (sFilename_mesh_in.lower().endswith('.geojson') or sFilename_mesh_in.lower().endswith('.json')):
        print(f"Warning: Mesh file may not be GeoJSON format: {sFilename_mesh_in}")

    if not (sFilename_flowline_in.lower().endswith('.geojson') or sFilename_flowline_in.lower().endswith('.json')):
        print(f"Warning: Flowline file may not be GeoJSON format: {sFilename_flowline_in}")

    if os.path.exists(sFilename_output_in):
        os.remove(sFilename_output_in)

    # Use global cache
    mesh_cache = _mesh_rtree_cache

    # Check if mesh data is cached
    cached_mesh_data = mesh_cache.get_mesh_data(sFilename_mesh_in)

    if cached_mesh_data is not None:
        print(f"Using cached mesh data for {sFilename_mesh_in}")
        return _intersect_cached_mesh_rtree_method(
            iMesh_type_in, cached_mesh_data, sFilename_flowline_in, sFilename_output_in
        )
    else:
        print(f"Caching mesh data for {sFilename_mesh_in}")
        return _cache_and_intersect(
            iMesh_type_in, sFilename_mesh_in, sFilename_flowline_in, sFilename_output_in, mesh_cache
        )


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


def _intersect_cached_mesh_rtree_method(iMesh_type_in, cached_mesh_data, sFilename_flowline_in,
                                       sFilename_output_in):
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
    if pDataset_flowline is None:
        raise RuntimeError(f"Could not open flowline dataset: {sFilename_flowline_in}")

    pDataset_out = None
    try:
        pLayer_flowline = pDataset_flowline.GetLayer(0)
        pSpatial_reference_flowline = pLayer_flowline.GetSpatialRef()
        nfeature_flowline = pLayer_flowline.GetFeatureCount()

        # Get cached mesh data
        index_mesh = cached_mesh_data['rtree_index']
        mesh_features = cached_mesh_data['mesh_features']
        pSpatial_reference_mesh = cached_mesh_data['spatial_reference']
        nfeature_mesh = len(mesh_features)

        print(f"Cached optimization analysis: {nfeature_mesh} mesh cells, {nfeature_flowline} flowlines")

        # Handle coordinate transformations
        comparison = pSpatial_reference_mesh.IsSame(pSpatial_reference_flowline)
        if comparison != 1:
            iFlag_transform = 1
            transform = osr.CoordinateTransformation(pSpatial_reference_flowline, pSpatial_reference_mesh)
        else:
            iFlag_transform = 0
            transform = None

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
                # Transform flowline geometry if needed before intersection check
                if iFlag_transform == 1:
                    pGeometry_flowline_transformed = pGeometry_flowline.Clone()
                    pGeometry_flowline_transformed.Transform(transform)
                    iFlag_intersect = pGeometry_flowline_transformed.Intersects(pGeometry_mesh)
                    # Use transformed geometry for intersection calculation
                    pGeometry_intersect = pGeometry_flowline_transformed.Intersection(pGeometry_mesh)
                else:
                    iFlag_intersect = pGeometry_flowline.Intersects(pGeometry_mesh)
                    pGeometry_intersect = pGeometry_flowline.Intersection(pGeometry_mesh)

                if iFlag_intersect:
                    # Get or create cell object - handle both POLYGON and MULTIPOLYGON parts
                    lCellID = data['cell_id']
                    if lCellID not in processed_cells:
                        pGeometrytype_mesh = pGeometry_mesh.GetGeometryName()
                        if pGeometrytype_mesh in ['POLYGON', 'MULTIPOLYGON']:
                            # For MULTIPOLYGON parts, use the original geometry coordinates if available
                            coords_to_use = data['coords_gcs']
                            # Ensure we're using the correct coordinates for cell creation
                            pCell = convert_gcs_coordinates_to_cell(iMesh_type_in, data['longitude'],
                                                                  data['latitude'], coords_to_use)
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

                    # Process intersection geometry (already calculated above)
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

        # Quality control
        pVertex_outlet = _perform_quality_control(pLayer_flowline, aFlowline_intersect_all, nfeature_flowline)

        end_time = time.time()
        print(f"Cached intersection completed in {end_time - start_time:.2f} seconds")
        print(f"Found {len(aFlowline_intersect_all)} total intersections")

        return aCell_intersect, aFlowline_intersect_all, pVertex_outlet

    finally:
        # Ensure proper cleanup of resources
        if pDataset_flowline is not None:
            pDataset_flowline = None
        if pDataset_out is not None:
            pDataset_out = None


def get_mesh_cache():
    """Get the global mesh cache instance"""
    return _mesh_rtree_cache


def clear_mesh_cache():
    """Clear the global mesh cache"""
    _mesh_rtree_cache.clear_cache()


def _cache_and_intersect(iMesh_type_in, sFilename_mesh_in, sFilename_flowline_in,
                          sFilename_output_in, mesh_cache):
    """
    Cache mesh data and perform intersection
    """
    pDriver_geojson = ogr.GetDriverByName("GeoJSON")

    # Open and cache mesh dataset
    pDataset_mesh = pDriver_geojson.Open(sFilename_mesh_in, 0)
    if pDataset_mesh is None:
        raise RuntimeError(f"Could not open mesh dataset: {sFilename_mesh_in}")

    try:
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
            sGeometryname = pGeometry_mesh.GetGeometryName()

            if sGeometryname == 'POLYGON':
                left_orig, right_orig, bottom, top = pGeometry_mesh.GetEnvelope()
                pBound = (left_orig, bottom, right_orig, top)
                index_mesh.insert(j, pBound)

                # Cache mesh feature data
                mesh_features[j] = {
                    'geometry': pGeometry_mesh.Clone(),
                    'coords_gcs': get_geometry_coordinates(pGeometry_mesh),
                    'cell_id': pFeature_mesh.GetField("cellid"),
                    'longitude': pFeature_mesh.GetField("longitude"),
                    'latitude': pFeature_mesh.GetField("latitude"),
                    'area': pFeature_mesh.GetField("area")
                }

            elif sGeometryname == 'MULTIPOLYGON':
                nPolygons = pGeometry_mesh.GetGeometryCount()

                for k in range(nPolygons):
                    pPolygon = pGeometry_mesh.GetGeometryRef(k)
                    left_orig, right_orig, bottom, top = pPolygon.GetEnvelope()
                    pBound = (left_orig, bottom, right_orig, top)

                    # Use integer-based unique ID for consistency
                    part_idx = nfeature_mesh + j * 1000 + k  # Create unique integer ID
                    index_mesh.insert(part_idx, pBound)

                    # Cache each polygon part separately
                    mesh_features[part_idx] = {
                        'geometry': pPolygon.Clone(),
                        'coords_gcs': get_geometry_coordinates(pPolygon),
                        'cell_id': pFeature_mesh.GetField("cellid"),
                        'longitude': pFeature_mesh.GetField("longitude"),
                        'latitude': pFeature_mesh.GetField("latitude"),
                        'area': pFeature_mesh.GetField("area"),
                        'original_index': j,  # Track original feature index
                        'part_number': k      # Track which part this is
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
            iMesh_type_in, cached_data, sFilename_flowline_in, sFilename_output_in
        )

    finally:
        # Ensure proper cleanup
        if pDataset_mesh is not None:
            pDataset_mesh = None

