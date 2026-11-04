"""Convert vector boundary files to standardized vector formats.

Reuses core vector operations from pyearth (convert_vector_format,
get_vector_driver_from_filename) to support any GDAL-compatible vector format
(Shapefile, GeoPackage, GeoJSON, etc.) rather than only GeoJSON.
"""

from __future__ import annotations

import logging
import os

from osgeo import gdal, ogr, osr
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_filename
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.toolbox.conversion.convert_vector_format import convert_vector_format

logger = logging.getLogger(__name__)

__all__ = ["convert_boundary_to_vector", "convert_boundary_to_geojson"]


def convert_boundary_to_vector(
    sFilename_in: str,
    sFilename_vector_out: str,
    iFlag_largest_polygon_in: int = 0,
    target_epsg: int | None = None,
) -> bool:
    """Convert any vector boundary dataset to a target vector format.

    Supports Shapefile (.shp), GeoJSON (.geojson, .json), GeoPackage (.gpkg),
    and all other vector formats supported by GDAL and pyearth.

    Args:
        sFilename_in: Path to input vector file.
        sFilename_vector_out: Path to output vector file. Format is inferred
            automatically from extension.
        iFlag_largest_polygon_in: If 1, extracts only the largest polygon from
            Polygon / MultiPolygon features. If 0 (default), converts the full dataset.
        target_epsg: Optional EPSG code for coordinate transformation. For GeoJSON,
            defaults to EPSG:4326.

    Returns:
        bool: True on success, False on failure.
    """
    if not os.path.isfile(sFilename_in):
        logger.error("This input file does not exist: %s", sFilename_in)
        return False

    sWorkspace_out = os.path.dirname(os.path.abspath(sFilename_vector_out))
    if sWorkspace_out and not os.path.exists(sWorkspace_out):
        os.makedirs(sWorkspace_out, exist_ok=True)

    if iFlag_largest_polygon_in == 1:
        pDataset_boundary = gdal.OpenEx(sFilename_in, gdal.OF_VECTOR | gdal.GA_ReadOnly)
        if pDataset_boundary is None:
            logger.error("Could not open input vector file: %s", sFilename_in)
            return False

        pLayer_boundary = pDataset_boundary.GetLayer(0)
        pSpatial_ref = pLayer_boundary.GetSpatialRef()
        if pSpatial_ref is None:
            pSpatial_ref = osr.SpatialReference()
            pSpatial_ref.ImportFromEPSG(4326)
        pSpatial_ref.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

        dArea_max = -1.0
        pGeometry_out = None

        pFeature_boundary = pLayer_boundary.GetNextFeature()
        while pFeature_boundary:
            pGeometry_tmp = pFeature_boundary.GetGeometryRef()
            if pGeometry_tmp is not None:
                pGeometrytype_boundary = pGeometry_tmp.GetGeometryName().upper()
                if pGeometrytype_boundary == "POLYGON":
                    aCoords_gcs = get_geometry_coordinates(pGeometry_tmp)
                    if aCoords_gcs is not None and len(aCoords_gcs) > 0:
                        aLon = [coord[0] for coord in aCoords_gcs]
                        aLat = [coord[1] for coord in aCoords_gcs]
                        dArea = calculate_polygon_area(aLon, aLat)
                    else:
                        dArea = pGeometry_tmp.GetArea()
                    if dArea > dArea_max:
                        dArea_max = dArea
                        pGeometry_out = pGeometry_tmp.Clone()
                elif pGeometrytype_boundary == "MULTIPOLYGON":
                    nPolygon = pGeometry_tmp.GetGeometryCount()
                    for i in range(nPolygon):
                        pGeometry = pGeometry_tmp.GetGeometryRef(i)
                        aCoords_gcs = get_geometry_coordinates(pGeometry)
                        if aCoords_gcs is not None and len(aCoords_gcs) > 0:
                            aLon = [coord[0] for coord in aCoords_gcs]
                            aLat = [coord[1] for coord in aCoords_gcs]
                            dArea = calculate_polygon_area(aLon, aLat)
                        else:
                            dArea = pGeometry.GetArea()
                        if dArea > dArea_max:
                            dArea_max = dArea
                            pGeometry_out = pGeometry.Clone()
                else:
                    logger.warning("Unsupported geometry type: %s", pGeometrytype_boundary)

            pFeature_boundary = pLayer_boundary.GetNextFeature()

        pDataset_boundary = None

        if pGeometry_out is None:
            logger.error("No valid polygon geometry found in %s", sFilename_in)
            return False

        if os.path.isfile(sFilename_vector_out):
            os.remove(sFilename_vector_out)

        pDriver = get_vector_driver_from_filename(sFilename_vector_out)
        pDataset_out = pDriver.CreateDataSource(sFilename_vector_out)
        pLayer_out = pDataset_out.CreateLayer(
            "boundary",
            geom_type=pGeometry_out.GetGeometryType(),
            srs=pSpatial_ref,
        )
        pLayer_out.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))
        pFeature_out = ogr.Feature(pLayer_out.GetLayerDefn())
        pFeature_out.SetGeometry(pGeometry_out)
        pFeature_out.SetField("id", 1)
        pLayer_out.CreateFeature(pFeature_out)
        pFeature_out = None
        pLayer_out = None
        pDataset_out = None

        # If target format is GeoJSON or target_epsg is requested, ensure CRS matches
        ext = os.path.splitext(sFilename_vector_out)[1].lower()
        if (ext in (".geojson", ".json") or target_epsg is not None) and (
            pSpatial_ref.GetAuthorityCode(None) != str(target_epsg or 4326)
        ):
            convert_vector_format(
                sFilename_vector_out,
                sFilename_vector_out,
                target_epsg=target_epsg or 4326,
            )
        return True

    # When largest polygon is not requested, convert directly using pyearth
    return convert_vector_format(
        sFilename_in,
        sFilename_vector_out,
        target_epsg=target_epsg,
    )


def convert_boundary_to_geojson(
    sFilename_in: str,
    sFilename_geojson_out: str,
    iFlag_largest_polygon_in: int = 0,
) -> bool:
    """Convert any vector boundary dataset to GeoJSON (EPSG:4326).

    Backward-compatible convenience wrapper for convert_boundary_to_vector.
    """
    return convert_boundary_to_vector(
        sFilename_in,
        sFilename_geojson_out,
        iFlag_largest_polygon_in=iFlag_largest_polygon_in,
        target_epsg=4326,
    )
