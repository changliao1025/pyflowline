"""Convert vector flowline files to standardized vector formats.

This module provides functionality to convert flowlines from any GDAL-supported
vector format (Shapefile, GeoPackage, GeoJSON, etc.) into a standardized vector file
with WGS84 geographic coordinates and line IDs.
"""

from __future__ import annotations

import logging
import os

from pyearthriver.classes.flowline import pyflowline
from pyearthriver.formats.export_flowline import export_flowline_to_vector
from pyearthriver.formats.load_flowline import read_flowlines_from_vector

logger = logging.getLogger(__name__)

__all__ = ["convert_flowline_to_vector", "convert_flowline_to_geojson"]


def convert_flowline_to_vector(
    sFilename_in: str,
    sFilename_vector_out: str,
    id_field: str | None = None,
) -> list[pyflowline]:
    """Convert any vector flowline file to a standardized vector file with line IDs.

    Reads flowlines from any supported vector format (ESRI Shapefile, GeoPackage,
    GeoJSON, etc.), transforms coordinates to WGS84 (EPSG:4326), populates flowline
    IDs and stream attributes, and exports to the requested vector format.

    Args:
        sFilename_in: Path to the input vector file.
        sFilename_vector_out: Path to the output vector file (.geojson, .shp, .gpkg).
        id_field: Optional field name to use as flowline ID.

    Returns:
        List of pyflowline objects that were exported.

    Raises:
        FileNotFoundError: If the input vector file does not exist.
        RuntimeError: If vector reading or exporting fails.
    """
    if not os.path.isfile(sFilename_in):
        msg = f"This input file does not exist: {sFilename_in}"
        logger.error(msg)
        raise FileNotFoundError(msg)

    # Ensure output directory exists
    sWorkspace_out = os.path.dirname(os.path.abspath(sFilename_vector_out))
    if sWorkspace_out and not os.path.exists(sWorkspace_out):
        os.makedirs(sWorkspace_out, exist_ok=True)

    # Ingest flowlines and spatial reference
    aFlowline_basin, pSpatial_ref = read_flowlines_from_vector(
        sFilename_in,
        id_field=id_field,
        return_spatial_ref=True,
    )

    # Collect populated attributes to preserve in output vector file
    aAttribute_field = []
    aAttribute_dtype = []
    aAttribute_data = []

    has_stream_order = any(
        f.iStream_order is not None and f.iStream_order > 0 for f in aFlowline_basin
    )
    has_stream_segment = any(
        f.iStream_segment is not None and f.iStream_segment > 0 for f in aFlowline_basin
    )
    has_drainage_area = any(
        f.dDrainage_area is not None and f.dDrainage_area > 0 for f in aFlowline_basin
    )
    has_nhdplus_id = any(f.lNHDPlusID is not None and f.lNHDPlusID > 0 for f in aFlowline_basin)

    if has_stream_order:
        aAttribute_field.append("stream_order")
        aAttribute_dtype.append("int")
        aAttribute_data.append(
            [f.iStream_order if f.iStream_order is not None else 0 for f in aFlowline_basin]
        )

    if has_stream_segment:
        aAttribute_field.append("stream_segment")
        aAttribute_dtype.append("int")
        aAttribute_data.append(
            [f.iStream_segment if f.iStream_segment is not None else 0 for f in aFlowline_basin]
        )

    if has_drainage_area:
        aAttribute_field.append("drainage_area")
        aAttribute_dtype.append("float")
        aAttribute_data.append(
            [f.dDrainage_area if f.dDrainage_area is not None else 0.0 for f in aFlowline_basin]
        )

    if has_nhdplus_id:
        aAttribute_field.append("nhdplus_id")
        aAttribute_dtype.append("int")
        aAttribute_data.append(
            [f.lNHDPlusID if f.lNHDPlusID is not None else 0 for f in aFlowline_basin]
        )

    # Export using pyearthriver's general export_flowline_to_vector
    export_flowline_to_vector(
        aFlowline_basin,
        sFilename_vector_out,
        pSpatial_reference_in=pSpatial_ref,
        aAttribute_field=aAttribute_field if aAttribute_field else None,
        aAttribute_dtype=aAttribute_dtype if aAttribute_dtype else None,
        aAttribute_data=aAttribute_data if aAttribute_data else None,
    )

    logger.info(
        "Successfully converted %d flowline reaches from %s to %s",
        len(aFlowline_basin),
        sFilename_in,
        sFilename_vector_out,
    )
    return aFlowline_basin


def convert_flowline_to_geojson(
    sFilename_in: str,
    sFilename_geojson_out: str,
    id_field: str | None = None,
) -> list[pyflowline]:
    """Convert any vector flowline file to a standardized GeoJSON file with line IDs.

    Backward-compatible convenience wrapper for convert_flowline_to_vector.
    """
    return convert_flowline_to_vector(sFilename_in, sFilename_geojson_out, id_field=id_field)
