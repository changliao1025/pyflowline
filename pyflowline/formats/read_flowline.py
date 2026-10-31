"""Module for reading river flowlines from vector spatial formats.

Supports Shapefile (.shp), GeoJSON (.geojson, .json), GeoPackage (.gpkg),
and any other vector format supported by pyearth and GDAL/OGR.
"""

import os
from typing import List, Optional, Tuple

from pyearthriver.classes.flowline import pyflowline
from pyearthriver.formats.load_flowline import read_flowlines_from_vector

from pyflowline.formats.export_flowline import export_flowline_to_geojson


def read_flowline_vector(
    sFilename_vector_in: str,
    sFilename_out: Optional[str] = None,
) -> Tuple[List[pyflowline], Optional[str]]:
    """Read stream flowlines from any supported vector format.

    Args:
        sFilename_vector_in: File path to the input vector file (.shp, .geojson, .gpkg, etc.).
        sFilename_out: Optional output file path. If provided, flowlines with
            resolved 'lineid' will be exported to this path.

    Returns:
        Tuple containing:
            - List of pyflowline objects.
            - Spatial reference WKT string (or None).

    Raises:
        FileNotFoundError: If the input file does not exist.
        RuntimeError: If GDAL/OGR cannot read the vector file.
    """
    if not os.path.exists(sFilename_vector_in):
        raise FileNotFoundError(f"This vector file does not exist: {sFilename_vector_in}")

    aFlowline, pProjection_geojson = read_flowlines_from_vector(
        sFilename_vector_in, return_spatial_ref=True
    )

    if sFilename_out is not None:
        aFlowlineID = [fl.lFlowlineID for fl in aFlowline]
        export_flowline_to_geojson(
            aFlowline,
            sFilename_out,
            aAttribute_field=["lineid"],
            aAttribute_dtype=["int"],
            aAttribute_data=[aFlowlineID],
        )

    return aFlowline, pProjection_geojson
