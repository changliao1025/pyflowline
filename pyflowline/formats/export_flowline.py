"""Export flowline objects and topological link networks to vector and JSON formats."""

from __future__ import annotations

import json
from typing import Any

from osgeo import osr
from pyearthmesh.classes.link import pycelllink
from pyearthriver.classes.edge import pyedge
from pyearthriver.formats.export_flowline import export_flowline_to_vector

__all__ = [
    "export_flowline_to_vector",
    "export_flowline_to_geojson",
    "export_flowline_to_shapefile",
    "export_flowline_info_to_json",
]


def export_flowline_to_geojson(
    aFlowline_in: list,
    sFilename_json_in: str,
    iFlag_projected_in: int | None = None,
    pSpatial_reference_in: osr.SpatialReference | str | None = None,
    aAttribute_field: list[str] | None = None,
    aAttribute_data: list[list[Any]] | None = None,
    aAttribute_dtype: list[str] | None = None,
) -> None:
    """Export flowline objects to GeoJSON or any vector format."""
    return export_flowline_to_vector(
        aFlowline_in=aFlowline_in,
        sFilename_vector_out=sFilename_json_in,
        iFlag_projected_in=iFlag_projected_in,
        pSpatial_reference_in=pSpatial_reference_in,
        aAttribute_field=aAttribute_field,
        aAttribute_data=aAttribute_data,
        aAttribute_dtype=aAttribute_dtype,
    )


def export_flowline_to_shapefile(
    iFlag_projected_in: int | None,
    aFlowline_in: list,
    pSpatial_reference_in: osr.SpatialReference | str | None,
    sFilename_shapefile_in: str,
    aAttribute_field: list[str] | None = None,
    aAttribute_data: list[list[Any]] | None = None,
    aAttribute_dtype: list[str] | None = None,
) -> None:
    """Export flowline objects to ESRI Shapefile or any vector format."""
    return export_flowline_to_vector(
        aFlowline_in=aFlowline_in,
        sFilename_vector_out=sFilename_shapefile_in,
        iFlag_projected_in=iFlag_projected_in,
        pSpatial_reference_in=pSpatial_reference_in,
        aAttribute_field=aAttribute_field,
        aAttribute_data=aAttribute_data,
        aAttribute_dtype=aAttribute_dtype,
    )


def export_flowline_info_to_json(aCell, aCell_intersect_in, aFlowline_in, sFilename_json_out):
    """Export the cell-flowline link network and topology to a JSON file."""
    ncell = len(aCell_intersect_in)
    nflowline = len(aFlowline_in)
    aLink = list()
    for i in range(1, nflowline + 1):
        pFlowline = aFlowline_in[i - 1]
        nEdge = pFlowline.nEdge
        for j in range(1, nEdge + 1):
            pEdge = pFlowline.aEdge[j - 1]
            pVertex_start = pEdge.pVertex_start
            pVertex_end = pEdge.pVertex_end
            pMpas_start = None
            pMpas_end = None
            for k in range(ncell):
                if aCell_intersect_in[k].pVertex_center == pVertex_start:
                    pMpas_start = aCell_intersect_in[k]
                if aCell_intersect_in[k].pVertex_center == pVertex_end:
                    pMpas_end = aCell_intersect_in[k]
            pEdge_link = pyedge(pVertex_start, pVertex_end)
            pLink = pycelllink(pMpas_start, pMpas_end, pEdge_link)
            aLink.append(pLink)

    with open(sFilename_json_out, "w", encoding="utf-8") as f:
        sJson = json.dumps([json.loads(ob.tojson()) for ob in aLink], indent=4)
        f.write(sJson)
