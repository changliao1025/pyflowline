"""Export vertex objects to vector formats (GeoJSON, Shapefile, GeoPackage, etc.)."""

from __future__ import annotations

from osgeo import osr
from pyearthriver.formats.export_vertex import export_vertex_to_vector

__all__ = [
    "export_vertex_to_vector",
    "export_vertex_to_geojson",
    "export_vertex_to_shapefile",
]


def export_vertex_to_geojson(
    aVertex_in: list,
    sFilename_json_in: str,
    iFlag_projected_in: int | None = None,
    pSpatial_reference_in: osr.SpatialReference | None = None,
    aAttribute_data: list[int] | None = None,
) -> None:
    """Export vertex objects to GeoJSON or any vector format."""
    return export_vertex_to_vector(
        aVertex_in=aVertex_in,
        sFilename_vector_out=sFilename_json_in,
        iFlag_projected_in=iFlag_projected_in,
        pSpatial_reference_in=pSpatial_reference_in,
        aAttribute_data=aAttribute_data,
    )


def export_vertex_to_shapefile(
    aVertex_in: list,
    sFilename_shapefile_out: str,
    pSpatial_reference_in: osr.SpatialReference | None,
    iFlag_projected_in: int | None,
    aAttribute_data: list[int] | None = None,
) -> None:
    """Export vertex objects to ESRI Shapefile or any vector format."""
    return export_vertex_to_vector(
        aVertex_in=aVertex_in,
        sFilename_vector_out=sFilename_shapefile_out,
        iFlag_projected_in=iFlag_projected_in,
        pSpatial_reference_in=pSpatial_reference_in,
        aAttribute_data=aAttribute_data,
    )
