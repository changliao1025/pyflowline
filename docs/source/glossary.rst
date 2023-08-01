########
Glossary
########


*****************
Structured mesh
*****************

In PyFlowline, structured mesh refers to meshes that have a repeating pattern or structure.

The following meshes are considered as structured:

1. Projected raster meshes (e.g. 100m by 100m)
2. GCS-based rectangle meshes  (e.g. 0.5 degree by 0.5 degree)
3. Hexagon meshes (e.g. 100m by edge)
4. DGGS meshes (e.g., DGGrid meshes)

*****************
Unstructured mesh
*****************

In PyFlowline, unstructured mesh refers to meshes that don't have a repeating pattern or structure and the cell size varies from cell to cell.

The following meshes are considered as unstructured:

1. Model for Prediction Across Scales (MPAS) meshes
2. Triangulated irregular network (TIN) meshes

************
Great circle
************

In mathematics, a great circle or orthodrome is the circular intersection of a sphere and a plane passing through the sphere's center point.

****
DGGS
****

A discrete global grid (DGG) is a mosaic that covers the entire Earth's surface. Mathematically it is a space partitioning: it consists of a set of non-empty regions that form a partition of the Earth's surface. In a usual grid-modeling strategy, to simplify position calculations, each region is represented by a point, abstracting the grid as a set of region-points. Each region or region-point in the grid is called a cell.

****
TIN
****

In computer graphics, a triangulated irregular network (TIN) is a representation of a continuous surface consisting entirely of triangular facets (a triangle mesh), used mainly as Discrete Global Grid in primary elevation modeling.

****
MPAS
****

Model for Prediction Across Scales.