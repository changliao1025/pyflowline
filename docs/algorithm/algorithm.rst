#########
Algorithm
#########


*************************
Overview
*************************

Pyflowline implemented a list of algorithms to generate conceptual river networks through the following steps.

1. Simplify vector river network
2. Generate or import meshes
3. Intersect mesh with flowline, and generate conceptual river network

*************************
Flowline simplification
*************************


==============================
Dam associate flowline burning
==============================

==============================
Flowline vertex extraction
==============================

The start and end vertices of a flowline define its type. 

1. If the start vertex has no upstream, this flowline is a headwater.
2. If the start or end vertex has only one upstream or downstream, it is a middle flowline and can be merged with others. 
3. If a start vertex has more than one upstream vertices, it is a river confluence.

==============================
Split flowline
==============================

==============================
Flow direction correction
==============================

Due to data quality issues, the existing flowline may have incorrect flow directions, which lead to multiple downstream flow directions. 
The corresponding node connection matrix has rows with multiple 1s. This algorithm scans from the outlet node and searches reversely; once such a row is detected, the corresponding flow direction is reversed.

==============================
Remove small river
==============================

To simplify the river networks, small rivers with lengths less than the user-provided threshold are removed. This algorithm only applies to headwater and should be called multiple times to achieve desired performance.


==============================
Remove braided flowlines
==============================

A braided loop occurs when a node has more than one downstream, even after the flow direction correction. This algorithm removes these loops by only keeping the first detected downstream of any node.


==============================
Flowline confluence extraction
==============================


==============================
Merge flowline
==============================
This algorithm merges flowlines, so there are only two types of flowlines:

1. headwaters

2. flowline between the confluences

If there are multiple flowlines within the same confluence bound, they are merged as one.

==============================
Flowline confluence definition
==============================


==============================
Stream segment index
==============================


==============================
Stream segment order
==============================

==============================
Split flowline by length
==============================

*************************
Mesh generation
*************************

==============================
Structured mesh
==============================

------------------
Latitude-longitude
------------------

------------------
Projected
------------------

------------------
Hexagon
------------------

==============================
Unstructured mesh
==============================

------------------
MPAS
------------------

------------------
TIN
------------------

*******************************************
Topological relationship reconstruction
*******************************************

==============================
Mesh and flowline intersection
==============================

==============================
Remove returning flowline
==============================

This algorithm simplifies the topology information for several unusual scenarios. For example, if a flowline leaves and reenters the same mesh cell through the same edge, this creates a loop in topology and will be simplified. 


==============================
Split flowline to edge
==============================

=======================================
Topological relationship reconstruction
=======================================