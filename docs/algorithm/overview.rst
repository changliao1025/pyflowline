* Overall model structure

The overall model workflow includes three components:

1. Preprocess vector river network
2. Generate or import meshes
3. Intersect mesh with flowline, and generate conceptual river network

* Connect disconnected flowline

If two flowlines are disconnected more than the threshold of floating data type error, this algorithm can be used to connect them by providing a nearest point and a searching radius. If a starting or ending vertex falls within the radius, a new flowline will be built based on the existing flowline and the provided point.


* Correct flow direction

Due to data quality issue, existing flowline may have incorrect flow direction, which leads to multiple downstream flow direction. The corresponding node connection matrix has rows with multiple 1s. This algorithm scans from the outlet node and searches reversely, once such a row was detected, the corresponding flow direction is flipped.


* Remove small river

To simplify river network, small river with length less than the user provided threshold is removed. This algorithm only applies to headwater and should be called multiple times to achieve desired performance.


* Remove braided loop

Braided loop occurs when a node has more than one downstream even after flow direction correction. This algorithm removes these loops by only keeping the first detected downstream of any node.


* Find critical vertex

The start and end vertices of a flowline define its type. 

1. If the start vertex has no upstream, this flowline is a headwater.
2. If the start or end vertex has only one upstream or downstream, it is a middle flowline and can be merged with others. 
3. If a start vertex has more than one upstream vertices, it is a river confluence.

The vertex type information is used to merge segmented flowlines.



* Merge flowline

This algorithm merge flowlines so there is only 2 types of flowlines:

1. headwater

2. flowline between confluence

If multiple flowlines are within the same confluence bound, they are merged as one.


* Mesh generation

This algorithm generates different types meshes based on the spatial extent and provided resolution.  

* Mesh and flowline intersection

This algorithm intersects any provided mesh with preprocessed flowline.

* Topology re-construction

Based on the intersection results, this algorithm build the upstream-downstream relationship using the shared flowline vertices.

* Topology simplification

This algorithm simpify the topology information for several unusual scenarios. For example, if a flowline leaves and re-enter the same mesh cell through the same edge, this creates a loop in topology and will be simplified. 
