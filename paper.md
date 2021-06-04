---
title: 'PyStream: A Python package for stream networks processing'

tags:
  - Python
  - Hydrology
  - geographic information system

authors:
  - name: Chang Liao
    orcid: 0000-0002-7348-8858    
    affiliation: 1

affiliations:
 - name: Atmospheric Sciences and Global Change, Pacific Northwest National Laboratory, Richland, WA, USA
   index: 1 
date: 20 May 2021

bibliography: paper.bib
---

# Summary

Computational hydrologic simulation often requires high quality river network information as inputs. Preparation of river network information can be challenging, especially if the hydrology model is based on unstructured mesh. In this study, we develope a Python package to automate river network preparation. It resolves issues including disconnected river, braided network. It also produces spatially-distributed cell-to-cell topology information, which can be used for precise river routing. The package supports various mesh types including traditional projected coordinate system, geographic coordinate system and unstructured mesh.

# Statement of need






# Algorithms and implementation

PyStream takes advantage of Python language's object oriented programming (OOP) architect. The river network and all of its elements (i.e. segment, reach, confluence.) are described as objects. These objects are processed throughout the package when applicable. 

![The data model. \label{fig:oop}](https://github.com/changliao1025/pystream/blob/main/pystream/figure/data_mode.png?raw=true)


A brief overview of the features provided by PyStream is list in Table 1.

1. It uses JSON/GeoJSON as the internal exchange file format. Most objects within the package can be imported and exported during runtime;
2. It relies on the open source Geospatial Data Abstraction Library (GDAL) and a few other Python packages to process geospatial data type;
3. It supports traditional projected coordinate system (PCS), geographic coordinate system (GCS), hexagon and other unstructured mesh for spatially distributed simulations;
4. It provides several algorithms to automate the river network simplification, generalization, etc.

* Overall model structure

The overall model workflow includes three components:

1. Preprocess flowline
2. Generate mesh based on the spatial extend of flowline and desired resolution
3. Intersect mesh with flowline, produce topology information

![The workflow of PyStream. \label{fig:workflow}](https://github.com/changliao1025/pystream/blob/main/pystream/figure/workflow.png?raw=true)

* Connect disconnected flowline

If two flowlines are disconnected more than the threshold of floating data type error, this algorithm can be used to connect them by providing a nearest point and a searching radius. If a starting or ending vertex falls within the radius, a new flowline will be built based on the existing flowline and the provided point.

![Disconnect flowline. \label{fig:disconnected}](https://github.com/changliao1025/pystream/blob/main/pystream/figure/disconnect_flowline.png?raw=true)

* Correct flow direction

Due to data quality, existing flowline may have incorrect flow direction, which leads to multiple downstream flow direction. The corresponding node connection matrix has rows with multiple 1s. This algorithm scans from the outlet node and searches reversely, once such a row was detected, the corresponding flow direction is flipped.

![Flow direction. \label{fig:direction}](https://github.com/changliao1025/pystream/blob/main/pystream/figure/flow_direction_matrix.png?raw=true)

* Remove small river

To simplify river network, small river with length less than the user provided threshold is removed. This algorithm only applies to headwater and should be called multiple times to achieve desired performance.

![Small river. \label{fig:small_river}](https://github.com/changliao1025/pystream/blob/main/pystream/figure/small_river.png?raw=true)

* Remove braided loop

Braided loop occurs when a node has more than one downstream even after flow direction correction. This algorithm removes these loops by only keeping the first detected downstream of any node.

![Remove loops. \label{fig:loops}](https://github.com/changliao1025/pystream/blob/main/pystream/figure/remove_loop_matrix.png?raw=true)

* Find critical vertex

The start and end vertices of a flowline define its type. 

1. If the start vertex has no upstream, this flowline is a headwater.
2. If the start or end vertex has only one upstream or downstream, it is a middle flowline and can be merged with others. 
3. If a start vertex has more than one upstream vertices, it is a river confluence.

The vertex type information is used to merge segmented flowlines.

![Vertex type. \label{fig:vertex}](https://github.com/changliao1025/pystream/blob/main/pystream/figure/find_vertex.png?raw=true)


* Merge flowline

This algorithm merge flowlines so there is only 2 types of flowlines:

1. headwater
2. flowline between confluence

![Merge flowline. \label{fig:merge}](https://github.com/changliao1025/pystream/blob/main/pystream/figure/merge_flowline.png?raw=true)

* Mesh generation

This algorithm generates different types meshes based on the spatial extent and provided resolution.  

* Mesh and flowline intersection

This algorithm intersects any provided mesh with preprocessed flowline.

* Topology re-construction

Based on the intersection results, this algorithm build the upstream-downstream relationship using the shared flowline vertices.

* Topology simplification

This algorithm simpify the topology information for several unusual scenarios. For example, if a flowline leaves and re-enter the same mesh cell through the same edge, this creates a loop in topology and will be simplified. 

![Topology simplification. \label{fig:topology_simplification}](https://github.com/changliao1025/pystream/blob/main/pystream/figure/simplification01.png?raw=true)




# Example results

A case study was performed for the Columbia River Basin (CRB).

Screenshot of before and after river networks at zoom-in regions are used to illustrate the effects of algorithms. Attribute tables is used when applicable.



* Remove loop

![Before remove loop. \label{fig:before_remove_loop}](https://github.com/changliao1025/pystream/blob/main/pystream/figure/before_loop.png?raw=true)

![After remove loop. \label{fig:after_remove_loop}](https://github.com/changliao1025/pystream/blob/main/pystream/figure/after_loop.png?raw=true)

* Merge flowline

![Before merge flowline. \label{fig:before_merge}](https://github.com/changliao1025/pystream/blob/main/pystream/figure/before_merge.png?raw=true)

![After merge flowline. \label{fig:after_merge}](https://github.com/changliao1025/pystream/blob/main/pystream/figure/after_merge.png?raw=true)

* Mesh generation

![Lan Lon mesh. \label{fig:lat_lon_mesh}](https://github.com/changliao1025/pystream/blob/main/pystream/figure/lat_lon.png?raw=true)

![Sqaure mesh. \label{fig:square_mesh}](https://github.com/changliao1025/pystream/blob/main/pystream/figure/square.png?raw=true)

![Hexagon mesh. \label{fig:hexagon_mesh}](https://github.com/changliao1025/pystream/blob/main/pystream/figure/hexagon.png?raw=true)

![Overlap mesh. \label{fig:meshed}](https://github.com/changliao1025/pystream/blob/main/pystream/figure/meshes.png?raw=true)


* Mesh and flowline intersection

![Lat Lon intersect. \label{fig:lat_lon_intersect}](https://github.com/changliao1025/pystream/blob/main/pystream/figure/lat_lon_intersect.png?raw=true)

![Square intersect. \label{fig:square_intersect}](https://github.com/changliao1025/pystream/blob/main/pystream/figure/square_intersect.png?raw=true)

![Hexagon intersect. \label{fig:hexagon_intersect}](https://github.com/changliao1025/pystream/blob/main/pystream/figure/hexagon_intersect.png?raw=true)

# Acknowledgement

The model described in this repository was supported by:

* Laboratory Directed Research and Development (LDRD) Program Quickstarter project at Pacific Northwest National Laboratory. 
* Earth System Model Development and Regional and Global Modeling and Analysis program areas of the U.S. Department of Energy, Office of Science, Office of Biological and Environmental Research as part of the multi-program, collaborative Integrated Coastal Modeling (ICoM) project.
* U.S. Department of Energy Office of Science Biological and Environmental Research through the Earth and Environmental System Modeling program as part of the Energy Exascale Earth System Model (E3SM) project. 

A portion of this research was performed using PNNL Research Computing at Pacific Northwest National Laboratory. 

PNNL is operated for DOE by Battelle Memorial Institute under contract DE-AC05-76RL01830.

# References

1. Liao, C., Tesfa, T., Duan, Z., & Leung, L. R. (2020). Watershed delineation on a hexagonal mesh grid. Environmental Modelling & Software, 128, 104702.
