---
title: 'pyflowline: a mesh independent river network generator for hydrologic models'

tags:
  - Python
  - Hydrology
  - geographic information system

authors:
  - name: Chang Liao
    orcid: 0000-0002-7348-8858    
    affiliation: 1
  
  - name: Matt Cooper
    orcid: 0000-0002-0165-209X
    affiliation: 1

affiliations:
 - name: Atmospheric Sciences and Global Change, Pacific Northwest National Laboratory, Richland, WA, USA
   index: 1 
date: 20 May 2021

bibliography: paper.bib
---

# Summary

Spatial-distributed hydrologic models often requires high quality river network information as inputs. Preparation of river network information can be challenging, especially when unstructured meshes are used. To resolve this challenge, we develope a Python package to automate river network generation. It resolves issues including disconnected river, braided network. It also produces spatially-distributed cell-to-cell topology information, which can be used for precise river routing. The package supports various mesh types including traditional projected coordinate system, geographic coordinate system and unstructured mesh.

# Statement of need






# Algorithms and implementation

pyflowline takes advantage of Python language's object oriented programming (OOP) architect. The river network and all of its elements (i.e. segment, reach, confluence.) are described as objects. These objects are processed throughout the package when applicable. 

![The data model. \label{fig:oop}](https://github.com/changliao1025/pyflowline/blob/main/docs/figure/basic_element.png?raw=true)


A brief overview of the features provided by pyflowline is list in Table 1.

1. It uses JSON/GeoJSON as the internal exchange file format. Most objects within the package can be imported and exported during runtime;
2. It relies on the open source Geospatial Data Abstraction Library (GDAL) and a few other Python packages to process geospatial data type;
3. It supports traditional projected coordinate system (PCS), geographic coordinate system (GCS), hexagon and other unstructured mesh for spatially distributed simulations;
4. It provides several algorithms to automate the river network simplification, generalization, etc.

* Overall model structure

The overall model workflow includes three components:

1. Preprocess flowline
2. Generate mesh based on the spatial extend of flowline and desired resolution
3. Intersect mesh with flowline, produce topology information


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


In this example, the flowline AH is represented by a list of edges after the intersection. Each edge defines a topology relationship between two cells. For example, edge AB defines flow direction from cell a to cell c.
Because the topology contains loops: cb-bc-cb-bc, it is simplified as acd, which means the final flow direction is from cell a to c, then to d.


# Example results

A case study was performed for the Columbia River Basin (CRB).

Screenshot of before and after river networks at zoom-in regions are used to illustrate the effects of algorithms. Attribute tables is used when applicable.



* Remove loop


* Merge flowline


* Mesh generation



* Mesh and flowline intersection


# Acknowledgement

The model described in this repository was supported by:

* Earth System Model Development and Regional and Global Modeling and Analysis program areas of the U.S. Department of Energy, Office of Science, Office of Biological and Environmental Research as part of the multi-program, collaborative Integrated Coastal Modeling (ICoM) project.

* Earth System Model Development and Regional and Global Modeling and Analysis program areas of the U.S. Department of Energy, Office of Science, Office of Biological and Environmental Research as part of the multi-program, collaborative Interdisciplinary Research for Arctic Coastal Environments (InteRFACE) project.

A portion of this research was performed using PNNL Research Computing at Pacific Northwest National Laboratory. 

PNNL is operated for DOE by Battelle Memorial Institute under contract DE-AC05-76RL01830.

# References

1. Liao, C., Tesfa, T., Duan, Z., & Leung, L. R. (2020). Watershed delineation on a hexagonal mesh grid. Environmental Modelling & Software, 128, 104702.

2. Liao, C., Zhou, T., Xu, D., Barnes, R., Bisht, G., Li, H. Y., ... & Leung, L. R. (2022). Advances in hexagon mesh-based flow direction modeling. Advances in Water Resources, 104099.