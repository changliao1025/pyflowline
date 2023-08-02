#####################
What is PyFlowline?
#####################

*********
Overview
*********

PyFlowline is a mesh-independent river network generator for hydrologic models.

River networks are landscape features typically represented using vector layers. However, most hydrologic models rely on regular grids to discretize the spatial domain and cannot directly ingest vector features into the model. As a result, hydrologic models usually implement a so-called stream-burning process to convert the vector-based river network into a mesh-based river network. 

However, all the existing stream-burning methods only support the structured meshes and there are also some other limitations. For example, existing stream-burning methods always treat the vector river networks as a binary mask and cannot describe the topology near river confluences and meanders. 

PyFlowline solves this issue by using a mesh-independent approach that intersects the vector river network and mesh to reconstruct the conceptual river network.

***********
Development
***********

PyFlowline is developed in an open-source, public repository hosted on Github: 
https://github.com/changliao1025/pyflowline

*********
Objective
*********

All the existing river network representation methods (except vector-based) only support the structured rectangle meshes.
As a result, if a spatially-distributed hydrologic model uses the unstructured mesh as the spatial discretization, there is no way to represent the river network.

To close this gap, PyFlowline was developed using a mesh-independent approach. At its core, PyFlowline uses the intersection between the vector river network and mesh to reconstruct the conceptual river network.


*****************
Target audience
*****************

PyFlowline is an advanced modeling tool for hydrologists and hydrologic modelers. 
Users of PyFlowline should be familiar with basic concepts in Geographic Information System (GIS), including vector and raster data, coordinate systems, and projections.

*****************
Important notice
*****************

1. PyFlowline is designed to run at regional to global scale, so all the datasets use the geographic coordinate system (GCS) with the WGS84 datum. See more details at https://pyflowline.readthedocs.io/en/latest/data/data.html

2. Visualization of the PyFlowline outputs is only experimental. This feature is not fully developed yet. There is ongoing effort to use the `PyEarth` python package to provide this feature. 
