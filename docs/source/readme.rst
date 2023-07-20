#####################
What is PyFlowline?
#####################

*********
Overview
*********

PyFlowline is a mesh-independent river network generator for hydrologic models.

Traditionally, the river network is represented by the vector-based river network. These vector-based river networks usually cannot be directly used by hydrologic models, which use the mesh as the spatial discretization. 

In most cases, the vector-based river network datasets are used in a so-called stream-burning process to convert the vector-based river network to the mesh-based river network. However, all the existing stream-burning methods only support the structured meshes and there are also some other limitations.

To close this gap, PyFlowline was developed using a mesh-independent approach. At its core, PyFlowline uses the intersection between the vector river network and mesh to reconstruct the conceptual river network. It also addresses several limitations of the existing stream-burning methods.

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
Important notice
*****************

1. PyFlowline is designed to run at regional to global scale, so all the datasets use the geographic coordinate system (GCS) with the WGS84 datum. See more details at https://pyflowline.readthedocs.io/en/latest/data/data.html

2. Visualization of the PyFlowline outputs is only experimental. This feature is not fully developed yet. There is ongoing effort to use the `PyEarth` python package to provide this feature. 
