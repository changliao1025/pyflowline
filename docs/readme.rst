#####################
What is PyFlowline?
#####################

*********
Overview
*********

PyFlowline is a mesh-independent river network generator for hydrologic models.

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






