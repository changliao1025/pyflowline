# Overview

# Dependency

* pystream use some functions from pyearth, which is a general purpose library for many Earth scient operation.


# Requirements

Pystream requires several input information for different operations. In general, these dataset are expected:

* DEM raster file: this file will be used to obtain the spatial extent of the study domain
* River flowline vector file: this file is used to guide the stream topology definition

# Steps

In most cases, you will run pystream in the following steps:
- pre-process river network
    * you need the approximate location of outlet
- generate mesh
    * you need to use the DEM as spatial extent
    * you can specify resolution
- intersect mesh with river flowline
    * you need outlet location again
    * both flowline from step 1 and mesh from step 2 are used to re-conconstruct the topology

