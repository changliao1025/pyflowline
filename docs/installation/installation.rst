############
Installation
############


********
Overview
********

This document provides the instruction to install the PyFlowline Python package.

************
Requirements
************

If the Conda system is used to install the PyFlowline package, conda will automatically install the dependency packages.

The dependency packages can also be installed manually.

* numpy
* gdal
* netCDF4
* shapely
* cython (optional)
* cartopy (optional, for visualization)
* matplotlib (optional, for visualization)

***********
Instruction 
***********

In most cases, you can install PyFlowline through the Conda system:

    conda install -c conda-forge pyflowline

To enable the `cython` feature, the user needs to build the cython code under the `./pyflowline/algorithms/cython` directory following the cython user guide (https://cython.readthedocs.io/en/latest/index.html).

In most cases, if a C/C++ compiler is available on the system, run the following command:

    python setup.py build_ext --inplace

