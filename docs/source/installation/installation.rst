############
Installation
############


********
Overview
********

This document provides the instruction to install the PyFlowline Python package.

Two different options are provided for different applications.

************
Requirements
************

We recommend the users to use the Conda system to install the PyFlowline package.

Conda can be installed on Linux, MacOS, and Windows systems. 
Please refer to the conda website for details on how to install Conda: 
https://docs.conda.io/en/latest/

After Conda is available on your system, you can create a conda environment for your application.
Then you use use either the Option A or B to install PyFlowline in the newly created environment.

==========
Option A
==========

In this option, you will use conda to install the released PyFlowline package, but not necessaily the latest version.
Conda will automatically install all the dependency packages.

    conda install -c conda-forge pyflowline


==========
Option B
==========

In this option, you have the oppotunity to install the `nightly` version.

First, you need to clone the PyFlowline package from GitHub directly.


install the dependency packages using the `requirements.txt`.

The following dependency packages will be installed.

* `numpy`
* `gdal`
* `netCDF4`
* `shapely`

    python setup.py build_ext --inplace

