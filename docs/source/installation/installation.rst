############
Installation
############


********
Overview
********

This document provides the instruction to install the PyFlowline Python package.

Two different options are provided below.

************
Requirements
************

We recommend to use the Conda system to install the PyFlowline package.

Conda can be installed on Linux, MacOS, and Windows systems. 
Please refer to the conda website for details on how to install Conda: 
https://docs.conda.io/en/latest/

After Conda is available on your system, you can create a conda environment for your application of PyFlowline.
Then use Option A or B to install PyFlowline in the newly created environment.

==========
Option A
==========

In this option, you will use conda to install the released PyFlowline package, but not necessarily the latest version.
Conda will automatically install all the dependency packages.

    conda install -c conda-forge pyflowline


==========
Option B
==========

In this option, you have the opportunity to manually install the `nightly` version.

First, you need to clone the PyFlowline package from GitHub directly.

Navigate into the downloaded folder and manually install the package using:

    python setup.py install

The following dependency packages will be installed during the process.

* `numpy`
* `gdal`
* `netCDF4`

=============
Visualization
=============

PyFlowline only provides experimental support for visualization through the optional `matplotlib` and `cartopy` packages.

You need to manually speicify these packages during the installation process

    conda install -c conda-forge pyflowline matplotlib cartopy

or install manually after the installation of PyFlowline:

    conda install -c conda-forge matplotlib cartopy


