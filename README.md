### PyFlowline

[![DOI](https://zenodo.org/badge/368338554.svg)](https://zenodo.org/badge/latestdoi/368338554)
[![Downloads](https://static.pepy.tech/badge/pyflowline)](https://pepy.tech/project/pyflowline)

PyFlowline: a mesh-independent river network generator for hydrologic models. 

### Quickstart

Please refer to the [quickstart documentation](https://pyflowline.readthedocs.io/en/latest/quickstart.html) for details on how to get started using the PyFlowline package.

PyFlowline is mesh independent, meaning you can apply it to both structured 

1. traditional rectangle projected mesh 
2. latitude-longitude 
3. hexagon
4. dggs ([dggrid](https://github.com/sahrk/DGGRID))

and unstructured mesh systems 

1. Model for Prediction Across Scales mesh ([MPAS](https://github.com/MPAS-Dev))
2. Triangulated Irregular Network (TIN) mesh

This package generates the mesh cell-based conceptual river networks using the following steps:

1. `Flowline simplification`: PyFlowline checks the vector dataset and corrects undesired flowlines, such as braided rivers.
2. `Mesh generation`: PyFlowline generates structured meshes (e.g., rectangle, hexagon) or imports user-provided unstructured meshes into the PyFlowline-compatible GEOJSON format.
3. `Topological relationship reconstruction`: PyFlowline reconstructs the topological relationship using the mesh and flowline intersections.

### Dependency

PyFlowline depends on the following packages

1. `numpy`
2. `gdal`
3. `netCDF4`

PyFlowline also has three optional dependency packages

1. `cython` for performance 
2. `matplotlib` for visualization
3. `cartopy` for visulization
4. `simplekml` for Google Earth KML support

### Installation

Please refer to the [official documentation](https://pyflowline.readthedocs.io/) for details on how to install the PyFlowline package.

### Application

We provide several examples in the `examples` folder to demonstrate the model capability. We also recommend starting with the `notebooks/mpas_example.ipynb` notebook, after following the Quickstart and Installation instructions.

### Acknowledgment

This work was supported by the Earth System Model Development program areas of the U.S. Department of Energy, Office of Science, Office of Biological and Environmental Research as part of the multi-program, collaborative Integrated Coastal Modeling (ICoM) project and the Interdisciplinary Research for Arctic Coastal Environments (InteRFACE) project.

This research was supported as part of the Next Generation Ecosystem Experiments-Tropics, funded by the U.S. Department of Energy, Office of Science, Office of Biological and Environmental Research at Pacific Northwest National Laboratory. The study was also partly supported by U.S. Department of Energy Office of Science Biological and Environmental Research through the Earth and Environmental System Modeling program as part of the Energy Exascale Earth System Model (E3SM) project.

### License

BSD 3-Clause License

Copyright © 2022, Battelle Memorial Institute

1. Battelle Memorial Institute (hereinafter Battelle) hereby grants permission to any person or entity lawfully obtaining a copy of this software and associated documentation files (hereinafter “the Software”) to redistribute and use the Software in source and binary forms, with or without modification. Such person or entity may use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software and may permit others to do so, subject to the following conditions:

* Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimers.

* Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

* Other than as used herein, neither the name Battelle Memorial Institute or Battelle may be used in any form whatsoever without the express written consent of Battelle.

2. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL BATTELLE OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

### References

Several publications describe the algorithms used in `PyFlowline` in detail. If you make use of `PyFlowline` in your work, please consider including a reference to the following:

* Liao. C. Cooper, M (2022) Pyflowline: a mesh-independent river network generator for hydrologic models. Zenodo.
https://doi.org/10.5281/zenodo.6407299

* Liao, C., Zhou, T., Xu, D., Cooper, M. G., Engwirda, D., Li, H.-Y., & Leung, L. R. (2023). Topological relationship-based flow direction modeling: Mesh-independent river networks representation. Journal of Advances in Modeling Earth Systems, 15, e2022MS003089. https://doi.org/10.1029/2022MS003089
