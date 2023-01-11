### Pyflowline

[![DOI](https://zenodo.org/badge/368338554.svg)](https://zenodo.org/badge/latestdoi/368338554)

Pyflowline: a mesh-independent river network generator for hydrologic models. 

PyFlowline is mesh independent, meaning you can apply it to both structured (e.g., traditional rectangle mesh, latitude-longitude, hexagon) and unstructured mesh systems (e.g., Triangulated Irregular Network (TIN) mesh and Model for Prediction Across Scales (MPAS) mesh).

This package generates the mesh cell-based conceptual river network using the following steps:
1. `Flowline simplification`: PyFlowline checks the vector dataset and corrects undesired flowlines, such as braided rivers. 
2. `Mesh generation`: PyFlowline generates structured meshes (e.g., rectangle, hexagon) or imports user-provided unstructured meshes into the PyFlowline-compatible GEOJSON format.
3. `Topological relationship reconstruction`: PyFlowline reconstructs the topological relationship using the mesh and flowline intersections. 

### Installation

In most cases, you can install PyFlowline through the Conda system:

    conda install -c conda-forge pyflowline

In rare cases, if the `Cython` and `AABB tree` features are needed for high-performance simulations, please refer to `Documentation` (https://pyflowline.readthedocs.io/) for details on how to enable them.

### Usage

We use the python scripts file under the `examples` directory to demonstrate the model capability.

The following steps are recommended:

1. Open the terminal or use your preferred Conda application to create a new Conda environment:

    * conda create --name pyflowline

2. Activate the newly created conda environment

    * conda activate pyflowline

3. Install dependency packages using conda

    * conda install -c conda-forge pyflowline

4. Clone this repository and set this environment as the workspace environment

5. Navigate to the examples and run the python scripts in your preferred Python IDE.
The visualization features require additional Python packages; please refer to `Documentation` (https://pyflowline.readthedocs.io/) for details on visualization.

### Acknowledgment

This work was supported by the Earth System Model Development program areas of the U.S. Department of Energy, Office of Science, Office of Biological and Environmental Research as part of the multi-program, collaborative Integrated Coastal Modeling (ICoM) project and the Interdisciplinary Research for Arctic Coastal Environments (InteRFACE) project. 
The data used for model simulations can be downloaded through the USGS website (https://www.usgs.gov/national-hydrography). 
The Pyflowline model can be accessed through the Python Conda system (https://anaconda.org/conda-forge/pyflowline). 

### License

Copyright © 2022, Battelle Memorial Institute

1. Battelle Memorial Institute (hereinafter Battelle) hereby grants permission to any person or entity lawfully obtaining a copy of this software and associated documentation files (hereinafter “the Software”) to redistribute and use the Software in source and binary forms, with or without modification. Such person or entity may use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software and may permit others to do so, subject to the following conditions:

* Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimers.

* Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

* Other than as used herein, neither the name Battelle Memorial Institute or Battelle may be used in any form whatsoever without the express written consent of Battelle.

2. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL BATTELLE OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

### References

There are several publications that describe the algorithms used in `PyFlowline` in detail. If you make use of `PyFlowline` in your work, please consider including a reference to the following:

* Liao. C. Cooper, M (2022) Pyflowline: a mesh independent river network generator for hydrologic models. Zenodo.
https://doi.org/10.5281/zenodo.6407299





    

