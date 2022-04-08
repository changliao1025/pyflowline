### `Pyflowline`

[![DOI](https://zenodo.org/badge/368338554.svg)](https://zenodo.org/badge/latestdoi/368338554)

Pyflowline: a mesh independent river network generator for hydrologic models. 

PyFlowline is mesh independent, meaning you can apply it to both structured (e.g., tradition rectangle mesh, latitude-longitude, hexagon) and unstructured mesh systems (e.g., Triangulated Irregular Network (TIN) mesh and MPAS mesh).

This package runs three steps in general:
1. pre-process the existing vector river flowlines. In this step, PyFlowline will check the vector dataset and correct undesired flaws such as braided rivers. The output is a much simplified flowline dataset.
2. mesh generation. In this step, PyFlowline will generate various structured meshes (e.g., rectangle, hexagon). Users can also use other unstructured meshes, which can be converted to the PyFlowline format.
3. Re-construct the topological relationship using the mesh and flowline intersections. In this step, PyFlowline will build the topological relationship between mesh cells using the vector flowline and mesh intersection.



### `Installation`

As of right now, the easiest way you can install PyFlowline is through the Conda platform:

1. install the package from the conda forge channel

    conda install -c conda-forge pyflowline

### `Usage`

We use the notebook.py example file under the the notebook directory to showcase the model workflow.
An additional Python package is required for the visualization purpose. 

The following steps are recommended:

1. Open the terminal or use your preferred Conda application to create a new Conda environment:

    * conda create --name pyflowline python=3.8

2. Activate the newly crated conda environment

    * conda activate pyflowline

3. Install dependency packages using conda

    * conda install -c conda-forge pyflowline

4. Install and setup the Python Jupyter Notebook

5. Clone this repository and set this environment as the workspace environment

6. Navigate to the notebook and run it in your preferred Python IDE.

Because of the Python package dependency issue, the visulization should use a different environment or using the QGIS.

### `Acknowledgement`

This work was supported by the Earth System Model Development program areas of the U.S. Department of Energy, Office of Science, Office of Biological and Environmental Research as part of the multi-program, collaborative Integrated Coastal Modeling (ICoM) project. The data used for model simulations can be downloaded through the USGS website (https://www.usgs.gov/national-hydrography). The Pyflowline model can be accessed through the Python Conda system (https://anaconda.org/conda-forge/pyflowline). 

### `License`

Copyright © 2022, Battelle Memorial Institute

1. Battelle Memorial Institute (hereinafter Battelle) hereby grants permission to any person or entity lawfully obtaining a copy of this software and associated documentation files (hereinafter “the Software”) to redistribute and use the Software in source and binary forms, with or without modification. Such person or entity may use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and may permit others to do so, subject to the following conditions:

* Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimers.

* Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

* Other than as used herein, neither the name Battelle Memorial Institute or Battelle may be used in any form whatsoever without the express written consent of Battelle.

2. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL BATTELLE OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

### `References`

There are a number of publications that describe the algorithms used in `HexWatershed` in detail. If you make use of `HexWatershed` in your work, please consider including a reference to the following:

* Liao, Chang, Tian Zhou, Donghui Xu, Richard Barnes, Gautam Bisht, Hong-Yi Li, Zeli Tan, et al. (02/2022AD) 2022. “Advances In Hexagon Mesh-Based Flow Direction Modeling”. Advances In Water Resources 160. Elsevier BV: 104099. 
https://doi.org/10.1016/j.advwatres.2021.104099.

* Liao, C., Tesfa, T., Duan, Z., & Leung, L. R. (2020). Watershed delineation on a hexagonal mesh grid. Environmental Modelling & Software, 128, 104702. https://doi.org/10.1016/j.envsoft.2020.104702

* Liao. C. (2022) Pyflowline: a mesh independent river network generator for hydrologic models. Zenodo.
https://doi.org/10.5281/zenodo.6407299





    

