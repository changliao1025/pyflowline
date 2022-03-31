# Pyflowline
The Pyflowline model is a Python package to generate conceptual river networks for hydrologic models. 

PyFlowline is mesh independent, meaning you can apply it to both structured (e.g., tradition rectangle mesh, latitude-longitude, hexagon) and unstructured mesh systems (e.g., Triangulated Irregular Network (TIN) mesh and MPAS mesh).

This package runs three steps in general:
1. pre-process the existing vector river flowlines. In this step, PyFlowline will check the vector dataset and correct undesired flaws such as braided rivers. The output is a much simplified flowline dataset.
2. mesh generation. In this step, PyFlowline will generate various structured meshes (e.g., rectangle, hexagon). Users can also use other unstructured meshes, which can be coverted to the PyFlowline format.
3. Re-construct the topological relationship using the mesh and flowline intersections. In this step, PyFlowline will build the topological relationship between mesh cells using the vector flowline and mesh intersection.



# Installation


As of right now, the easiest way you can install PyFlowline is through the Conda platform:

1. install the package from the conda forge channel
    conda install -c conda-forge pyflowline

# Usage
We use the notebook.py example file under the the notebook directory to showcase the model workflow.
An additional Python package is required for the visualization purpose. 

The follow steps are recommended:

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

# Acknowledgement

This work was supported by the Earth System Model Development program areas of the U.S. Department of Energy, Office of Science, Office of Biological and Environmental Research as part of the multi-program, collaborative Integrated Coastal Modeling (ICoM) project. The data used for model simulations can be downloaded through the USGS website (https://www.usgs.gov/national-hydrography). The Pyflowline model can be accessed through the Python Package Index service (https://pypi.org/project/pyflowline/). 

# Citation

* Liao, Chang, Tian Zhou, Donghui Xu, Richard Barnes, Gautam Bisht, Hong-Yi Li, Zeli Tan, et al. (02/2022AD) 2022. “Advances In Hexagon Mesh-Based Flow Direction Modeling”. Advances In Water Resources 160. Elsevier BV: 104099. doi:10.1016/j.advwatres.2021.104099.

* Liao, C., Tesfa, T., Duan, Z., & Leung, L. R. (2020). Watershed delineation on a hexagonal mesh grid. Environmental Modelling & Software, 128, 104702. https://doi.org/10.1016/j.envsoft.2020.104702



    

