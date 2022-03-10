# Overview 
This repository is a tutorial to demonstrate the capability and workflow of the PyFlowline model. In this demo, we will use the Model for Prediction Across Scales (MPAS) mesh (see below) as an example.

# Pyflowline
The Pyflowline model is a Python package to generate conceptual river networks for hydrologic models. PyFlowline is mesh independent, meaning you can apply it to almost any mesh system including the tradition rectangle mesh, Triangulated Irregular Network (TIN) mesh and MPAS mesh.


# Installation
The full deployment of PyFlowline is still under development. It can be installed through either Pythin PyPI or the Conda system, which is recommended because of the dependency packages.

As of right now, you can install PyFlowline using the following steps:

1. install the dependency packages through Conda 

2. install PyFlowline through the PyPI:
    pip install pyflowline

3. (Optional) Install the Python JupterNote to run this tutorial.


# Usage
We use the notebook.py example file under the the notebook directory to showcase the model workflow.
An additional Python package is required for the visualization purpose. 

The follow steps are recommended:

1. Open the terminal or use your preferred Conda application to create a new Conda environment:

    * conda create --name pyflowline python=3.8

2. Activate the newly crated conda environment

    * conda activate pyflowline

3. Install dependency packages using conda

    * conda install -c conda-forge numpy

    * conda install -c conda-forge shapely

    * conda install -c conda-forge netCDF4

    * conda install -c conda-forge gdal

4. Install PyFlowline

    pip install pyflowline

5. Install and setup the Python Jupyter Notebook

6. Clone this repository and set this environment as the workspace environment

7. Navigate to the notebook and run it in your preferred Python IDE.

Because of the Python package dependency issue, the visulization should use a different environment or using the QGIS.

# Acknowledgement

This work was supported by the Earth System Model Development program areas of the U.S. Department of Energy, Office of Science, Office of Biological and Environmental Research as part of the multi-program, collaborative Integrated Coastal Modeling (ICoM) project. The data used for model simulations can be downloaded through the USGS website (\url{https://www.usgs.gov/national-hydrography}). The Pyflowline model can be accessed through the Python Package Index service (\url{https://pypi.org/project/pyflowline/}). 

# Citation

* Liao, Chang, Tian Zhou, Donghui Xu, Richard Barnes, Gautam Bisht, Hong-Yi Li, Zeli Tan, et al. (02/2022AD) 2022. “Advances In Hexagon Mesh-Based Flow Direction Modeling”. Advances In Water Resources 160. Elsevier BV: 104099. doi:10.1016/j.advwatres.2021.104099.

* Liao, C., Tesfa, T., Duan, Z., & Leung, L. R. (2020). Watershed delineation on a hexagonal mesh grid. Environmental Modelling & Software, 128, 104702. https://doi.org/10.1016/j.envsoft.2020.104702



    

