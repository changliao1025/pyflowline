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

3. (Optional) Install the visualization package through Conda:
    conda install -c conda-forge

4. (Optional) Install the Python JupterNote to run this tutorial.


# Usage
We use the notebook.py example file under the the notebook directory to showcase the model workflow.
An additional Python package is required for the visualization purpose. 

The follow steps are recommended:
1. Open the terminal or use your preferred Conda application to create a new Conda environment:
    conda create --name pyflowline python=3.8
2. Activate the newly crated conda environment
    conda activate pyflowline
3. Install dependency packages using conda
    conda install -c conda-forge numpy
    conda install -c conda-forge shapely
    conda install -c conda-forge netCDF4
    conda install -c conda-forge gdal
4. Install PyFlowline
    pip install pyflowline
5. Install and setup the Python Jupyter Notebook
6. Install the visualization Package
    conda install -c conda-forge matplotlib   
    conda install -c conda-forge cartopy 
    conda install -c conda-forge geopandas 
7. Clone this repository and set this environment as the workspace environment
8. Navigate to the notebook and run it in your preferred Python IDE.




    

