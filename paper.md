---
title: 'pyflowline: a mesh-independent river network generator for hydrologic models'

tags:
  - Python
  - hydrologic model
  - river networks
  - mesh
  - geographic information system

authors:
  - name: Chang Liao
    orcid: 0000-0002-7348-8858    
    affiliation: 1

  - name: Matt G Cooper
    orcid: 0000-0002-0165-209X 
    affiliation: 1

affiliations:
 - name: Atmospheric Sciences and Global Change, Pacific Northwest National Laboratory, Richland, WA, USA
   index: 1 
date: 23 Mar 2023

bibliography: paper.bib
---

# Summary

River networks are crucial in hydrologic and Earth system models. Accurately representing river networks in spatially-distributed hydrologic models requires considering the model's spatial discretization and computational mesh. However, current methods of generating river networks for hydrologic models do not typically support unstructured meshes. Unstructured meshes offer numerous advantages over traditional, structured meshes. To overcome this limitation, we developed PyFlowline, a Python package that generates mesh-independent river networks. With PyFlowline, hydrologic modelers can generate conceptual river networks and their topological relationships for both structured and unstructured meshes.

# Statement of need

Generating a mesh cell-based conceptual river network for a given vector-based river network and arbitrary mesh geometry is a major challenge in hydrologic model configuration. Existing methods are typically limited to structured rectangular meshes such as 30m x 30m cartesian grids for high-resolution watershed-scale modeling or 0.5 degree x 0.5 degree geographic grids for global climate modeling. Structured meshes have fixed cartesian or geographic grid-cell sizes, which can be limiting for coupled land-river-ocean continuum simulations, which are often based on unstructured meshes [@Engwirda:2021]. However, unstructured meshes offer a flexible structure with variable grid-cell sizes and shapes. This flexibility makes them ideal for adapting to complex geometry such as river channels and coastlines.Thus, unstructured meshes are increasingly being adopted in hydrologic modeling.

Although unstructured meshes offer flexibility, the additional effort required to generate a computational river network that aligns with the mesh topology limits their adoption. A mesh-independent river network representation method that supports both unstructured and structured mesh-based hydrologic models could address this limitation. PyFlowline is a Python package that provides a framework for generating river networks for hydrologic models, meeting the identified need. Using an object-oriented programming approach, PyFlowline represents river network elements and mesh cell relationships. It relies on open-source Python libraries like GDAL and Cython for data input/output and spatial data operations.

The computational geometry algorithms used in PyFlowline are designed and implemented using a unified spherical framework, making it suitable for regional and global-scale simulations. Moreover, PyFlowline is mesh-independent, supporting both structured and unstructured meshes. It can quickly adopt other mesh types, such as triangulated irregular networks (TIN) or discrete global grid systems (DGGs). PyFlowline is a core component of the HexWatershed model, a mesh-independent flow direction model. Several scientific studies focused on coupled Earth system models [@Feng:2022; @Liao:2022] have utilized PyFlowline. A workshop tutorial has also been provided online and in person to support its implementation


# Model features

PyFlowline uses Python's OOP architecture to describe river networks using three essential elements: segments, reaches, and confluences. When applicable, river networks are processed as objects throughout the package.

![The data model. \label{fig:oop}](https://github.com/changliao1025/pyflowline/blob/main/docs/figures/basic_element.png?raw=true)

Pyflowline provides the following features:

PyFlowline provides several key features, including

1. Support for both structured and unstructured meshes, with JSON as the default file I/O format. For geospatial datasets such as vector river networks, GEOJSON is used.
2. Regional and global-scale processing capabilities through the use of fast Cython- and (global-only) AABB tree-based algorithms.
3. Built-in visualization functions based on the Python Matplotlib package, making it easy to visualize and analyze the PyFlowline model outputs.

# State of the field

PyFlowline is the only modeling software that can generate river networks on unstructured meshes. Model documentation is hosted at https://pyflowline.readthedocs.io/en/latest/, including a case study for the Susquehanna River Basin in the Mid-Atlantic region of the United States.


# Acknowledgment

The model described in this repository was supported by the following:

* the Earth System Model Development and Regional and Global Model Analysis program areas of the U.S. Department of Energy, Office of Science, Biological and Environmental Research program as part of the multi-program, collaborative Integrated Coastal Modeling (ICoM) project.

* the Earth System Model Development and Regional and Global Model Analysis program areas of the U.S. Department of Energy, Office of Science, Biological and Environmental Research program as part of the multi-program, collaborative Interdisciplinary Research for Arctic Coastal Environments (InteRFACE) project.

A portion of this research was performed using PNNL Research Computing at Pacific Northwest National Laboratory. 

PNNL is operated for DOE by Battelle Memorial Institute under contract DE-AC05-76RL01830.


# References

