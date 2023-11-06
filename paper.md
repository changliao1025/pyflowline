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

River networks are crucial in hydrologic and Earth system models. Accurately representing river networks in hydrologic models requires considering the model's spatial resolution and computational mesh. However, current river network representation methods often have several limitations (1)  vector-based; (2) they perform poorly at coarse resolution (3) they do not support unstructured meshes. To overcome these limitations, we developed PyFlowline, a Python package that generates mesh-independent river networks. With PyFlowline, hydrologic modelers can generate conceptual river networks at various spatial resolutions for both structured and unstructured computational meshes. The generated river network datasets can be used by hydrologic models across scales.

# Statement of need

For hydrologic modelers, river networks are a key input for hydrologic models. 
While some hydrologic models accept vector-based river networks [@mizukami_2016_GMD; @Schwenk:2021], others only accept mesh cell-based representations, or require a generation method from a vector-based river network. 
Currently, generating a mesh cell-based river network from a given vector-based river network and arbitrary computational mesh is a major challenge.
Existing methods are typically limited to structured rectangular meshes, such as 30m x 30m cartesian grids for high-resolution watershed-scale modeling or 0.5 degree x 0.5 degree geographic grids for global climate modeling. In PyFlowline, we define structured meshes (e.g., lat-lon, raster files with projections, and hexagon) as those with fixed cell sizes and shapes and unstructured meshes as those with variable cell sizes and shapes.

Structured mesh-based methods use fixed cartesian or geographic cell sizes, which have several limitations: (1) they cannot resolve fine-scale river networks at coarse cell resolutions (>1km), and (2) they cannot be seamlessly coupled with other unstructured mesh-based hydrologic models such as oceanic models [@Engwirda:2021]. In contrast, unstructured meshes offer a flexible structure with variable grid-cell sizes and shapes. This flexibility makes them ideal for adapting to complex geometry such as river channels and coastlines. Besides, unstructured meshes provide the flexibility to couple different hydrologic models under a unified framework.
Thus, unstructured meshes are increasingly being adopted in hydrologic modeling.

Although unstructured meshes offer these flexibilities, additional efforts are required to generate conceptual river networks that capture real-world river networks across different spatial scales.

A mesh-independent river network representation method that preserves topological relationships across scales could address this limitation. 
PyFlowline is a Python package that provides a framework for generating river networks for hydrologic models, meeting the identified need. Using an object-oriented programming approach, PyFlowline represents river network elements and mesh cell relationships. It relies on open-source Python libraries like GDAL and Cython for data input/output and spatial data operations.

The computational geometry algorithms used in PyFlowline are designed and implemented using a unified spherical framework, making it suitable for regional and global-scale simulations. PyFlowline uses topological relationships to capture the river networks so they are preserved even at coarse spatial resolutions.
Moreover, PyFlowline is mesh-independent, supporting both structured and unstructured meshes. It can quickly adopt other mesh types, such as triangulated irregular networks (TIN) or discrete global grid systems (DGGs) [@Sahr:2011]. PyFlowline is a core component of the HexWatershed model, a mesh-independent flow direction model. Several scientific studies focused on coupled Earth system models [@Feng:2022; @Liao:2023; @Liao:2023b] have utilized PyFlowline. A workshop tutorial has also been provided online and in person to support its implementation.


# Model features

PyFlowline uses Python's OOP architecture to describe river networks using three essential elements: segments, reaches, and confluences. When applicable, river networks are processed as objects throughout the package.

![The data model. A vertex class object represents a point on the Earth surface. It have three coordiantes. An edge class object represents a directed line between two points. Besides, it has a length attribute. A flowline class object represents a list of connected lines. \label{fig:oop}](https://github.com/changliao1025/pyflowline/blob/main/docs/figures/basic_element.png?raw=true)

PyFlowline provides several key features, including

1. Support for both structured and unstructured meshes, with JSON as the default file I/O format. For geospatial datasets such as vector river networks, GEOJSON is used.
2. Regional and global-scale processing capabilities through the use of fast Cython- and (global-only) AABB tree-based algorithms.
3. Built-in visualization functions (experimental) based on the Python Matplotlib package [@LiaoPyearth:2022], making it easy to visualize and analyze the PyFlowline model outputs.

# State of the field

Existing river network representation methods often fall into three categories, each with associated strengths and weaknesses:

1. Vector-based. Hydrologic models that use this method can represent fine-scale river networks, but cannot couple river and land without one-to-one mapping between river segments and land model cells [@Schwenk:2021];
2. High-resolution DEM-based. Flow networks derived from structured rectangle-grid DEMs are widely available, but resolving fine-scale river networks requires grids with very high spatial resolution (e.g., 30m x 30m or finer) [@Esri:2011];
3. Upscaling-based methods address the scale mismatch between coarse-resolution grids and fine-scale river networks, but only support structured geographic grids (e.g., 0.5 degree x 0.5 degree) at coarse resolutions [@Wu:2012]. This method often cannot provide global coverage, including Greenland and the Antarctic.

PyFlowline is the only modeling software that provides these unique features:

1. It can generate river networks on unstructured meshes; 
2. It uses topological relationships to capture river networks precisely; 
3. It can be applied at both high and coarse resolutions; 
4. It can provide global coverage, including Greenland and the Antarctic.

Model documentation is hosted at https://pyflowline.readthedocs.io/en/latest/, including a case study for the Susquehanna River Basin in the Mid-Atlantic region of the United States.

# Acknowledgment

The model described in this repository was supported by the following:

* the Earth System Model Development and Regional and Global Model Analysis program areas of the U.S. Department of Energy, Office of Science, Biological and Environmental Research program as part of the multi-program, collaborative Integrated Coastal Modeling (ICoM) project.

* the Earth System Model Development and Regional and Global Model Analysis program areas of the U.S. Department of Energy, Office of Science, Biological and Environmental Research program as part of the multi-program, collaborative Interdisciplinary Research for Arctic Coastal Environments (InteRFACE) project.

* the Next Generation Ecosystem Experiments-Tropics project, funded by the U.S. Department of Energy, Office of Science, Office of Biological and Environmental Research at Pacific Northwest National Laboratory. 

A portion of this research was performed using PNNL Research Computing at Pacific Northwest National Laboratory. 

PNNL is operated for DOE by Battelle Memorial Institute under contract DE-AC05-76RL01830.

# References

