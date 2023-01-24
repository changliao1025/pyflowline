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
  
  - name: Matt Cooper
    orcid: 0000-0002-0165-209X
    affiliation: 1

affiliations:
 - name: Atmospheric Sciences and Global Change, Pacific Northwest National Laboratory, Richland, WA, USA
   index: 1 
date: 11 Jan 2023

bibliography: paper.bib
---

# Summary

**River networks** are important features in hydrologic models. 

For spatially-distributed hydrologic models, the representation of river networks also needs to consider the model's spatial discretization, i.e., meshes.

Existing methods generally do not support unstructured meshes.

Therefore we developed a mesh-independent river network generation Python package, i.e., `PyFlowline` to close the gap. 

This package can be used to generate the **conceptual river networks** and their topologic relationship.

It also supports both **structured** and **unstructured** meshes.

# Statement of need

For a given **vector river network** and any **mesh**, generating the mesh cell-based conceptual river network remains challenging. 

Existing methods can only accept structured rectangle meshes, and cannot be used if the hdyrologic models use the unstructured meshes.

As a result, there is a need to develop a mesh-independent river network representation method for unstructured mesh-based hydrologic models.

`PyFlowline` is a Python package to generate river networks for hydrologic models. It uses the object-oriented programming (OOP) approach to represent the river networks and mesh cell relationships. It relies on several existing open source Python libraries, including the Geospatial Data Abstraction Library (`GDAL`) and `Cython` for data I/O and computational geometry.

It was designed using a unified framework that all the computational geometry algorithms are implemented on a sphere, which allows both regional and global scale simulations. It was designed using a mesh-independent approach so that both structured and unstructured meshes are supported. Other mesh structures such as triangular irregular networks (`TIN`) or discrete global grid systems (`DGGs`) can be easily adopted.

`PyFlowline` is a core component within the `HexWatershed` model, which is a mesh-independent flow direction model. `PyFlowline` has supported several scientific studies forcing on coupled Earth system models [@Feng:2022; @Liao:2022]. A workshop tutorial was also provided online and in-person. 

# Model features

Pyflowline uses Python's OOP architecture to describe the river network using three basic elements (i.e., segment, reach, confluence.) and processes them as objects throughout the package when applicable. 

![The data model. \label{fig:oop}](https://github.com/changliao1025/pyflowline/blob/main/docs/figures/basic_element.png?raw=true)

Pyflowline provides the following features:

1. It uses JSON as the default file I/O format. For spatial datasets, i.e., vector river network, GEOJSON is used.
2. It supports both structured and unstructured meshes.
3. It supports both regional scale and global scale (through AABB tree and Cython) simulations.
4. It provides built-in visualization functions based on the Python matplotlib package.

# Example

A case study was performed for the Susquehanna River Basin (SRB).

Screenshots of before and after river networks at zoom-in regions are used to illustrate the capability of the model.

# Acknowledgment

The model described in this repository was supported by the following:

* Earth System Model Development and Regional and Global Modeling and Analysis program areas of the U.S. Department of Energy, Office of Science, Office of Biological and Environmental Research as part of the multi-program, collaborative Integrated Coastal Modeling (ICoM) project.

* Earth System Model Development and Regional and Global Modeling and Analysis program areas of the U.S. Department of Energy, Office of Science, Office of Biological and Environmental Research as part of the multi-program, collaborative Interdisciplinary Research for Arctic Coastal Environments (InteRFACE) project.

A portion of this research was performed using PNNL Research Computing at Pacific Northwest National Laboratory. 

PNNL is operated for DOE by Battelle Memorial Institute under contract DE-AC05-76RL01830.

# References

* Liao. C. Cooper, M (2022) Pyflowline: a mesh independent river network generator for hydrologic models. Zenodo.
https://doi.org/10.5281/zenodo.6407299
