---
title: 'pyflowline: a mesh-independent river network generator for hydrologic models'

tags:
  - Python
  - Hydrology
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

Spatially-distributed hydrologic models often require high-quality river network information as inputs. Generation of river network information can be challenging, especially when the hydrologic models use unstructured meshes. We developed a mesh-independent river network generation Python package to address this challenge. It resolves issues, including disconnected and braided river networks. It also produces spatially-distributed cell-to-cell topology information, which can be used for river routing. The package supports various mesh types, including traditional rectangle and unstructured meshes.

# Statement of need
For a given vector river network and spatial discretization (mesh), generating the mesh cell-based conceptual river network remains challenging. Most existing methods can only accept structured rectangle meshes.
As a result, there is a need to develop a mesh-independent river network representation method for unstructured mesh-based hydrologic models.

# Model features
Pyflowline uses Python's object-oriented programming (OOP) architecture to describe the river network and its elements (i.e., segment, reach, confluence.) as objects processed throughout the package when applicable. 

![The data model. \label{fig:oop}](https://github.com/changliao1025/pyflowline/blob/main/docs/figures/basic_element.png?raw=true)


Pyflowline provides the following features:

1. It uses JSON as the file I/O format. For spatial datasets, i.e., vector river network, GEOJSON is used.
2. It supports both structured and unstructured meshes.
3. It supports both regional scale and global scale (through AABB tree and Cython) simulations.


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
