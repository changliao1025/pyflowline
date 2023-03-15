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
date: 11 Jan 2023

bibliography: paper.bib
---

# Summary

River networks are important hydrologic features in both hydrologic and earth system models.
For spatially-distributed hydrologic models, the representation of river networks must also consider the model's spatial discretization, i.e., meshes.
Existing methods generally do not support unstructured meshes.
Therefore we developed a mesh-independent river network generation Python package, i.e., PyFlowline, to close this gap.
Hydrologic modelers can use this package to generate conceptual river networks and their topologic relationships.
It supports both structured and unstructured meshes.

# Statement of need

For a given vector-based river network and any mesh, generating the mesh cell-based conceptual river network remains challenging.

Existing methods can only accept structured rectangle meshes (e.g., 30 m x 30 m or 0.5 degrees x 0.5 degrees) and cannot be used if the hydrologic models use unstructured meshes.

As a result, there is a need to develop a mesh-independent river network representation method that supports unstructured mesh-based hydrologic models.

PyFlowline is a Python package to generate river networks for hydrologic models. It uses the object-oriented programming (OOP) approach to represent the river network elements and mesh cell relationships. It relies on several existing open-source Python libraries, including the Geospatial Data Abstraction Library (GDAL) and Cython, for data I/O and spatial data operations.

We use a unified spherical framework to design and implement all the computational geometry algorithms, allowing regional and global scale simulations. It is mesh-independent so that both structured and unstructured meshes are supported. Other mesh types, such as triangulated irregular networks (TIN) or discrete global grid systems (DGGs), can be quickly adopted.

PyFlowline is a core component within the HexWatershed model, which is a mesh-independent flow direction model. PyFlowline has supported several scientific studies forcing on coupled Earth system models [@Feng:2022; @Liao:2022; @Cooper:2022]. A workshop tutorial was also provided online and in person.

# Model features

Pyflowline uses Python's OOP architecture to describe the river network using three basic elements (i.e., segment, reach, confluence.) and processes them as objects throughout the package when applicable. 

![The data model. \label{fig:oop}](https://github.com/changliao1025/pyflowline/blob/main/docs/figures/basic_element.png?raw=true)

Pyflowline provides the following features:

1. It uses JSON as the default file I/O format. For geospatial datasets, i.e., vector river network, GEOJSON is used.
2. It supports both structured and unstructured meshes.
3. It supports regional and global scales (through Cython and AABB tree (global-only)) simulations.
4. It provides built-in visualization functions based on the Python Matplotlib package.

# Example

A case study was provided for the Susquehanna River Basin (SRB) in the Mid-Atlantic region of the United States.
You can either use the `Python` scripts within the `examples` folder or the `notebooks` to test the model.

Model documentation is hosted at https://pyflowline.readthedocs.io/en/latest/.

# Acknowledgment

The model described in this repository was supported by the following:

* Earth System Model Development and Regional and Global Modeling and Analysis program areas of the U.S. Department of Energy, Office of Science, Office of Biological and Environmental Research as part of the multi-program, collaborative Integrated Coastal Modeling (ICoM) project.

* Earth System Model Development and Regional and Global Modeling and Analysis program areas of the U.S. Department of Energy, Office of Science, Office of Biological and Environmental Research as part of the multi-program, collaborative Interdisciplinary Research for Arctic Coastal Environments (InteRFACE) project.

A portion of this research was performed using PNNL Research Computing at Pacific Northwest National Laboratory. 

PNNL is operated for DOE by Battelle Memorial Institute under contract DE-AC05-76RL01830.

# References

* Liao. C. Cooper, M (2022) Pyflowline: a mesh-independent river network generator for hydrologic models. Zenodo.
https://doi.org/10.5281/zenodo.6407299

* Liao, C., Zhou, T., Xu, D., Cooper, M. G., Engwirda, D., Li, H.-Y., & Leung, L. R. (2023). Topological relationship-based flow direction modeling: Mesh-independent river networks representation. Journal of Advances in Modeling Earth Systems, 15, e2022MS003089. https://doi.org/10.1029/2022MS003089
