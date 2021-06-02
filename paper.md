---
title: 'PyStream: A Python package for stream networks processing'

tags:
  - Python
  - Hydrology
  - geographic information system

authors:
  - name: Chang Liao
    orcid: 0000-0002-7348-8858    
    affiliation: 1

affiliations:
 - name: Atmospheric Sciences and Global Change, Pacific Northwest National Laboratory, Richland, WA, USA
   index: 1 
date: 20 May 2021

bibliography: paper.bib
---

# Summary

Computational hydrologic simulation often requires high quality river network information as inputs. Preparation of river network information can be challenging, especially if the hydrology model is based on unstructured mesh. In this study, we develope a Python package to automate river network preparation. It resolves issues including disconnected river, braided network. It also produces spatially-distributed cell-to-cell topology information, which can be used for precise river routing. The package supports various mesh types including traditional projected coordinate system, geographic coordinate system and unstructured mesh.

# Statement of need






# Algorithms and implementation

PyStream takes advantage of Python language's object oriented programming (OOP) architect. The river network and all its elements (i.e. segment, reach, confluence.) are described as objects. These objects are processed throughout the package when applicable. A brief overview of the features provided by PyStream is list in Table 1.
1. It uses JSON/GeoJSON as the internal exchange file format. Most objects within the package can be imported and exported during runtime;
2. It relies on the open source Geospatial Data Abstraction Library (GDAL) and a few other Python packages to process geospatial data type;
3. It supports traditional projected coordinate system (PCS), geographic coordinate system (GCS), hexagon and other unstructured mesh for spatially distributed simulations;
4. It provides several algorithms to automate the river network simplification, generalization, etc.

* Overall model structure
![The workflow of PyStream. \label{fig:workflow}](https://github.com/changliao1025/pystream/blob/master/figure/workflow.png?raw=true)

# Example results


* Braided loop removal


# Acknowledgement

The model described in this repository was supported by:

* Laboratory Directed Research and Development (LDRD) Program Quickstarter project at Pacific Northwest National Laboratory. 
* Earth System Model Development and Regional and Global Modeling and Analysis program areas of the U.S. Department of Energy, Office of Science, Office of Biological and Environmental Research as part of the multi-program, collaborative Integrated Coastal Modeling (ICoM) project.
* U.S. Department of Energy Office of Science Biological and Environmental Research through the Earth and Environmental System Modeling program as part of the Energy Exascale Earth System Model (E3SM) project. 

A portion of this research was performed using PNNL Research Computing at Pacific Northwest National Laboratory. 

PNNL is operated for DOE by Battelle Memorial Institute under contract DE-AC05-76RL01830.

# References

1. Liao, C., Tesfa, T., Duan, Z., & Leung, L. R. (2020). Watershed delineation on a hexagonal mesh grid. Environmental Modelling & Software, 128, 104702.
