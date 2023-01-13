##########
Data model
##########

*********
Basic
*********

River networks are represented using three basic elements: vertex, edge, and flowline.

.. image:: ../figures/basic_element.png
  :width: 400
  :alt: Basic element


****************************************************
Spatial references and computational geometry
****************************************************

All the internal data elements use the geographic coordinate system (GCS).

All the computational geometry algorithms are based on GCS:

+------------------------+------------+----------+----------+
|                        | Input      | Output   | Algorithm|
|                        |            |          |          |
+========================+============+==========+==========+
| Location               | vertex(lon, lat)    |  vertex(lon, lat)         |  |
+------------------------+------------+----------+----------+
| Distance               | vertex A, B        | Distance (m)      |  Great circle    |
+------------------------+------------+----------+----------+
| Area                   | vertex A, B, C, ... D        | Distance (m2)   |   Spheric area       |
+------------------------+------------+----------+----------+