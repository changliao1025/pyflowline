##########
Data model
##########

*********
Basic
*********

River networks are represented using three basic elements: vertex, edge, and flowline.

.. image:: ../../figures/basic_element.png
  :width: 400
  :alt: Basic elements

Within PyFlowline, these three elements are combined with several other data structures.

.. image:: ../../figures/structure_pyflowline.png
  :width: 400
  :alt: PyFlowline structure. 
  
This figure illustrates a domain containing two watersheds/basins. Each basin has an outlet. Within each basin, there are several subbasins and confluences. The lower right is a zoom-in view of a flowline.

****************************************************
Spatial references and computational geometry
****************************************************

All the internal data elements use the geographic coordinate system (GCS).

All the computational geometry algorithms are based on GCS:

+------------------------+-----------------------+-------------------+--------------+
|                        | Input                 | Output            | Algorithm    |
|                        |                       |                   |              |
+========================+=======================+===================+==============+
| Location               | vertex(lon, lat)      |  vertex(lon, lat) |              |
+------------------------+-----------------------+-------------------+--------------+
| Distance               | vertex A, B           | Distance (m)      | Great circle |
+------------------------+-----------------------+-------------------+--------------+
| Area                   | vertex A, B, C, ... D | Distance (m2)     | Spheric area |
+------------------------+-----------------------+-------------------+--------------+


*********
File I/O
*********

==============================
Inputs
==============================


PyFlowline uses two configuration files to manage all the input information. Within this configuration file, it stores major model input parameters and paths to input files. 

These two configuration files have a parent-child relationship:
1. The parent configuration file stores parameters for the whole domain, and
2. The child configuration file stores parameters for every single watershed.

Both configuration files are in JSON format.

An example parent JSON file is provided below:

::

    {
        "sFilename_model_configuration": "/qfs/people/liao313/workspace/python/pyflowline/pyflowline/config/hexwatershed_susquehanna_mpas.json",
        "sWorkspace_data": "/people/liao313/data",    
        "sWorkspace_output": "/compyfs/liao313/04model/pyflowline/susquehanna",
        "sWorkspace_project": "/hexwatershed/susquehanna",
        "sWorkspace_bin": "/people/liao313/bin",
        "sRegion": "susquehanna",
        "sModel": "pyflowline",
        "sJob": "hex",   
        "iFlag_standalone": 1,      
        "iFlag_create_mesh": 1,
        "iFlag_save_mesh" :1 ,
        "iFlag_simplification": 1,
        "iFlag_intersect": 1,
        "iFlag_flowline":1,
        "iFlag_use_mesh_dem":1,
        "iFlag_global": 0,
        "iFlag_multiple_outlet": 0,
        "iFlag_rotation": 0, 
        "iCase_index": 1,
        "iMesh_type": 1,    
        "dLongitude_left": -79,
        "dLongitude_right": -74.5,
        "dLatitude_bot": 39.20,
        "dLatitude_top": 42.8,
        "dResolution_degree": 5000,
        "dResolution_meter": 5000,    
        "sDate": "20220110",        
        "sMesh_type": "mpas",       
        "sFilename_spatial_reference": "/qfs/people/liao313/workspace/python/pyhexwatershed_icom/data/susquehanna/input/boundary_proj_buff.shp",
        "sFilename_dem": "/qfs/people/liao313/workspace/python/pyhexwatershed_icom/data/susquehanna/input/dem_buff_ext.tif",     
        "sFilename_mesh_netcdf": "/qfs/people/liao313/data/icom/mesh/delaware_lnd_60_30_5_2_v1/lnd_cull_mesh.nc",    
        "sFilename_basins": "/qfs/people/liao313/workspace/python/pyflowline_icom/examples/susquehanna/pyflowline_susquehanna_basins.json"
    }

+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| Parameter                      | Data type  | Usage                                   | Default value  | Note                                |
|                                |            |                                         |                |                                     |
+================================+============+=========================================+================+=====================================+
| sFilename_model_configuration  | string     | The filename of the configuration file  | None           | It will be automatically generated  |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| sWorkspace_data                | string     | The workspace of data                   | None           | Unused                              |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| sWorkspace_output              | string     | The output workspace                    | None           | The output folder                   |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| sWorkspace_project             | string     | The project workspace                   | None           | Unused                              |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| sWorkspace_bin                 | string     | The workspace for binary executable     | None           | Reserved for HexWatershed model     |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| sRegion                        | string     | Study region                            | None           | None                                |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| sModel                         | string     | Model name                              | pyflowline     | None                                |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| sJob                           | string     | HPC batch job name                      | pyflowline     | None                                |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| iFlag_standalone               | int        | Flag to run pyflowlone standalone       |  1             | 0 when called by hexwatershed       |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| iFlag_create_mesh              | int        | Flag to create mesh                     |  1             | None                                |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| iFlag_save_mesh                | int        | Flag to save mesh                       |  1             | None                                |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| iFlag_simplification           | int        | Flag to simplification                  |  1             | None                                |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| iFlag_intersect                | int        | Flag to intersect                       |  1             | None                                |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| iFlag_flowline                 | int        | Flag for flowline                       |  1             | None                                |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| iFlag_use_mesh_dem             | int        | Flag to use DEM data                    |  0             | Not used                            |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| iFlag_global                   | int        | Flag to run on global scale             |  0             | None                                |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| iFlag_multiple_outlet          | int        | Flag to run with multi-outlet           |  0             | None                                |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| iFlag_rotation                 | int        | Flag for hexagon rotation               |  0             | None                                |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| iCase_index                    | int        | Index of case                           |  1             | None                                |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| iMesh_type                     | int        | Type of mesh                            |  1             | None                                |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| dLongitude_left                | float      | Boundary                                |  -180          | None                                |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| dLongitude_right               | float      | Boundary                                |  +180          | None                                |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| dLatitude_bot                  | float      | Boundary                                |  -90           | None                                |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| dLatitude_top                  | float      | Boundary                                |  +90           | None                                |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| dResolution_degree             | float      | Resolution in degree                    |  1             | None                                |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| dResolution_meter              | float      | Resolution in meter                     |  5000          | None                                |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| sDate                          | string     | Date of simulation                      |  None          | None                                |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| sMesh_type                     | string     | Mesh type                               |  None          | None                                |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| sFilename_spatial_reference    | string     | Spatial reference                       |  None          | None                                |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| sFilename_dem                  | string     | DEM file                                |  None          | Reserved for HexWatershed model     |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| sFilename_mesh_netcdf          | string     | Netcdf mesh file                        |  None          |                                     |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| sFilename_basins               | string     | Filename of child JSON file             |  None          | None                                |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+


An example child JSON file is provided below:

::

    [
    {
        "dLatitude_outlet_degree": 39.4620,
        "dLongitude_outlet_degree": -76.0093,    
        "dAccumulation_threshold": 100000,
        "dThreshold_small_river": 10000,
        "iFlag_dam": 0,
        "iFlag_debug":1,
        "iFlag_disconnected": 0,
        "lBasinID": 1,
        "sFilename_dam": "/qfs/people/liao313/data/hexwatershed/susquehanna/auxiliary/ICoM_dams.csv",
        "sFilename_flowline_filter": "/qfs/people/liao313/workspace/python/pyhexwatershed_icom/data/susquehanna/input/flowline.geojson",
        "sFilename_flowline_raw": "/qfs/people/liao313/data/hexwatershed/susquehanna/vector/hydrology/allflowline.shp",
        "sFilename_flowline_topo": "/qfs/people/liao313/data/hexwatershed/susquehanna/auxiliary/flowline.csv"
    }
    ]

+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| Parameter                      | Data type  | Usage                                   | Default value  | Note                                |
|                                |            |                                         |                |                                     |
+================================+============+=========================================+================+=====================================+
| dLatitude_outlet_degree        | string     | The latitude of outlet                  | None           |                                     |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| dLongitude_outlet_degree       | string     | The longitude of outlet                 |                |                                     |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| dAccumulation_threshold        | string     | The flow accumulation threshold         |                |                                     |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| dThreshold_small_river         | string     | The small river threshold               |                |                                     |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| iFlag_dam                      | string     | Flag for dam burning                    |  0             |                                     |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| iFlag_debug                    | string     | Flag to turn on debug info              |  0             |                                     |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| iFlag_disconnected             | string     | Flag for disconnected flowline          |  0             |                                     |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| lBasinID                       | string     | Basin/watershed ID                      |  0             |                                     |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| sFilename_dam                  | int        | Filename of dam file                    |  1             | Only used for dam burning           |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| sFilename_flowline_filter      | int        | Filename of original flowline file      |                | GeoJSON format                      |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| sFilename_flowline_raw         | int        | Filename of flowline including dam      |                | Only used for dam burning           |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+
| sFilename_flowline_topo        | int        | Filename of dam topology                |                | Only used for dam burning           |
+--------------------------------+------------+-----------------------------------------+----------------+-------------------------------------+


==============================
Outputs
==============================

After the PyFlowline simulation, the output workspace has a structure like this:


::

    pyflowlinecase 
    ├── 00000001          
    │   ├── basin_info.json
    │   └── conceptual_flowline.geojson
    │   └── ...
    ├── 00000002          
    │   ├── basin_info.json
    │   └── conceptual_flowline.geojson
    │   └── ...
    ├── mpas_mesh_info.json          
    ├── mpas.geojson
    ├── run_pyflowline.py          
    ├── submit.job
    ├── stdout.out
    └── stderr.err

At the root directory, three files `submit.job`, `stdout.out`, `stderr.err` are HPC associated files.

The `run_pyflowline.py` is the python script that was ran by the HPC job. If you are running on a local machine, you can run this script directly.

The `mpas_mesh_info.json` is the model output that has all the information.

The `mpas.geojson` is the model generated mesh file in the GEOJSON format.

The sub-folders `00000001` et. al, are results for every watershed. Within each watershed sub-folder, there are both json and geojson result files.