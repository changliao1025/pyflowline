###########
Application
###########

********
Overview
********

An example is provided within the `examples` folder. This example contains a case study in the `Susquehanna` watershed with several python scripts (corresponding to four mesh types). For example, the MPAS mesh type-based case is explained here.

****************
Model simulation
****************


================
Step 1
================

The example 'run_simulation_mpas.py' script import a few packages and functions.

::

    import os, sys
    from pathlib import Path
    from os.path import realpath
    from pyflowline.pyflowline_read_model_configuration_file import pyflowline_read_model_configuration_file

The 'pyflowline_read_model_configuration_file' function reads in a JSON configuration file and loads all the necessary model parameters. 


================
Step 2
================

The script setups some paths, which should be adjusted based on a real case.

::   
    
    sPath_parent = str(Path(__file__).parents[2]) # data is located two dir's up
    sPath_data = realpath( sPath_parent +  '/data/susquehanna' )
    sWorkspace_input =  str(Path(sPath_data)  /  'input')
    sWorkspace_output=  str(Path(sPath_data)  /  'output')
    sWorkspace_output=  '/compyfs/liao313/04model/pyflowline/susquehanna'

================
Step 3
================

Check the configuration file:

::   

    sFilename_configuration_in = realpath( sPath_parent +  '/examples/susquehanna/pyflowline_susquehanna_mpas.json' )
    if os.path.isfile(sFilename_configuration_in):
        pass
    else:
        print('This configuration does not exist: ', sFilename_configuration_in )

================
Step 4
================

Set up case information and read the configuration file.

::   

    iCase_index = 17
    Mesh = 'mpas'
    Date='20220901'


    oPyflowline = pyflowline_read_model_configuration_file(sFilename_configuration_in, \
       iCase_index_in=iCase_index, sDate_in=sDate)
    oPyflowline.aBasin[0].dLatitude_outlet_degree=39.462000
    oPyflowline.aBasin[0].dLongitude_outlet_degree=-76.009300

================
Step 5
================

Setup the model and run the three major steps.

::   

    oPyflowline.setup()
    oPyflowline.flowline_simplification()
    Cell = oPyflowline.mesh_generation()
    oPyflowline.reconstruct_topological_relationship(aCell)

================
Step 6
================
Analyze and export the model outputs.

::   

    oPyflowline.analyze()
    oPyflowline.evaluate()
    oPyflowline.export()

================
Step 7
================

Optionally, the user can also visualize the model outputs using the following method.

::

    aExtent_full = [-78.5,-75.5, 39.2,42.5]
    sFilename =  'filtered_flowline.png'
    oPyflowline._plot(sFilename, sVariable_in = 'flowline_filter', aExtent_in =aExtent_full  )
    
    sFilename =  'conceptual_flowline_with_mesh.png'
    oPyflowline._plot(sFilename,  iFlag_title=1 ,sVariable_in='overlap',   aExtent_in =aExtent_full )  