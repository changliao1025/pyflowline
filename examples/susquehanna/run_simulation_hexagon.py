import os, sys
from pathlib import Path
from os.path import realpath
from datetime import datetime

#%% setup case information

sMesh = 'hexagon'
iCase_index = 1
dResolution_meter = 50000 # mesh resolution
sDate = datetime.now().strftime('%Y%m%d')

# Temporary patch to avoid bug in mesh generation function
iFlag_mesh_generation = False

#%% set up workspace path
sPath_parent = str(Path(__file__).parents[2]) # data is located two dir's up
import sys
sys.path.append(sPath_parent)

#%% set up workspace folders
sWorkspace_data = realpath( sPath_parent + '/data/susquehanna' )
sWorkspace_input = str(Path(sWorkspace_data) / 'input')
sWorkspace_output = str(Path(sWorkspace_data) / 'output')

#%% set the configuration file names.
sFilename_configuration_in = realpath( sWorkspace_input + '/' + 'pyflowline_susquehanna_hexagon.json' )
if os.path.isfile(sFilename_configuration_in):
    pass
else:
    print(f"The model configuration does not exist: \n{sFilename_configuration_in}\n")

# set the basin configuration file name
sFilename_basins_in = realpath( sWorkspace_input + '/' + 'pyflowline_susquehanna_basins.json' )
if os.path.isfile(sFilename_configuration_in):
    pass
else:
    print(f"The basin configuration does not exist: \n{sFilename_basins_in}\n")

#%% set the input data filenames
sFilename_flowline_filter = realpath( sWorkspace_input + '/' + 'flowline.geojson' )

sFilename_mesh_boundary = realpath( sWorkspace_input + '/' + 'boundary_wgs.geojson' )

#%% set the configuration file parameters

# To change a parameter in the json configuration file, pass the configuration file filename followed by a single key-value pair. For the basin configuration file only, pass the final argument iFlag_basin_in=1. Import the function used to change the parameters.
from pyflowline.configuration.change_json_key_value import change_json_key_value as pyflowline_change_json_key_value

# Set the path to the data folder
pyflowline_change_json_key_value(sFilename_configuration_in, 'sWorkspace_data', sWorkspace_data)

# Set the path to the input folder
pyflowline_change_json_key_value(sFilename_configuration_in, 'sWorkspace_input', sWorkspace_input)

# Set the path to the output folder
pyflowline_change_json_key_value(sFilename_configuration_in, 'sWorkspace_output', sWorkspace_output)

# Set the model configuration file name
pyflowline_change_json_key_value(sFilename_configuration_in, 'sFilename_model_configuration', sFilename_configuration_in)

# Set the basin configuration file name
pyflowline_change_json_key_value(sFilename_configuration_in, 'sFilename_basins', sFilename_basins_in)

# Set the domain boundary file name
pyflowline_change_json_key_value(sFilename_configuration_in, 'sFilename_mesh_boundary', sFilename_mesh_boundary)

# pyflowline_change_json_key_value(sFilename_configuration_in, 'sFilename_spatial_reference', sFilename_mesh_boundary)

# Note that a DEM can also be supplied. Set the path to the DEM with the sFilename_dem parameter. Although the DEM is not required by `pyflowline`, users will need a DEM to run `hexwatershed`. If a DEM is supplied, a boundary file is not required - the DEM boundary will be used instead. However, if a DEM is not supplied, then the model requires a domain boundary to create the mesh.

# If a boundary file is used, set the iFlag_mesh_boundary flag true
pyflowline_change_json_key_value(sFilename_configuration_in, 'iFlag_mesh_boundary', 1)

# Set the path to the flowline in the basin configuration file
pyflowline_change_json_key_value(sFilename_basins_in, 'sFilename_flowline_filter', sFilename_flowline_filter, iFlag_basin_in=1)

#%% read the configuration file to create a pyflowline case
from pyflowline.configuration.read_configuration_file import pyflowline_read_configuration_file
oPyflowline = pyflowline_read_configuration_file(sFilename_configuration_in,iCase_index_in=iCase_index, dResolution_meter_in=dResolution_meter, sDate_in=sDate)

#%% change model parameters as needed before running the model
oPyflowline.aBasin[0].dLatitude_outlet_degree=39.462000
oPyflowline.aBasin[0].dLongitude_outlet_degree=-76.009300

#%% run the model
oPyflowline.pyflowline_setup()
oPyflowline.pyflowline_flowline_simplification();

#%% generate the mesh.

# There is currently a bug in the mesh generation function when a boundary file is used instead of a DEM. Skip this step until that bug is fixed.

if iFlag_mesh_generation is True:

    aCell = oPyflowline.pyflowline_mesh_generation();
    oPyflowline.pyflowline_reconstruct_topological_relationship(aCell)

    # Export the model results for hexwatershed
    oPyflowline.pyflowline_export()
    oPyflowline.pyflowline_export_config_to_json()
    iCase_index= iCase_index+1

else:
    pass

# %%
print('Finished hexagon demo.')
print('Try opening the mesh and flowline files in QGIS for visualization.')

# %%
