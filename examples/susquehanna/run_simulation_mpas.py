import os, sys
from pathlib import Path
from datetime import datetime

#%% Define the case information (configuration file parameters)
sMesh = 'mpas'
sRegion = 'susquehanna'
iCase_index = 1
iFlag_simulation = 1
iFlag_visualization = 1
sDate = datetime.now().strftime('%Y%m%d')

#%% Define workspace paths and configuration file path parameters
oPath_parent = Path(__file__).parents[2] # data is located two dir's up
sys.path.append(str(oPath_parent))

# Define the full path to the input and output folders.
oFolder_input = oPath_parent.joinpath(
    'data', sRegion, 'input')
oFolder_output = oPath_parent.joinpath(
    'data', sRegion, 'output')

# Define the full path to the domain ("parent") configuration file.
oFilename_domain_config = oFolder_input.joinpath(
    'pyflowline_susquehanna_mpas.json')

# Define the full path to the individual basin ("child") configuration file.
oFilename_basins_config = oFolder_input.joinpath(
    'pyflowline_susquehanna_basins.json')

# Define the full path to the MPAS mesh file
oFilename_mesh_netcdf = oFolder_input.joinpath(
    'lnd_cull_mesh.nc')

# Define the full path to the input flowline file.
oFilename_flowline = oFolder_input.joinpath(
    'flowline.geojson')

# Define the full path to the boundary file used to clip the mesh.
oFilename_mesh_boundary = oFolder_input.joinpath(
    'boundary_wgs.geojson')

# Confirm that the domain configuration file exists.
if os.path.isfile(oFilename_domain_config):
    pass
else:
    print('The domain configuration file does not exist: \n', oFilename_domain_config)

#%% Update the domain (parent) configuration file

# The following parameters need to be set in the configuration file:
#   sWorkspace_output: full/path/to/output
#   sFilename_mesh_netcdf: full/path/to/lnd_cull_mesh.nc
#   sFilename_mesh_boundary: full/path/to/boundary_wgs.geojson
# 	sFilename_basins: full/path/to/pyflowline_susquehanna_basins.json.

# The json file will be overwritten, you may want to make a copy of it first.

from pyflowline.configuration.change_json_key_value import change_json_key_value as pyflowline_change_json_key_value

# Set the path to the output folder
# Pass the configuration filename followed by a single key-value pair.
pyflowline_change_json_key_value(
    oFilename_domain_config,
    'sWorkspace_output', str(oFolder_output))

# Set the path to the mpas mesh file
pyflowline_change_json_key_value(
    oFilename_domain_config,
    'sFilename_mesh_netcdf', str(oFilename_mesh_netcdf))

# Set the path to the boundary file used to clip the mesh
pyflowline_change_json_key_value(
    oFilename_domain_config,
    'sFilename_mesh_boundary', str(oFilename_mesh_boundary))

# Set the path to the basin ("child") configuration file
pyflowline_change_json_key_value(
    oFilename_domain_config,
    'sFilename_basins', str(oFilename_basins_config))

# These parameters are not strictly necessary, but will reduce the potential for errors or confusion.
pyflowline_change_json_key_value(
    oFilename_domain_config,
    'sFilename_model_configuration', str(oFilename_domain_config))

pyflowline_change_json_key_value(
    oFilename_domain_config,
    'sWorkspace_data', str(oPath_parent.joinpath('data')))

#%% Update the basin ("child") configuration file

# Set the path to the flowline file.
pyflowline_change_json_key_value(
    oFilename_basins_config,
    'sFilename_flowline_filter', str(oFilename_flowline),
    iFlag_basin_in=1) # Set iFlag_basin_in=1 when changing the basin configuration file.

#%% Read the configuration file
from pyflowline.configuration.read_configuration_file import pyflowline_read_configuration_file

oPyflowline = pyflowline_read_configuration_file(
    oFilename_domain_config,
    iCase_index_in=iCase_index,
    sDate_in=sDate)

#%% Check the model parameters
oPyflowline.pyflowline_print()

#%% Now we can change some model parameters

# There are two ways to change the model parameters: 1) use a function, or 2) assign a value directly. Use the function:
oPyflowline.pyflowline_change_model_parameter(
    'sWorkspace_output', str(oFolder_output))

# To change a parameter for a basin instead of the whole model domain, use the iFlag_basin_in option, this will change the value for all of the basins in the basin configuration file.
oPyflowline.pyflowline_change_model_parameter(
    'dLatitude_outlet_degree', 39.462000,
    iFlag_basin_in=1)

oPyflowline.pyflowline_change_model_parameter(
    'dLongitude_outlet_degree', -76.009300,
    iFlag_basin_in=1)

# The second way is to assign a value directly
oPyflowline.aBasin[0].dLatitude_outlet_degree=39.462000
oPyflowline.aBasin[0].dLongitude_outlet_degree=-76.009300

#%% Export the config file after changing parameters

# If desired, the config file can be exported to disk after changing parameters. By default, the configuration files are written to the output folder.
# oPyflowline.pyflowline_export_config_to_json()

#%% Now we can build the flowline 

# Set up the flowlinecase and optionally plot the flowline.
oPyflowline.pyflowline_setup()

if iFlag_visualization == 1:
    oPyflowline.plot(sVariable_in='flowline_filter', sFilename_output_in='filter_flowline.png')
    pass

#%% Run the flowline simplification step.

# [see](https://pyflowline.readthedocs.io/en/latest/algorithm/algorithm.html#flowline-simplification)
if iFlag_simulation == 1:
    oPyflowline.pyflowline_flowline_simplification();
    pass

#%%
if iFlag_visualization == 1:

    # aExtent_meander = [-76.5, -76.2, 41.6, 41.9]
    # oPyflowline.plot(sVariable_in='flowline_simplified', sFilename_output_in='flowline_simplified.png')
    # oPyflowline.plot(sVariable_in='flowline_simplified',
    #                  sFilename_output_in='flowline_simplified_zoom.png',
    #                  aExtent_in=aExtent_meander)
    pass

#%%
if iFlag_simulation == 1:
    aCell = oPyflowline.pyflowline_mesh_generation()

#%%
if iFlag_visualization == 1:
    # This will fail unless iFlag_simulation is true or has been run
    # oPyflowline.plot(sVariable_in='mesh', sFilename_output_in='mesh.png' )
    pass
#%%
if iFlag_simulation == 1:
    oPyflowline.pyflowline_reconstruct_topological_relationship(aCell)
    # oPyflowline.pyflowline_export()

#%%
if iFlag_visualization == 1:
    # oPyflowline.plot(sVariable_in='overlap', sFilename_output_in='mesh_w_flowline.png')
    pass

# %%
print('Finished')
print('Try opening the mesh and flowline files in QGIS for visualization.')
