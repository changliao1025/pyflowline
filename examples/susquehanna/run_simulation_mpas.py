import os, sys
from pathlib import Path
from os.path import realpath

from pyflowline.configuration.change_json_key_value import pyflowline_change_json_key_value
from pyflowline.configuration.read_configuration_file import pyflowline_read_configuration_file

#%% Set up workspace path
sPath_parent = str(Path(__file__).parents[2]) # data is located two dir's up
sys.path.append(sPath_parent)

#%% Define the case information
sDomainName = 'susquehanna'
iCase_index = 1
iFlag_simulation = 0
iFlag_visualization = 1
sMesh = 'mpas'
sDate = '20230701'

sFilename_domain_configuration = 'pyflowline_susquehanna_mpas.json'
sFilename_basins_configuration = 'pyflowline_susquehanna_basins.json'
sFilename_flowline = 'flowline.geojson'
sFilename_mesh_netcdf = 'lnd_cull_mesh.nc'
sFilename_mesh_boundary = 'boundary_wgs.geojson'

#%% Define the configuration file parameters (paths to inputs and outputs)

# Define the full path to the input and output folders.
sFolder_input = realpath(os.path.join(
    sPath_parent, 'data', sDomainName, 'input'))
sFolder_output = realpath(os.path.join(
    sPath_parent, 'data', sDomainName, 'output'))

# Define the full path to the domain ("parent") configuration file.
sFilename_domain_configuration = realpath(os.path.join(
    sFolder_input, sFilename_domain_configuration))

# Define the full path to the individual basin ("child") configuration file.
sFilename_basins_configuration = realpath(os.path.join(
    sFolder_input, sFilename_basins_configuration))

# Define the full path to the MPAS mesh file
sFilename_mesh_netcdf = realpath(os.path.join(
    sFolder_input, sFilename_mesh_netcdf))

# Define the full path to the input flowline file.
sFilename_flowline = realpath(os.path.join(
    sFolder_input, sFilename_flowline))

# Define the full path to the boundary file used to clip the mesh.
sFilename_mesh_boundary = realpath(os.path.join(
    sFolder_input, sFilename_mesh_boundary))

# Confirm that the domain configuration file exists.
if os.path.isfile(sFilename_domain_configuration):
    pass
else:
    print('The domain configuration file does not exist: ', sFilename_domain_configuration)

#%% Update the domain ("parent") configuration file

# As noted in the pyflowline docs, some parameters in the configuration file need to be set before we can create the flowline object. The following parameters need to be set in the domain configuration file (pyflowline_susquehanna_mpas.json):
# 
#   sWorkspace_output: full/path/to/output
#   sFilename_mesh_netcdf: full/path/to/lnd_cull_mesh.nc
#   sFilename_mesh_boundary: full/path/to/boundary_wgs.geojson
# 	sFilename_basins: full/path/to/pyflowline_susquehanna_basins.json.
#
# The json file will be overwritten, you may want to make a copy of it first.

# The first argument is the configuration file name, followed by one key-value pair, and then the iFlag_basin_in flag. If iFlag_basin_in is None, then the function updates the "parent" (domain) configuration file. Otherwise, it updates the "child" (basin) configuration file.

# Set the path to the output folder
pyflowline_change_json_key_value(
    sFilename_domain_configuration, 
    'sWorkspace_output', sFolder_output)

# Set the path to the mpas file we just downloaded
pyflowline_change_json_key_value(
    sFilename_domain_configuration, 
    'sFilename_mesh_netcdf', sFilename_mesh_netcdf) 

# Set the path to the boundary file used to clip the mesh
pyflowline_change_json_key_value(
    sFilename_domain_configuration, 
    'sFilename_mesh_boundary', sFilename_mesh_boundary) 

# Set the path to the individual-basin ("child") configuration file
pyflowline_change_json_key_value(
    sFilename_domain_configuration, 
    'sFilename_basins', sFilename_basins_configuration) 

#%% Update the basin configuration file

# For the basin configuration file (pyflowline_susquehanna_basins.json), the following parameters need to be set:
# 
# 	sFilename_flowline_filter: full/path/to/flowline.geojson

# Set the path to the user-provided flowline. Note that when changing the basin ("child") configuration file, set iFlag_basin_in=1.
pyflowline_change_json_key_value(
    sFilename_basins_configuration, 
    'sFilename_flowline_filter', sFilename_flowline, 
    iFlag_basin_in=1)

#%% Read the configuration file
oPyflowline = pyflowline_read_configuration_file(
    sFilename_domain_configuration, 
    iCase_index_in=iCase_index, 
    sDate_in=sDate)

# Check the model parameters
oPyflowline.pyflowline_print()

#%% Now we can change some model parameters

# There are two ways to change the model parameters: 1) use a function, or 2) assign a value directly.
oPyflowline.pyflowline_change_model_parameter(
    'sWorkspace_output', sFolder_output)

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

#%%
oPyflowline.pyflowline_setup()

#%%
if iFlag_visualization == 1:
    oPyflowline.plot(sVariable_in='flowline_filter', sFilename_output_in='filter_flowline.png')
    pass

#%%
if iFlag_simulation == 1:
    oPyflowline.pyflowline_flowline_simplification()
    pass

#%%
if iFlag_visualization == 1:
     
    aExtent_meander = [-76.5, -76.2, 41.6, 41.9]
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

#%%
if iFlag_visualization == 1:
    # oPyflowline.plot(sVariable_in='overlap', sFilename_output_in='mesh_w_flowline.png')
    pass

# oPyflowline.pyflowline_export()

print('Finished')

# %%
