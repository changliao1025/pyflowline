#####################
Quickstart
#####################

Users can run a PyFlowline simulation in the following steps:

1. Create a new Python environment using Conda, and activate the new environment.
2. Install the package using `conda install -c conda-forge pyflowline`. Conda will automatically install all the required dependencies.
3. Clone the latest PyFlowline repository from https://github.com/changliao1025/pyflowline. 
4. Download the additional large MPAS mesh file `lnd_cull_mesh.nc` from https://github.com/changliao1025/pyflowline/releases/tag/0.2.0 and move it under the `data/susquehanna/input` folder.
5. Open the `examples/susquehanna/pyflowline_susquehanna_mpas.json` file and change `sWorkspace_output` to the full path to the directory where you want to save the output (e.g. `/full/path/to/pyflowline/data/susquehanna/output`), change `"sFilename_mesh_netcdf"` to the full path to `lnd_cull_mesh.nc`, `"sFilename_mesh_boundary"` to the full path to `data/susquehanna/input/mesh_boundary_buffer.geojson`, and `"sFilename_basins"` to the full path to `examples/susquehanna/pyflowline_susquehanna_basins.json`.
6. Open the `examples/susquehanna/pyflowline_susquehanna_basins.json` file and change `"sFilename_flowline_filter"` to the full path to `data/susquehanna/input/flowline.geojson`. Ignore the other settings in these json files for now.
7. Open the preferred Python IDE (Visual Studio Code recommended) and run the  `examples/susquehanna/run_simulation_mpas.py` Python script. Optionally, you can also run the `notebooks/mpas_example.ipynb` notebook.
8. You should produce a list of model outputs in the `data/susquehanna/output` folder or the user-specified output folder.

If you encounter any issues, refer to the FAQ or submit a GitHub issue (https://github.com/changliao1025/pyflowline/issues).