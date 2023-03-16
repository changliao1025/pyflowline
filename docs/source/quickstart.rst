#####################
Quickstart
#####################

Users can run a PyFlowline simulation in the following steps:

1. Create a new Python environment using Conda, and activate new environment
2. Clone the latest PyFlowline repository from https://github.com/changliao1025/pyflowline. 
3. Install the dependency packages using conda.
4. Download the additional large files (MPAS mesh (https://github.com/changliao1025/pyflowline/releases/tag/0.2.0)) and move them under the `data/susquehanna/input` folder.
5. Change the `sFilename_mesh_netcdf`, `sFilename_basins`, and `sFilename_flowline_filter` to the actual paths,
6. Open the preferred Python IDE and run the  `examples/susquehanna/run_simulation_mpas.py` Python script. Optionally, you can also run the `notebooks/pyflowline.ipynb` notebook.
7. You should produce a list of model outputs in the `data/susquehanna/output` folder.

If you encounter any issues, refer to the FAQ or submit a GitHub issue (https://github.com/changliao1025/pyflowline/issues).