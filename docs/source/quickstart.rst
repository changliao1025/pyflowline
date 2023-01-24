#####################
Quick Start
#####################

Users can run a PyFlowline simulation in the following steps:

1. Create a new Python environment using Conda, and install the following packages: 
    * `numpy`
    * `gdal`
    * `netCDF4`
    * `shapely`
    * `cython` (optional)
    * `cartopy` (optional, for visualization)
    * `matplotlib` (optional, for visualization)
2. Clone the latest PyFlowline repository from https://github.com/changliao1025/pyflowline. Or Install the PyFlowline through Conda for a released version.
3. Download the additional large files (DEM and MPAS mesh) and move them under the `data/susquehanna/input` folder.
4. Change the `sFilename_mesh_netcdf` and `sFilename_basins` to the actual paths,
5. Open the preferred Python IDE and run the  `examples/susquehanna/run_simulation_mpas.py` Python script. Optionally, you can also run the `notebooks/pyflowline.ipynb` notebook.
6. You should produce a list of model outputs in the `data/susquehanna/output` folder.

If you encountered any issues, refer to the FAQ or submit a Github issue (https://github.com/changliao1025/pyflowline/issues).