
import subprocess
import os
import time

import numpy as np
import xarray as xr
from pyearth.system.python.get_python_environment import get_python_environment
from pyflowline.mesh.mpas.mpas_tools.mpasmsh import jigsaw_mesh_to_netcdf, inject_edge_tags

from mpas_tools.mesh.conversion import convert
from mpas_tools.io import write_netcdf
from jigsawpy import savevtk


def savetin(sWorkspace_jigsaw_out, geom, mesh):

    ttic = time.time()

    print("")
    print("Convert jigsaw mesh to netcdf...")
    inject_edge_tags(mesh)
    # adapted from BUILD_MESH.py
    sFilename_triangles = os.path.join(sWorkspace_jigsaw_out, "tmp", "mesh_triangles.nc")
    if (geom.mshID.lower() == "ellipsoid-mesh"):
        print("Forming mesh_triangles.nc")
        jigsaw_mesh_to_netcdf(
            mesh=mesh,
            on_sphere=True,
            sphere_radius=np.mean(geom.radii) * 1e3,
            output_name=sFilename_triangles)

    print("Forming base_mesh.nc")
    sFilename_base_mesh = os.path.join(sWorkspace_jigsaw_out, "out", "base_mesh.nc")
    write_netcdf( convert(xr.open_dataset(
            sFilename_triangles)),
            fileName=sFilename_base_mesh)

    sFilename_executable = 'paraview_vtk_field_extractor.py'
    #find the full sWorkspace_jigsaw_out to the python executable
    for folder in os.environ['PATH'].split(os.pathsep):
        sFilename_dummy = os.path.join(folder, sFilename_executable)
        if os.path.isfile(sFilename_dummy):
            iFlag_found_binary = 1
            sPath_executable = sFilename_dummy
            break

    sFilename_base_mesh_vtk = os.path.join(sWorkspace_jigsaw_out, "out", "base_mesh_vtk")
    args = [sPath_executable,
            "--ignore_time",
            "-l",
            "-d", "maxEdges=0",
            "-v", "allOnCells",
            "-f", sFilename_base_mesh,
            "-o", sFilename_base_mesh_vtk]
    print("")
    print("running:", " ".join(args))


    sConda_env_path , sConda_env_name = get_python_environment()
    sPython = sConda_env_path + "/bin/python3"
    args = [sPython] + args
    sEnv = os.environ.copy()
    subprocess.check_call(args, env=sEnv)

    ttoc = time.time()

    print("CPUSEC =", (ttoc - ttic))

    return  sFilename_triangles