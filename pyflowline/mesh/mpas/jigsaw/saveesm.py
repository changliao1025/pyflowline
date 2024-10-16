
import subprocess
import os
import time
import numpy as np
import xarray as xr
from geometric_features import GeometricFeatures
from pyearth.system.python.retrieve_python_environment import retrieve_python_environment
from pyflowline.mesh.mpas.mpas_tools.mpasmsh import jigsaw_mesh_to_netcdf, inject_edge_tags, subtract_critical_passages, mask_reachable_ocean

import mpas_tools
from mpas_tools.mesh.conversion import convert, cull
from mpas_tools.logging import check_call
from mpas_tools.io import write_netcdf

HERE = os.path.abspath(os.path.dirname(__file__))

def saveesm(path, geom, mesh):
    """
    SAVEESM: export a jigsaw mesh obj. to MPAS-style output.

    1. Writes "mesh_triangles.nc" and "base_mesh.nc" files.
    2. (Optionally) injects elevation + floodplain data.
    3. Calls MPAS-Tools + Geometric-Data to cull mesh into
       ocean/land partitions.
    4. Writes "culled_mesh.nc" (ocean) and "invert_mesh.nc"
       (land) MPAS-spec. output files.

    Data is written to "../path/out/" and/or "../path/tmp/".

    """
    # Authors: Darren Engwirda

    ttic = time.time()

    print("")
    print("Running MPAS mesh-tools...")

    inject_edge_tags(mesh)

    # adapted from BUILD_MESH.py

    if (geom.mshID.lower() == "ellipsoid-mesh"):
        print("Forming mesh_triangles.nc")
        sFilename_triangles = os.path.join(path, "tmp", "mesh_triangles.nc")
        jigsaw_mesh_to_netcdf(
            mesh=mesh,
            on_sphere=True,
            sphere_radius=np.mean(geom.radii) * 1e3,
            output_name=sFilename_triangles)

    if (geom.mshID.lower() == "euclidean-mesh"):
        print("Forming mesh_triangles.nc")
        sFilename_triangles = os.path.join(path, "tmp", "mesh_triangles.nc")
        jigsaw_mesh_to_netcdf(
            mesh=mesh,
            on_sphere=False,
            output_name=sFilename_triangles)

    print("Forming base_mesh.nc")
    sFilename_base_mesh = os.path.join(path, "out", "base_mesh.nc")
    write_netcdf( convert(xr.open_dataset(
            sFilename_triangles)),
            fileName=sFilename_base_mesh)

    sFilename_executable = 'paraview_vtk_field_extractor.py'
    #find the full path to the python executable
    for folder in os.environ['PATH'].split(os.pathsep):
        sFilename_dummy = os.path.join(folder, sFilename_executable)
        if os.path.isfile(sFilename_dummy):
            iFlag_found_binary = 1
            sPath_executable = sFilename_dummy
            break

    sFilename_base_mesh_vtk = os.path.join(path, "out", "base_mesh_vtk")
    args = [sPath_executable,
            "--ignore_time",
            "-l",
            "-d", "maxEdges=0",
            "-v", "allOnCells",
            "-f", sFilename_base_mesh,
            "-o", sFilename_base_mesh_vtk]
    print("")
    print("running:", " ".join(args))

    sConda_env_path , sConda_env_name = retrieve_python_environment()
    sPython = sConda_env_path + "/bin/python3"
    args = [sPython] + args
    sEnv = os.environ.copy()
    subprocess.check_call(args, env=sEnv)

    # adapted from CULL_MESH.py

    # required for compatibility with MPAS
    netcdfFormat = "NETCDF3_64BIT"

    gf = GeometricFeatures(
        cacheLocation="{}".format(os.path.join(
            HERE, "..", "data", "geometric_data")))

    # start with the land coverage from Natural Earth
    fcLandCoverage = gf.read(
        componentName="natural_earth", objectType="region",
        featureNames=["Land Coverage"])

    # save the feature collection to a geojson file
    sFilename_land_coverage = os.path.join( path, "tmp", "land_coverage.geojson")
    fcLandCoverage.to_geojson(sFilename_land_coverage)

    # Create the land mask based on the land coverage,
    # i.e. coastline data.
    dsBaseMesh = xr.open_dataset(  sFilename_base_mesh)
    sFilename_land_mask = os.path.join( path, "tmp", "land_mask.nc")

    _land_mask_from_geojson(    mesh_filename=sFilename_base_mesh,
                                geojson_filename=sFilename_land_coverage,
                                mask_filename=sFilename_land_mask)
    dsLandMask = xr.open_dataset(sFilename_land_mask)

    sFilename_culled_mesh = os.path.join(path, "out", "culled_graph.info")

    dsCulledMesh = cull(
        dsBaseMesh, dsMask=dsLandMask,
        graphInfoFileName=sFilename_culled_mesh)

    sFilename_invert_mesh = os.path.join(path, "out", "invert_graph.info")
    dsInvertMesh = cull(
        dsBaseMesh, dsInverse=dsLandMask,
        graphInfoFileName= sFilename_invert_mesh)

    sFilename_culled_mesh = os.path.join(path, "out", "culled_mesh.nc")
    write_netcdf(
        dsCulledMesh, sFilename_culled_mesh, netcdfFormat)

    sFilename_invert_mesh = os.path.join(path, "out", "invert_mesh.nc")
    write_netcdf(
        dsInvertMesh, sFilename_invert_mesh, netcdfFormat)

    sFilename_culled_mesh_vtk = os.path.join(path, "out", "culled_mesh_vtk")
    args = [sPath_executable,
            "--ignore_time",
            "-d", "maxEdges=",
            "-v", "allOnCells",
            "-f", sFilename_culled_mesh,
            "-o", sFilename_culled_mesh_vtk]

    print("")
    print("running", " ".join(args))

    args = [sPython] + args
    subprocess.check_call(args, env=os.environ.copy())

    sFilename_invert_mesh_vtk = os.path.join(path, "out", "invert_mesh_vtk")
    args = [sPath_executable,
            "--ignore_time",
            "-d", "maxEdges=",
            "-v", "allOnCells",
            "-f", sFilename_invert_mesh,
            "-o", sFilename_invert_mesh_vtk]
    print("running", " ".join(args))
    args = [sPython] + args
    subprocess.check_call(args, env=os.environ.copy())

    ttoc = time.time()

    print("CPUSEC =", (ttoc - ttic))

    return sFilename_culled_mesh, sFilename_invert_mesh



def _land_mask_from_geojson(    mesh_filename, geojson_filename, mask_filename):
    gf = GeometricFeatures()

    # start with the land coverage from Natural Earth
    fcLandCoverage = gf.read(componentName='natural_earth',
                             objectType='region',
                             featureNames=['Land Coverage'])

    # remove the region south of 60S so we can replace it based on ice-sheet
    # topography
    fcSouthMask = gf.read(componentName='ocean', objectType='region',
                          featureNames=['Global Ocean 90S to 60S'])

    fcLandCoverage = fcLandCoverage.difference(fcSouthMask)

    # Add "land" coverage from either the full ice sheet or just the grounded
    # part

    fcAntarcticLand = gf.read(
            componentName='bedmachine', objectType='region',
            featureNames=['AntarcticIceCoverage'])

    fcLandCoverage.merge(fcAntarcticLand)

    # save the feature collection to a geojson file
    fcLandCoverage.to_geojson(geojson_filename)

    # these defaults may have been updated from config options -- pass them
    # along to the subprocess
    netcdf_format = mpas_tools.io.default_format
    netcdf_engine = 'scipy'

    # Create the land mask based on the land coverage, i.e. coastline data
    process_count = 1
    args = ['compute_mpas_region_masks',
            '-m', mesh_filename,
            '-g', geojson_filename,
            '-o', mask_filename,
            '-t', 'cell',
            '--process_count', f'{process_count}',
            '--format', netcdf_format,
            '--engine', netcdf_engine]
    check_call(args, logger=None)

