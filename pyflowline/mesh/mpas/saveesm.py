
import subprocess
import os
import time
import numpy as np
import xarray

from geometric_features import GeometricFeatures

from from pyflowline.mesh.mpas.mpasmsh import jigsaw_mesh_to_netcdf, inject_edge_tags, subtract_critical_passages, mask_reachable_ocean

from mpas_tools.mesh.conversion import convert, mask, cull
from mpas_tools.io import write_netcdf
from mpas_tools.ocean.inject_preserve_floodplain import    inject_preserve_floodplain
from mpas_tools.ocean.coastline_alteration import    widen_transect_edge_masks, add_critical_land_blockages,    add_land_locked_cells_to_mask

HERE = os.path.abspath(os.path.dirname(__file__))

def saveesm(path, geom, mesh,
            preserve_floodplain=False,
            floodplain_elevation=20.0,
            do_inject_elevation=False,
            with_cavities=False,
            lat_threshold=43.00,
            with_critical_passages=True):
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
        jigsaw_mesh_to_netcdf(
            mesh=mesh,
            on_sphere=True,
            sphere_radius=np.mean(geom.radii) * 1e3,
            output_name=os.path.join(
                path, "tmp", "mesh_triangles.nc"))

    if (geom.mshID.lower() == "euclidean-mesh"):
        print("Forming mesh_triangles.nc")
        jigsaw_mesh_to_netcdf(
            mesh=mesh,
            on_sphere=False,
            output_name=os.path.join(
                path, "tmp", "mesh_triangles.nc"))

    print("Forming base_mesh.nc")
    write_netcdf(
        convert(xarray.open_dataset(
            os.path.join(
                path, "tmp", "mesh_triangles.nc"))),
        fileName=os.path.join(
            path, "out", "base_mesh.nc"))

    """
    if do_inject_elevation:
        print("Injecting cell elevations")
        inject_elevation(
            cell_elev=mesh.value,
            mesh_file=os.path.join(
                path, "out", "base_mesh.nc"))
    """

    if preserve_floodplain:
        print("Injecting floodplain flag")
        inject_preserve_floodplain(
            mesh_file=os.path.join(
                path, "out", "base_mesh.nc"),
            floodplain_elevation=floodplain_elevation)

    args = ["paraview_vtk_field_extractor.py",
            "--ignore_time",
            "-l",
            "-d", "maxEdges=0",
            "-v", "allOnCells",
            "-f", os.path.join(
                path, "out", "base_mesh.nc"),
            "-o", os.path.join(
                path, "out", "base_mesh_vtk")]
    print("")
    print("running:", " ".join(args))
    subprocess.check_call(args, env=os.environ.copy())

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

    # remove the region south of 60S so we can replace
    # it based on ice-sheet topography
    fcSouthMask = gf.read(
        componentName="ocean", objectType="region",
        featureNames=["Global Ocean 90S to 60S"])

    fcLandCoverage = \
        fcLandCoverage.difference(fcSouthMask)

    # add land coverage from either the full ice sheet
    # or just the grounded part
    if with_cavities:
        fcAntarcticLand = gf.read(
            componentName="bedmap2", objectType="region",
            featureNames=["AntarcticGroundedIceCoverage"])
    else:
        fcAntarcticLand = gf.read(
            componentName="bedmap2", objectType="region",
            featureNames=["AntarcticIceCoverage"])

    fcLandCoverage.merge(fcAntarcticLand)

    # save the feature collection to a geojson file
    fcLandCoverage.to_geojson(
        os.path.join(
            path, "tmp", "land_coverage.geojson"))

    # Create the land mask based on the land coverage,
    # i.e. coastline data.
    dsBaseMesh = xarray.open_dataset(
        os.path.join(path, "out", "base_mesh.nc"))
    dsLandMask = mask(dsBaseMesh, fcMask=fcLandCoverage)

    dsLandMask = add_land_locked_cells_to_mask(
        dsLandMask, dsBaseMesh,
        latitude_threshold=lat_threshold, nSweeps=20)

    if with_critical_passages:
        # merge transects for critical passages into
        # critical_passages.geojson
        fcCritPassages = gf.read(
            componentName="ocean", objectType="transect",
            tags=["Critical_Passage"])

        # create masks from the transects
        dsCritPassMask = \
            mask(dsBaseMesh, fcMask=fcCritPassages)

        # Alter critical passages to be at least two
        # cells wide, to avoid sea ice blockage.
        dsCritPassMask = widen_transect_edge_masks(
            dsCritPassMask, dsBaseMesh,
            latitude_threshold=lat_threshold)

        dsLandMask = subtract_critical_passages(
            dsLandMask, dsCritPassMask)

        # merge transects for critical land blockages
        # into critical_land_blockages.geojson
        fcCritBlockages = gf.read(
            componentName="ocean", objectType="transect",
            tags=["Critical_Land_Blockage"])

        # create masks from the transects for critical
        # land blockages
        dsCritBlockMask = \
            mask(dsBaseMesh, fcMask=fcCritBlockages)

        dsLandMask = add_critical_land_blockages(
            dsLandMask, dsCritBlockMask)

    # create seed points for a flood fill of the ocean
    # use all points in the ocean directory, on the
    # assumption that they are, in fact *in* the ocean
    fcSeed = gf.read(
        componentName="ocean",
        objectType="point", tags=["seed_point"])

    # update the land mask to ensure all ocean cells really
    # are "reachable" from the rest of the global ocean
    dsLandMask = mask_reachable_ocean(
        dsMesh=dsBaseMesh,
        dsMask=dsLandMask, fcSeed=fcSeed)

    # cull the (ocean) mesh based on the land mask, and a
    # cull the (land) mesh using the inverse mask

    if preserve_floodplain:
    # with "preserve_floodplains", the (ocean) mesh will
    # contain overlap with the (land) mesh, otherwise the
    # two are "perfectly" disjoint
        dsCulledMesh = cull(
            dsBaseMesh, dsMask=dsLandMask,
            dsPreserve=dsBaseMesh,
            graphInfoFileName=os.path.join(
                path, "out", "culled_graph.info"))

        dsInvertMesh = cull(
            dsBaseMesh, dsInverse=dsLandMask,
            graphInfoFileName=os.path.join(
                path, "out", "invert_graph.info"))

    else:
        dsCulledMesh = cull(
            dsBaseMesh, dsMask=dsLandMask,
            graphInfoFileName=os.path.join(
                path, "out", "culled_graph.info"))

        dsInvertMesh = cull(
            dsBaseMesh, dsInverse=dsLandMask,
            graphInfoFileName=os.path.join(
                path, "out", "invert_graph.info"))

    write_netcdf(
        dsCulledMesh, os.path.join(
            path, "out", "culled_mesh.nc"), netcdfFormat)

    write_netcdf(
        dsInvertMesh, os.path.join(
            path, "out", "invert_mesh.nc"), netcdfFormat)

    args = ["paraview_vtk_field_extractor.py",
            "--ignore_time",
            "-d", "maxEdges=",
            "-v", "allOnCells",
            "-f", os.path.join(
                path, "out", "culled_mesh.nc"),
            "-o", os.path.join(
                path, "out", "culled_mesh_vtk")]
    print("")
    print("running", " ".join(args))
    subprocess.check_call(args, env=os.environ.copy())

    args = ["paraview_vtk_field_extractor.py",
            "--ignore_time",
            "-d", "maxEdges=",
            "-v", "allOnCells",
            "-f", os.path.join(
                path, "out", "invert_mesh.nc"),
            "-o", os.path.join(
                path, "out", "invert_mesh_vtk")]
    print("running", " ".join(args))
    subprocess.check_call(args, env=os.environ.copy())

    ttoc = time.time()

    print("CPUSEC =", (ttoc - ttic))

    return


def inject_elevation(mesh_file):

    return
