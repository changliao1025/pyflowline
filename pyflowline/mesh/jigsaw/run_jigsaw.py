import os
from pathlib import Path
import shutil
import math
import copy
import time
import numpy as np
from osgeo import gdal
import jigsawpy
from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import gdal_read_geotiff_file

import mpas_tools.mesh.creation.mesh_definition_tools as mdt
from pyflowline.mesh.jigsaw.inpoly2 import inpoly2
from pyflowline.mesh.jigsaw.loadgeo import loadgeo
from pyflowline.mesh.jigsaw.utility import addpoint, addline, addpoly, innerto
#from pyflowline.classes.vertex import pyvertex

def compute_mask(aData, value, dLongitude_upper_left, dLatitude_upper_left, pixelWidth, pixelHeight, ncolumn_space, nrow_space):
    dummy_index = np.where(aData == value)
    irows, icols = dummy_index
    dLon = dLongitude_upper_left + icols * pixelWidth
    dLat = dLatitude_upper_left + irows * pixelHeight
    iX = ((dLon + 180.0) / (360.0 / ncolumn_space)).astype(int)
    iY = ((dLat + 90) / (180.0 / nrow_space)).astype(int)
    # Ensure indices are within range
    iX = np.clip(iX, 0, ncolumn_space - 1)
    iY = np.clip(iY, 0, nrow_space - 1)
    # Create a mask with the same shape as the grid
    mask = np.full((nrow_space, ncolumn_space), False, dtype=bool)
    # Use 2D indexing to set the mask
    mask[iY, iX] = True
    return mask

def coarsen_mask(mask, down):
#-- down-sample the mask by combining every set of DOWN pixels into one
    rows = mask.shape[0] // down
    cols = mask.shape[1] // down

    mtmp = np.full((rows, cols), False, dtype=bool)
    for jpos in range(down):
        for ipos in range(down):
            iend = mask.shape[0] - down + ipos + 1
            jend = mask.shape[1] - down + jpos + 1
            mtmp = np.logical_or(
                mtmp,
            mask[ipos:iend:down, jpos:jend:down])

    return mtmp

def retrieve_coordinates_from_array(aData_in, dValue_in, pixelWidth, pixelHeight):
    dummy_index = np.where(aData_in == dValue_in)
    irows, icols = dummy_index
    dLon = -180.0 + icols * pixelWidth
    dLat = 90.0 + irows * pixelHeight
    #clip to the range
    dLon = np.clip(dLon, -180.0, 180.0)
    dLat = np.clip(dLat, -90.0, 90.0)
    return dLon, dLat

#def run_jigsaw_mpas_workflow(sWorkspace_jigsaw_in,
#                             aConfig_in = None,
#                             projector=[0.0, 0.0],
#                             aFilename_river_in=None,
#                             aFilename_watershed_boundary_in= None,
#                             aFilename_lake_boundary_in = None,
#                             aFilename_coastline_in = None):
#    """
#    run_jigsaw_mpas_workflow: main call to the ICoM mesh-gen. infrastructure.
#
#    1. Call user-defined MESH-PATH/COMPOSE.py to assemble
#       mesh geometry, spacing pattern, initial conditions
#       and mesh-gen. parameters.
#    2. Call JIGSAW to build the triangulation.
#    4. Call MPAS meshtools to make the MPAS mesh data files.
#
#    """
#    # Authors: Darren Engwirda and Chang Liao
#
#
#
##---------- call JIGSAW to build the initial triangular mesh
#
#    geom, gprj, mesh, mprj =  runjgsw(sWorkspace_jigsaw_in, projector,
#                                      aConfig_in=aConfig_in,
#                                      aFilename_river_in=aFilename_river_in,
#                                      aFilename_watershed_boundary_in= aFilename_watershed_boundary_in,
#                                      aFilename_lake_boundary_in = aFilename_lake_boundary_in,
#                                      aFilename_coastline_in = aFilename_coastline_in)
#
#
#
##-------------------------------------- write output for ESM
#
#    iFlag_mpas_tool = 1
#    if iFlag_mpas_tool ==1:
#        from pyflowline.mesh.jigsaw.saveesm import saveesm
#        sFilename_culled_mesh, sFilename_invert_mesh = saveesm(sWorkspace_jigsaw_in, geom, mesh)
#    else:
#        #we will a new function to convert jigsaw mesh to mpas mesh11
#        print('The algorithm is not completed yet')
#        pass
#
#    return sFilename_invert_mesh

def run_jigsaw(sWorkspace_jigsaw_in, projector,
            aConfig_in=None,
            aFilename_river_network_in=None,
            aFilename_watershed_boundary_in= None,
            aFilename_lake_boundary_in = None,
            aFilename_coastline_in = None):
    """
    RUNJGSW: main call to JIGSAW to build the triangulation.

    MESH-PATH should point to a user-defined mesh directory,
    containing the COMPOSE.py template.

    Firstly, MESH-PATH/COMPOSE.py is called to build
    user-defined geometry, initial conditions, mesh spacing
    and mesh configuration information.

    The boolean flags MAKE-BOOL control whether mesh
    information is built from scratch, or if an exitsing an
    existing file is to be used. For example, setting
    MAKE-BOOL.SPAC = FALSE relies on an existing spacing
    pattern be available in MESH-PATH/tmp/.

    This information is written to MESH-PATH/tmp/ to be
    accessed by subsequent calls to JIGSAW.

    Finally, JIGSAW is run to build the triangular mesh,
    calling either the multi-level (TETRIS) or single-level
    (JIGSAW) algorithms. Cells are assigned ID-tags via
    the polygon/regions defined in GEOM.BOUNDS.

    Returns full-dimensional and 2d-projected msh_t objects.

    """
    # Authors: Darren Engwirda

    mesh = jigsawpy.jigsaw_msh_t()
    mprj = jigsawpy.jigsaw_msh_t()
    gprj = jigsawpy.jigsaw_msh_t()

#------------------------------------ setup via user COMPOSE
    FULL_SPHERE_RADIUS = +6371.0 #unit check
    if aConfig_in is not None:
        if "ncolumn_space" in aConfig_in:
            ncolumn_space = aConfig_in["ncolumn_space"]
        else:
            ncolumn_space = 360 #1 degree spacing

        if "nrow_space" in aConfig_in:
            nrow_space = aConfig_in["nrow_space"]
        else:
            nrow_space = 180 #1 degree spacing

        if "dSpac_value" in aConfig_in:
            dSpac_value = aConfig_in["dSpac_value"]
        else:
            dSpac_value = 50.0

        if "iFlag_geom" in aConfig_in:
            iFlag_geom = aConfig_in["iFlag_geom"]
        else:
            iFlag_geom = False

        if "iFlag_spac" in aConfig_in:
            iFlag_spac = aConfig_in["iFlag_spac"]
        else:
            iFlag_spac = False

        if "iFlag_init" in aConfig_in:
            iFlag_init = aConfig_in["iFlag_init"]
        else:
            iFlag_init = False

        if "iFlag_opts" in aConfig_in:
            iFlag_opts = aConfig_in["iFlag_opts"]
        else:
            iFlag_opts = False

        if "iFlag_ocean" in aConfig_in:
            iFlag_ocean = aConfig_in["iFlag_ocean"]
        else:
            iFlag_ocean = False

        if "iFlag_land" in aConfig_in:
            iFlag_land = aConfig_in["iFlag_land"]
        else:
            iFlag_land = False

        if "iFlag_geom_dam" in aConfig_in:
            iFlag_geom_dam = aConfig_in["iFlag_geom_dam"]
        else:
            iFlag_geom_dam = False

        if "iFlag_spac_dam" in aConfig_in:
            iFlag_spac_dam = aConfig_in["iFlag_spac_dam"]
        else:
            iFlag_spac_dam = False

        if "iFlag_geom_city" in aConfig_in:
            iFlag_geom_city = aConfig_in["iFlag_geom_city"]
        else:
            iFlag_geom_city = False

        if "iFlag_spac_city" in aConfig_in:
            iFlag_spac_city = aConfig_in["iFlag_spac_city"]
        else:
            iFlag_spac_city = False

        if "iFlag_geom_river_network" in aConfig_in:
            iFlag_geom_river_network = aConfig_in["iFlag_geom_river_network"]
        else:
            iFlag_geom_river_network = False

        if "iFlag_spac_river_network" in aConfig_in:
            iFlag_spac_river_network = aConfig_in["iFlag_spac_river_network"]
        else:
            iFlag_spac_river_network = False

        if "iFlag_geom_coastline" in aConfig_in:
            iFlag_geom_coastline = aConfig_in["iFlag_geom_coastline"]
        else:
            iFlag_geom_coastline = False

        if "iFlag_spac_coastline" in aConfig_in:
            iFlag_spac_coastline = aConfig_in["iFlag_spac_coastline"]
        else:
            iFlag_spac_coastline = False

        if "iFlag_geom_watershed_boundary" in aConfig_in:
            iFlag_geom_watershed_boundary = aConfig_in["iFlag_geom_watershed_boundary"]
        else:
            iFlag_geom_watershed_boundary = False

        if "iFlag_spac_watershed_boundary" in aConfig_in:
            iFlag_spac_watershed_boundary = aConfig_in["iFlag_spac_watershed_boundary"]
        else:
            iFlag_spac_watershed_boundary = False

        if "iFlag_geom_lake_boundary" in aConfig_in:
            iFlag_geom_lake_boundary = aConfig_in["iFlag_geom_lake_boundary"]
        else:
            iFlag_geom_lake_boundary = False

        if "iFlag_spac_lake_boundary" in aConfig_in:
            iFlag_spac_lake_boundary = aConfig_in["iFlag_spac_lake_boundary"]
        else:
            iFlag_spac_lake_boundary = False

        if "FULL_SPHERE_RADIUS" in aConfig_in:
            FULL_SPHERE_RADIUS = aConfig_in["FULL_SPHERE_RADIUS"]
        else:
            FULL_SPHERE_RADIUS = +6371.0

        if "geom_mshID" in aConfig_in:
            geom_mshID = aConfig_in["geom_mshID"]
        else:
            geom_mshID = "ellipsoid-mesh"

        if "spac_mshID" in aConfig_in:
            spac_mshID = aConfig_in["spac_mshID"]
        else:
            spac_mshID = "ellipsoid-grid" # "ellipsoid-grid"

        if "hfun_scal" in aConfig_in:
            hfun_scal = aConfig_in["hfun_scal"]
        else:
            hfun_scal = "absolute"

        if "hfun_hmax" in aConfig_in:
            hfun_hmax = aConfig_in["hfun_hmax"]
        else:
            hfun_hmax = float("inf")       # null spacing lim

        if "hfun_hmin" in aConfig_in:
            hfun_hmin = aConfig_in["hfun_hmin"]
        else:
            hfun_hmin = float(+0.00)

        if "mesh_dims" in aConfig_in:
            mesh_dims = aConfig_in["mesh_dims"]
        else:
            mesh_dims = +2                 # 2-dim. simplexes

        if "bisection" in aConfig_in:
            bisection = aConfig_in["bisection"]
        else:
            bisection = -1                 # call heutristic!

        if "optm_qlim" in aConfig_in:
            optm_qlim = aConfig_in["optm_qlim"]
        else:
            optm_qlim = +9.5E-01           # tighter opt. tol

        if "optm_iter" in aConfig_in:
            optm_iter = aConfig_in["optm_iter"]
        else:
            optm_iter = +32

        if "optm_qtol" in aConfig_in:
            optm_qtol = aConfig_in["optm_qtol"]
        else:
            optm_qtol = +1.0E-05

        if "dResolution_land" in aConfig_in:
            dResolution_land = float(aConfig_in["dResolution_land"])
        else:
            dResolution_land = 45.0

        if "dResolution_dam" in aConfig_in:
            dResolution_dam = float(aConfig_in["dResolution_dam"])
        else:
            dResolution_dam = 3.0

        if "dResolution_city" in aConfig_in:
            dResolution_city = float(aConfig_in["dResolution_city"])
        else:
            dResolution_city = 3.0


        if "dResolution_river_network" in aConfig_in:
            dResolution_river_network = float(aConfig_in["dResolution_river_network"])
        else:
            dResolution_river_network = 3.0

        if "dResolution_coastline" in aConfig_in:
            dResolution_coastline = float(aConfig_in["dResolution_coastline"])
        else:
            dResolution_coastline = 3.0

        if "dResolution_lake_boundary"  in aConfig_in:
            dResolution_lake_boundary = float(aConfig_in["dResolution_lake_boundary"])
        else:
            dResolution_lake_boundary = 3.0

        if "dResolution_watershed_boundary" in aConfig_in:
            dResolution_watershed_boundary = float(aConfig_in["dResolution_watershed_boundary"])
        else:
            dResolution_watershed_boundary = 3.0

        if "dhdx_lim" in aConfig_in:
            dhdx_lim = aConfig_in["dhdx_lim"]
        else:
            dhdx_lim = 0.25
    else:
        ncolumn_space = 360
        nrow_space = 180
        dSpac_value = 50.0
        iFlag_geom = False
        iFlag_spac = False
        iFlag_init = False
        iFlag_opts = False
        iFlag_ocean = False
        iFlag_land = False
        #point
        iFlag_geom_dam = False
        iFlag_geom_city = False
        #polyline
        iFlag_geom_river_network = False
        iFlag_geom_coastline = False
        #polygon
        iFlag_geom_watershed_boundary = False
        iFlag_geom_lake_boundary = False

        #spacing
        iFlag_spac_dam = False
        iFlag_spac_city = False
        iFlag_spac_river_network = False
        iFlag_spac_coastline = False
        iFlag_spac_watershed_boundary = False
        iFlag_spac_lake_boundary = False

        hfun_scal = "absolute"
        hfun_hmax = float("inf")       # null spacing lim
        hfun_hmin = float(+0.00)
        mesh_dims = +2                 # 2-dim. simplexes
        bisection = -1                 # call heutristic!
        optm_qlim = +9.5E-01           # tighter opt. tol
        optm_iter = +32
        optm_qtol = +1.0E-05
        dhdx_lim = 0.25                   # |dH/dx| thresh smaller values will make h(x) more smooth

    geom = jigsawpy.jigsaw_msh_t()
    if iFlag_geom:
        geom.mshID = geom_mshID
        geom.radii = np.full(3, FULL_SPHERE_RADIUS, dtype=geom.REALS_t)

        poly = jigsawpy.jigsaw_msh_t()
        nPolygon = 1
        print("BUILDING MESH GEOM.")

        #iFlag_dam = False
        #iFlag_city = False
        if iFlag_geom_dam :
            sFilename_dam_vector = aConfig_in["sFilename_dam_vector"]
            if not os.path.exists(sFilename_dam_vector):
                print("The file does not exist", sFilename_dam_vector)
            else:
                loadgeo(sFilename_dam_vector, poly)
                poly.point["coord"] *= np.pi / +180.
                addpoint(geom, poly, +1)
            pass

        if iFlag_geom_city:
            sFilename_city_vector = aConfig_in["sFilename_city_vector"]
            if not os.path.exists(sFilename_city_vector):
                print("The file does not exist", sFilename_city_vector)
            else:
                loadgeo(sFilename_city_vector, poly)
                poly.point["coord"] *= np.pi / +180.
                addpoint(geom, poly, +2)
            pass

        if iFlag_geom_river_network:
            sFilename_river_network_vector = aConfig_in["sFilename_river_network_vector"]
            if not os.path.exists(sFilename_river_network_vector):
                print("The file does not exist", sFilename_river_network_vector)
            else:
                loadgeo(sFilename_river_network_vector, poly)
                poly.point["coord"] *= np.pi / +180.
                addline(geom, poly, +1)
            pass

        if iFlag_geom_coastline:
            #skip the coastline for geom,
            #it will be only used for spac
            iFlag_coastline_geom = 1
            if iFlag_coastline_geom == 1:
                sFilename_coastline_vector = aConfig_in["sFilename_coastline_vector"]
                if not os.path.exists(sFilename_coastline_vector):
                    print("The file does not exist", sFilename_coastline_vector)
                else:
                    loadgeo(sFilename_coastline_vector, poly)
                    poly.point["coord"] *= np.pi / +180.
                    addpoly(geom, poly, nPolygon) #coastline is actually a polygon, not line
                    nPolygon += 1
                pass

        if iFlag_geom_watershed_boundary:
            iFlag_watershed_boundary_geom = 0
            if iFlag_watershed_boundary_geom == 1:
                sFilename_watershed_boundary_vector = aConfig_in["sFilename_watershed_boundary_vector"]
                if not os.path.exists(sFilename_watershed_boundary_vector):
                    print("The file does not exist", sFilename_watershed_boundary_vector)
                else:
                    loadgeo(sFilename_watershed_boundary_vector, poly)
                    poly.point["coord"] *= np.pi / +180.
                    addpoly(geom, poly, nPolygon)
                    nPolygon += 1
            pass

        if iFlag_geom_lake_boundary:
            sFilename_lake_boundary_vector = aConfig_in["sFilename_lake_boundary_vector"]
            if not os.path.exists(sFilename_lake_boundary_vector):
                print("The file does not exist", sFilename_lake_boundary_vector)
            else:
                loadgeo(sFilename_lake_boundary_vector, poly)
                poly.point["coord"] *= np.pi / +180.
                addpoly(geom, poly, nPolygon)
                nPolygon += 1
            pass

    spac = jigsawpy.jigsaw_msh_t()
    if iFlag_spac:
        print("Compute global h(x)... for spac")
        print(ncolumn_space, nrow_space, dResolution_land)
        spac.mshID = spac_mshID
        spac.radii = np.full( 3, FULL_SPHERE_RADIUS, dtype=spac.REALS_t)
        #for space setting base matrix, we prioritize the river network
        #if river network is not available, we will use the coastline
        if iFlag_spac_river_network:
            sFilename_river_network_raster = aConfig_in["sFilename_river_network_raster"]
            if not os.path.exists(sFilename_river_network_raster):
                pass
            else:
                dummy = gdal_read_geotiff_file(sFilename_river_network_raster)
            aData_river_network = dummy['dataOut']
            dOriginX = dummy['originX']
            dOriginY = dummy['originY']
            pixelHeight = dummy['pixelHeight'] #30/3600.0
            pixelWidth = dummy['pixelWidth']
            missingValue = dummy['missingValue']
            aRiver_network_mask = compute_mask(aData_river_network, 1, dOriginX, dOriginY, pixelWidth, pixelHeight, ncolumn_space, nrow_space)
            #we need to downsample the mask if the resolution is too high
            if ncolumn_space > 40000 and nrow_space > 20000:
                aRiver_network_mask = coarsen_mask(aRiver_network_mask, down=4)
                nrow_space = aRiver_network_mask.shape[0]
                ncolumn_space = aRiver_network_mask.shape[1]
        else:
            if iFlag_spac_coastline:
                sFilename_coastline_raster = aConfig_in["sFilename_coastline_raster"]
                #check if it exists
                if not os.path.exists(sFilename_coastline_raster):
                    print("The file does not exist")
                    print("The vector file will be converted to raster file")
                else:
                    #check spatial reference, reproject if needed
                    #sProjection =
                    dummy = gdal_read_geotiff_file(sFilename_coastline_raster) #1km, 30/3600.0
                aLand_ocean_mask = dummy['dataOut']
                dOriginX = dummy['originX']
                dOriginY = dummy['originY']
                pixelHeight = dummy['pixelHeight'] #30/3600.0
                pixelWidth = dummy['pixelWidth']
                missingValue = dummy['missingValue']
                aCoastline_mask = compute_mask(aLand_ocean_mask, 1, dOriginX, dOriginY, pixelWidth, pixelHeight, ncolumn_space, nrow_space)
                #we need to downsample the mask if the resolution is too high
                if ncolumn_space > 40000 and nrow_space > 20000:
                    aCoastline_mask = coarsen_mask(aCoastline_mask, down=4)
                    nrow_space = aCoastline_mask.shape[0]
                    ncolumn_space = aCoastline_mask.shape[1]
            else:
                print("You need to provide either coastline or river network to set the base spacing")

        #the low resolution grid
        #xgrid_low = np.linspace( -1. * np.pi, +1. * np.pi, 360)
        #ygrid_low = np.linspace( -.5 * np.pi, +.5 * np.pi, 180)

        #the high resolution grid
        xgrid_high = np.linspace( -1. * np.pi, +1. * np.pi, ncolumn_space)
        ygrid_high = np.linspace( -.5 * np.pi, +.5 * np.pi, nrow_space)

        spac.xgrid = xgrid_high
        spac.ygrid = ygrid_high
        spac.value = np.full( (nrow_space, ncolumn_space), dSpac_value, dtype=spac.REALS_t)

        if iFlag_ocean : #this requires the coastline
            if iFlag_spac_coastline:
                print("Compute global ocean h(x)...")
                vals = mdt.EC_CellWidthVsLat(spac.ygrid * 180. / np.pi)
                vals = np.reshape(vals, (spac.ygrid.size, 1))
                aOcean_value = np.array(np.tile( vals, (1, spac.xgrid.size)), dtype=spac.REALS_t)
                aOcean_mask = compute_mask(aLand_ocean_mask, missingValue, dOriginX, dOriginY, pixelWidth, pixelHeight, ncolumn_space, nrow_space)
                spac.value[aOcean_mask] =  np.minimum(aOcean_value[aOcean_mask], spac.value[aOcean_mask])
            else:
                print("You need to provide coastline to set the ocean spacing")
                pass
            pass

        if iFlag_spac_coastline:
            print("Compute global h(x)... for coastline")
            spac.value[aCoastline_mask] =  np.minimum(dResolution_coastline, spac.value[aCoastline_mask])
            pass

        if iFlag_land : #land feature, maybe mountain, etc
            print("Compute global h(x)... for land")
            aLand_mask = compute_mask(aLand_ocean_mask, 2, dOriginX, dOriginY, pixelWidth, pixelHeight, ncolumn_space, nrow_space)
            spac.value[aLand_mask] =  np.minimum(dResolution_land, spac.value[aLand_mask])
            pass

        if iFlag_spac_dam:
            print("Compute global h(x)... for dam")
            sFilename_dam_raster = aConfig_in["sFilename_dam_raster"]
            if not os.path.exists(sFilename_dam_raster):
                pass
            else:
                dummy = gdal_read_geotiff_file(sFilename_dam_raster)
            aData = dummy['dataOut']
            dOriginX = dummy['originX']
            dOriginY = dummy['originY']
            pixelWidth = dummy['pixelWidth']
            pixelHeight = dummy['pixelHeight']
            nDam = np.count_nonzero(aData == 1)
            print("The number of dams is", nDam)
            aDam_mask = compute_mask(aData, 1, dOriginX, dOriginY, pixelWidth, pixelHeight, ncolumn_space, nrow_space)
            spac.value[aDam_mask] =  np.minimum(dResolution_dam, spac.value[aDam_mask])

        if iFlag_spac_city:
            print("Compute global h(x)... for city")
            sFilename_city_raster = aConfig_in["sFilename_city_raster"]
            if not os.path.exists(sFilename_city_raster):
                pass
            else:
                dummy = gdal_read_geotiff_file(sFilename_city_raster)
            aData = dummy['dataOut']
            dOriginX = dummy['originX']
            dOriginY = dummy['originY']
            pixelWidth = dummy['pixelWidth']
            pixelHeight = dummy['pixelHeight']
            nCity = np.count_nonzero(aData == 1)
            print("The number of cities is", nCity)
            aCity_mask = compute_mask(aData, 1, dOriginX, dOriginY, pixelWidth, pixelHeight, ncolumn_space, nrow_space)
            spac.value[aCity_mask] =  np.minimum(dResolution_city, spac.value[aCity_mask])

        if iFlag_spac_river_network:
            print("Compute global h(x)... for river network")
            #sFilename_river_network_raster = aConfig_in["sFilename_river_network_raster"]
            #if not os.path.exists(sFilename_river_network_raster):
            #    pass
            #else:
            #    dummy = gdal_read_geotiff_file(sFilename_river_network_raster)
            #aData_river_network = dummy['dataOut']
            #dOriginX = dummy['originX']
            #dOriginY = dummy['originY']
            #pixelHeight = dummy['pixelHeight'] #30/3600.0
            #pixelWidth = dummy['pixelWidth']
            #aRiver_network_mask = compute_mask(aData_river_network, 1, dOriginX, dOriginY, pixelWidth, pixelHeight, ncolumn_space, nrow_space)
            spac.value[aRiver_network_mask] = np.minimum(dResolution_river_network, spac.value[aRiver_network_mask])
            pass

        if iFlag_spac_watershed_boundary:
            print("Compute global h(x)... for watershed boundary")
            sFilename_watershed_boundary_raster = aConfig_in["sFilename_watershed_boundary_raster"]
            if not os.path.exists(sFilename_watershed_boundary_raster):
                pass
            else:
                dummy = gdal_read_geotiff_file(sFilename_watershed_boundary_raster)
            aData = dummy['dataOut']
            dOriginX = dummy['originX']
            dOriginY = dummy['originY']
            pixelHeight = dummy['pixelHeight'] #30/3600.0
            pixelWidth = dummy['pixelWidth']
            aWatershed_mask = compute_mask(aData, 1, dOriginX, dOriginY, pixelWidth, pixelHeight, ncolumn_space, nrow_space)
            spac.value[aWatershed_mask] = np.minimum(dResolution_watershed_boundary, spac.value[aWatershed_mask])
            pass

        if iFlag_spac_lake_boundary:
            print("Compute global h(x)... for lake boundary")
            sFilename_lake_boundary_raster = aConfig_in["sFilename_lake_boundary_raster"]
            if not os.path.exists(sFilename_lake_boundary_raster):
                pass
            else:
                dummy = gdal_read_geotiff_file(sFilename_lake_boundary_raster)
            aData = dummy['dataOut']
            dOriginX = dummy['originX']
            dOriginY = dummy['originY']
            pixelHeight = dummy['pixelHeight'] #30/3600.0
            pixelWidth = dummy['pixelWidth']
            aLake_mask = compute_mask(aData, 1, dOriginX, dOriginY, pixelWidth, pixelHeight, ncolumn_space, nrow_space)
            spac.value[aLake_mask] =  np.minimum(dResolution_lake_boundary, spac.value[aLake_mask])
            pass

        spac.slope = np.full(spac.value.shape, dhdx_lim, dtype=spac.REALS_t)
        pass
    else:
        spac.mshID = spac_mshID
        spac.radii = np.full( 3, FULL_SPHERE_RADIUS, dtype=spac.REALS_t)
        xgrid_low = np.linspace( -1. * np.pi, +1. * np.pi, 360)
        ygrid_low = np.linspace( -.5 * np.pi, +.5 * np.pi, 180)
        spac.xgrid = xgrid_low
        spac.ygrid = ygrid_low
        spac.value = np.full( (nrow_space, ncolumn_space), dSpac_value, dtype=spac.REALS_t)
        pass

    opts = jigsawpy.jigsaw_jig_t()
    opts.verbosity=1
    if iFlag_opts:
        opts.hfun_scal = hfun_scal
        opts.hfun_hmax = hfun_hmax
        opts.hfun_hmin = hfun_hmin
        opts.mesh_dims = mesh_dims
        opts.bisection = bisection
        opts.optm_qlim = optm_qlim
        opts.optm_iter = optm_iter
        opts.optm_qtol = optm_qtol

    init = jigsawpy.jigsaw_msh_t()
    if iFlag_init:

        pass
#------------------------------------ setup files for JIGSAW

    sWorkspace_tmp = os.path.join(sWorkspace_jigsaw_in, "tmp")
    if os.path.exists(sWorkspace_tmp):
        #delete the folder first
        shutil.rmtree(sWorkspace_tmp)


    Path(sWorkspace_tmp).mkdir(parents=True, exist_ok=True)

    opts.geom_file = os.path.join(
        sWorkspace_jigsaw_in, "tmp", "geom.msh")

    opts.jcfg_file = os.path.join(
        sWorkspace_jigsaw_in, "tmp", "opts.jig")

    opts.init_file = os.path.join(
        sWorkspace_jigsaw_in, "tmp", "init.msh")

    opts.hfun_file = os.path.join(
        sWorkspace_jigsaw_in, "tmp", "spac.msh")

    opts.mesh_file = os.path.join(
        sWorkspace_jigsaw_in, "tmp", "mesh.msh")

    opts.hfun_tags = "precision = 9"    # less float prec.

    jigsawpy.savemsh(opts.geom_file, geom,
                     opts.geom_tags)

    jigsawpy.savemsh(opts.hfun_file, spac,
                     opts.hfun_tags)

    jigsawpy.savemsh(opts.init_file, init,
                     opts.init_tags)

    if iFlag_spac:
        # call jigsaw's Eikonal solver to impose the gradient constraint on h(x)
        jigsawpy.cmd.marche(opts, spac)
        pass

#------------------------------------ make mesh using JIGSAW

    if (not hasattr(opts, "bisection")):

        opts.bisection = +0

    if (opts.bisection < +0):           # bisect heuristic

        rbar = np.mean(geom.radii)
        hbar = np.mean(spac.value)

        nlev = round(math.log2(
            rbar / math.sin(.4 * math.pi) / hbar)
        )

        nlev = nlev - 1

        ttic = time.time()

        jigsawpy.cmd.tetris(opts, nlev - 0, mesh)

        ttoc = time.time()

        print("CPUSEC =", (ttoc - ttic))
        print("BISECT =", +nlev)

    elif (opts.bisection > +0):         # bisect specified

        nlev = opts.bisection

        ttic = time.time()

        jigsawpy.cmd.tetris(opts, nlev - 0, mesh)

        ttoc = time.time()

        print("CPUSEC =", (ttoc - ttic))
        print("BISECT =", +nlev)

    else:                               # do non-recursive

        ttic = time.time()

        jigsawpy.cmd.jigsaw(opts, mesh)

        ttoc = time.time()

        print("CPUSEC =", (ttoc - ttic))

#------------------------------------ form local projections

    gprj = copy.deepcopy(geom)          # local 2d objects
    mprj = copy.deepcopy(mesh)

    if (mesh.vert3.size > +0):
        project(geom, mesh, gprj, mprj, projector)

#------------------------------------ assign IDtag's to cell

    if (geom.bound is not None and
        geom.bound.size > +0):      # tags per polygon

        imin = np.amin(geom.bound["IDtag"])
        imax = np.amax(geom.bound["IDtag"])

        for itag in range(
                imin + 0, imax + 1):

            tagcell(geom, mesh, gprj, mprj, itag)

#------------------------------------ check mesh for quality

    cost = jigsawpy.triscr2(            # quality metrics!
        mesh.point["coord"],
        mesh.tria3["index"])

    print("TRISCR =", np.min(cost), np.mean(cost))

    cost = jigsawpy.pwrscr2(
        mesh.point["coord"],
        mesh.power,
        mesh.tria3["index"])

    print("PWRSCR =", np.min(cost), np.mean(cost))

    tbad = jigsawpy.centre2(
        mesh.point["coord"],
        mesh.power,
        mesh.tria3["index"])

    print("OBTUSE =",
          +np.count_nonzero(np.logical_not(tbad)))

    ndeg = jigsawpy.trideg2(
        mesh.point["coord"],
        mesh.tria3["index"])

    print("TOPOL. =",
          +np.count_nonzero(ndeg==+6) / ndeg.size)

#------------------------------------ save mesh for Paraview

    sWorkspace_tmp = os.path.join(
        sWorkspace_jigsaw_in, "out")
    Path(sWorkspace_tmp).mkdir(parents=True, exist_ok=True)

    jigsawpy.savevtk(os.path.join(
        sWorkspace_jigsaw_in, "out", "geom.vtk"), geom)
    jigsawpy.savevtk(os.path.join(
        sWorkspace_jigsaw_in, "out", "spac.vtk"), spac)
    jigsawpy.savevtk(os.path.join(
        sWorkspace_jigsaw_in, "out", "init.vtk"), init)
    jigsawpy.savevtk(os.path.join(
        sWorkspace_jigsaw_in, "out", "mesh.vtk"), mesh)
    jigsawpy.savevtk(os.path.join(
        sWorkspace_jigsaw_in, "out", "geom_prj.vtk"), gprj)
    jigsawpy.savevtk(os.path.join(
        sWorkspace_jigsaw_in, "out", "mesh_prj.vtk"), mprj)

    return geom, gprj, mesh, mprj


def project(geom, mesh, gprj, mprj, pmid):
    """
    PROJECT: projection of GEOM/MESH objects to a 2d. plane.

    Modifies GPRJ, MPRJ objects "inplace".

    """
    # Authors: Darren Engwirda

    proj = jigsawpy.jigsaw_prj_t()

    mprj.point = np.full(   mesh.point.size, +0, dtype=mprj.VERT2_t)

    mprj.point["coord"] =   jigsawpy.R3toS2(
        geom.radii, mesh.point["coord"])

    proj.prjID = "stereographic"
    proj.radii = np.mean(geom.radii)
    proj.xbase = pmid[0] * np.pi / +180.
    proj.ybase = pmid[1] * np.pi / +180.

    jigsawpy.project(gprj, proj, "fwd")
    jigsawpy.project(mprj, proj, "fwd")

    return


def tagcell(geom, mesh, gprj, mprj, itag):
    """
    TAGCELL: assign ID-tags to mesh cells based-on polygons
    defined in GEOM.BOUND.

    Modifies MESH, MPRJ objects "inplace".

    """
    # Authors: Darren Engwirda

    this = geom.bound["IDtag"] == itag
    cell = geom.bound["index"][this]

    loop = geom.edge2["index"][cell, :]

#------------------------------------ compute cell "centres"
    ball = jigsawpy.tribal2(
        mprj.point["coord"], mesh.tria3["index"])

#------------------------------------ calc. in-polygon tests
    tmsk, _ = inpoly2(
        ball[:, :-1], gprj.point["coord"], loop)

    vmsk, _ = inpoly2(mprj.point["coord"],
                      gprj.point["coord"], loop)

    vbnd = np.full(
        mesh.point.size, False, dtype=np.bool)

    vbnd[mesh.edge2["index"].flatten()] = True

    vmsk[vbnd] = False

    tbnd = np.logical_and.reduce((      # all vert. on bnd
        vbnd[mesh.tria3["index"][:, 0]],
        vbnd[mesh.tria3["index"][:, 1]],
        vbnd[mesh.tria3["index"][:, 2]]))

    tmsk[np.logical_not(tbnd)] = False

    mask = np.logical_or.reduce((       # vert. or cell in
        vmsk[mesh.tria3["index"][:, 0]],
        vmsk[mesh.tria3["index"][:, 1]],
        vmsk[mesh.tria3["index"][:, 2]],
        tmsk))

#------------------------------------ assign IDtags to cells
    mesh.tria3["IDtag"][mask == 1] = itag
    mprj.tria3["IDtag"][mask == 1] = itag

    return
