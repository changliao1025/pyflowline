import os, sys
import numpy as np
import osgeo
from osgeo import ogr, osr, gdal, gdalconst
from pyearth.gis.gdal.gdal_function import obtain_raster_metadata
from pyearth.gis.gdal.gdal_function import reproject_coordinates
from pyearth.gis.projection.degree_to_meter import degree_to_meter

from pystream.mesh.hexagon.create_hexagon_mesh import create_hexagon_mesh
from pystream.mesh.square.create_latlon_mesh import create_latlon_mesh
from pystream.mesh.square.create_square_mesh import create_square_mesh
from pystream.mesh.jigsaw.create_mpas_mesh import create_mpas_mesh
from pystream.mesh.tin.create_tin_mesh import create_tin_mesh


def create_mesh(oModel_in):


    #we can use the dem extent to setup 
    dResolution = oModel_in.dResolution
    sFilename_dem = oModel_in.sFilename_dem
    sFilename_spatial_reference = oModel_in.sFilename_spatial_reference
    sFilename_mesh = oModel_in.sFilename_mesh

    sWorkspace_simulation_case = oModel_in.sWorkspace_simulation_case


    dPixelWidth, dOriginX, dOriginY, nrow, ncolumn, pSpatialRef, pProjection, pGeotransform = obtain_raster_metadata(sFilename_dem)

    spatial_reference_source = pSpatialRef
    spatial_reference_target = osr.SpatialReference()  
    spatial_reference_target.ImportFromEPSG(4326)

    dY_bot = dOriginY - (nrow+1) * dPixelWidth
    dLongitude_left,  dLatitude_bot= reproject_coordinates(dOriginX, dY_bot,pSpatialRef,spatial_reference_target)
    dX_right = dOriginX + (ncolumn +1) * dPixelWidth

    dLongitude_right, dLatitude_top= reproject_coordinates(dX_right, dOriginY,pSpatialRef,spatial_reference_target)
    dLatitude_mean = 0.5 * (dLatitude_top + dLatitude_bot)
    dResolution_meter = degree_to_meter(dLatitude_mean, dResolution )
    dX_left = dOriginX
    dY_top = dOriginY
    dArea = np.power(dResolution_meter,2.0)

    ncolumn= int( (dLongitude_right - dLongitude_left) / dResolution )
    nrow= int( (dLatitude_top - dLatitude_bot) / dResolution )
   
    ncolumn= int( (dX_right - dX_left) / dResolution_meter )
    nrow= int( (dY_top - dY_bot) / dResolution_meter )
   
    #hexagon edge
    dLength_edge = np.sqrt(  2.0 * dArea / (3.0* np.sqrt(3.0))  )
    dLength_shift = 0.5 * dLength_edge * np.sqrt(3.0)
    dX_spacing = dLength_edge * 1.5
    dY_spacing = dLength_edge * np.sqrt(3.0)

    ncolumn= int( (dX_right - dX_left) / dX_spacing )
    nrow= int( (dY_top - dY_bot) / dY_spacing )
   
    iMesh_type = oModel_in.iMesh_type
    if iMesh_type =='hexagon': #hexagon
        aHexagon = create_hexagon_mesh(dX_left, dY_bot, dResolution_meter, ncolumn, nrow, sFilename_mesh, sFilename_spatial_reference)
        return aHexagon
    else:
        if iMesh_type =='sqaure': #sqaure
            aSquare = create_square_mesh(dX_left, dY_bot, dResolution_meter, ncolumn, nrow, sFilename_mesh, sFilename_spatial_reference)
            return aSquare
        else:
            if iMesh_type =='latlon': #latlon
                aLatlon = create_latlon_mesh(dLongitude_left, dLatitude_bot, dResolution, ncolumn, nrow, sFilename_mesh)
                return aLatlon
            else:
                if iMesh_type =='mpas': #mpas
                    aMpas = create_mpas_mesh(oModel_in)
                    return aMpas
                else:
                    if iMesh_type =='tin': #tin
                        aTin = create_tin_mesh(oModel_in)
                        return aTin
                    else:
                        print('Unsupported mesh type?')
                        return
    
    

