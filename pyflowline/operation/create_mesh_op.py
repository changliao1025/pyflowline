import os, sys
import numpy as np
import osgeo
from osgeo import ogr, osr, gdal, gdalconst
from pyearth.gis.gdal.gdal_function import obtain_raster_metadata
from pyearth.gis.gdal.gdal_function import reproject_coordinates
from pyearth.gis.projection.degree_to_meter import degree_to_meter
from pyearth.gis.projection.meter_to_degree import meter_to_degree
from pyflowline.mesh.hexagon.create_hexagon_mesh import create_hexagon_mesh
from pyflowline.mesh.latlon.create_latlon_mesh import create_latlon_mesh
from pyflowline.mesh.square.create_square_mesh import create_square_mesh
from pyflowline.mesh.mpas.create_mpas_mesh import create_mpas_mesh
from pyflowline.mesh.tin.create_tin_mesh import create_tin_mesh


def create_mesh_op(opyflowline_in):


    #we can use the dem extent to setup 
    iMesh_type = opyflowline_in.iMesh_type
    iFlag_rotation = opyflowline_in.iFlag_rotation
    
    dResolution = opyflowline_in.dResolution
    dResolution_meter = opyflowline_in.dResolution_meter
    

    sFilename_dem = opyflowline_in.sFilename_dem
    sFilename_spatial_reference = opyflowline_in.sFilename_spatial_reference
    sFilename_mesh = opyflowline_in.sFilename_mesh

    dPixelWidth, dOriginX, dOriginY, nrow, ncolumn, pSpatialRef_dem, pProjection, pGeotransform\
         = obtain_raster_metadata(sFilename_dem)

    spatial_reference_source = pSpatialRef_dem
    spatial_reference_target = osr.SpatialReference()  
    spatial_reference_target.ImportFromEPSG(4326)

    dY_bot = dOriginY - (nrow+1) * dPixelWidth
    dLongitude_left,  dLatitude_bot= reproject_coordinates(dOriginX, dY_bot,pSpatialRef_dem,spatial_reference_target)
    dX_right = dOriginX + (ncolumn +1) * dPixelWidth

    dLongitude_right, dLatitude_top= reproject_coordinates(dX_right, dOriginY,pSpatialRef_dem,spatial_reference_target)
    dLatitude_mean = 0.5 * (dLatitude_top + dLatitude_bot)


    if dResolution_meter < 0:
        #not used
        pass
    else:
        dResolution = meter_to_degree(dResolution_meter, dLatitude_mean)
        

    dX_left = dOriginX
    dY_top = dOriginY
   
    
    if iMesh_type ==1: #hexagon

        #hexagon edge
        dResolution_meter = degree_to_meter(dLatitude_mean, dResolution )
        dArea = np.power(dResolution_meter,2.0)
        dLength_edge = np.sqrt(  2.0 * dArea / (3.0* np.sqrt(3.0))  )
        if iFlag_rotation ==0:            
            dX_spacing = dLength_edge * np.sqrt(3.0)
            dY_spacing = dLength_edge * 1.5
            ncolumn= int( (dX_right - dX_left) / dX_spacing )
            nrow= int( (dY_top - dY_bot) / dY_spacing ) 
        else:            
            dX_spacing = dLength_edge * 1.5
            dY_spacing = dLength_edge * np.sqrt(3.0)    
            ncolumn= int( (dX_right - dX_left) / dX_spacing )+1
            nrow= int( (dY_top - dY_bot) / dY_spacing )

        aHexagon = create_hexagon_mesh(iFlag_rotation, dX_left, dY_bot, dResolution_meter, ncolumn, nrow, \
            sFilename_mesh, sFilename_spatial_reference)
        return aHexagon
    else:
        if iMesh_type ==2: #sqaure
            ncolumn= int( (dX_right - dX_left) / dResolution_meter )
            nrow= int( (dY_top - dY_bot) / dResolution_meter )
            
            aSquare = create_square_mesh(dX_left, dY_bot, dResolution_meter, ncolumn, nrow, \
                sFilename_mesh, sFilename_spatial_reference)
            return aSquare
        else:
            if iMesh_type ==3: #latlon
                dResolution_meter = degree_to_meter(dLatitude_mean, dResolution)
                dArea = np.power(dResolution_meter,2.0)
                dLatitude_top    = opyflowline_in.dLatitude_top   
                dLatitude_bot    = opyflowline_in.dLatitude_bot   
                dLongitude_left  = opyflowline_in.dLongitude_left 
                dLongitude_right = opyflowline_in.dLongitude_right
                ncolumn= int( (dLongitude_right - dLongitude_left) / dResolution )
                nrow= int( (dLatitude_top - dLatitude_bot) / dResolution )
                aLatlon = create_latlon_mesh(dLongitude_left, dLatitude_bot, dResolution, ncolumn, nrow, \
                    sFilename_mesh)
                return aLatlon
            else:
                if iMesh_type ==4: #mpas
                    iFlag_use_mpas_dem = opyflowline_in.iFlag_use_mpas_dem
                    sFilename_mesh_netcdf = opyflowline_in.sFilename_mesh_netcdf
                    dLatitude_top    = opyflowline_in.dLatitude_top   
                    dLatitude_bot    = opyflowline_in.dLatitude_bot   
                    dLongitude_left  = opyflowline_in.dLongitude_left 
                    dLongitude_right = opyflowline_in.dLongitude_right
                    aMpas = create_mpas_mesh(iFlag_use_mpas_dem, \
                        sFilename_mesh_netcdf, \
                            dLatitude_top, dLatitude_bot, dLongitude_left, dLongitude_right,\
                        sFilename_mesh)
                    return aMpas
                else:
                    if iMesh_type ==5: #tin
                        #tin edge
                        dArea = np.power(dResolution_meter,2.0)
                        dLength_edge = np.sqrt(  4.0 * dArea /  np.sqrt(3.0) )  
                        dX_shift = 0.5 * dLength_edge
                        dY_shift = 0.5 * dLength_edge * np.sqrt(3.0) 
                        dX_spacing = dX_shift * 2
                        dY_spacing = dY_shift
                        ncolumn= int( (dX_right - dX_left) / dX_shift )
                        nrow= int( (dY_top - dY_bot) / dY_spacing ) 
                        aTin = create_tin_mesh(dX_left, dY_bot, dResolution_meter, ncolumn, nrow,sFilename_mesh, sFilename_spatial_reference)
                        return aTin
                    else:
                        print('Unsupported mesh type?')
                        return
    
    

