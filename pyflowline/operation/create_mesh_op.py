import os, sys
import numpy as np
import osgeo
from osgeo import ogr, osr, gdal, gdalconst
from pyflowline.algorithm.auxiliary.gdal_function import obtain_raster_metadata
from pyflowline.algorithm.auxiliary.reproject_coordinates import reproject_coordinates
from pyflowline.algorithm.auxiliary.gdal_function  import degree_to_meter
from pyflowline.algorithm.auxiliary.gdal_function  import meter_to_degree
from pyflowline.mesh.hexagon.create_hexagon_mesh import create_hexagon_mesh
from pyflowline.mesh.latlon.create_latlon_mesh import create_latlon_mesh
from pyflowline.mesh.square.create_square_mesh import create_square_mesh
from pyflowline.mesh.mpas.create_mpas_mesh import create_mpas_mesh
from pyflowline.mesh.tin.create_tin_mesh import create_tin_mesh


def create_mesh_op(oPyflowline_in):


    #we can use the dem extent to setup 
    iFlag_global =  oPyflowline_in.iFlag_global
    iMesh_type = oPyflowline_in.iMesh_type
    iFlag_save_mesh = oPyflowline_in.iFlag_save_mesh
    iFlag_rotation = oPyflowline_in.iFlag_rotation
    
    dResolution = oPyflowline_in.dResolution
    dResolution_meter = oPyflowline_in.dResolution_meter
    

    sFilename_dem = oPyflowline_in.sFilename_dem
    sFilename_spatial_reference = oPyflowline_in.sFilename_spatial_reference
    sFilename_mesh = oPyflowline_in.sFilename_mesh

    if iMesh_type !=4: #hexagon

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
    else:
        pass
   
    
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
                dLatitude_top    = oPyflowline_in.dLatitude_top   
                dLatitude_bot    = oPyflowline_in.dLatitude_bot   
                dLongitude_left  = oPyflowline_in.dLongitude_left 
                dLongitude_right = oPyflowline_in.dLongitude_right
                ncolumn= int( (dLongitude_right - dLongitude_left) / dResolution )
                nrow= int( (dLatitude_top - dLatitude_bot) / dResolution )
                aLatlon = create_latlon_mesh(dLongitude_left, dLatitude_bot, dResolution, ncolumn, nrow, \
                    sFilename_mesh)
                return aLatlon
            else:
                if iMesh_type == 4: #mpas
                    iFlag_use_mesh_dem = oPyflowline_in.iFlag_use_mesh_dem
                    sFilename_mesh_netcdf = oPyflowline_in.sFilename_mesh_netcdf
                    dLatitude_top    = oPyflowline_in.dLatitude_top   
                    dLatitude_bot    = oPyflowline_in.dLatitude_bot   
                    dLongitude_left  = oPyflowline_in.dLongitude_left 
                    dLongitude_right = oPyflowline_in.dLongitude_right
                    aMpas = create_mpas_mesh(iFlag_global, iFlag_use_mesh_dem, iFlag_save_mesh, \
                            dLatitude_top, dLatitude_bot, dLongitude_left, dLongitude_right,\
                                sFilename_mesh_netcdf,      sFilename_mesh)
                    return aMpas
                else:
                    if iMesh_type ==5: #tin this one need to be updated because central location issue
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
    
    

