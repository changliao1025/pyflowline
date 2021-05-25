import os, sys
import numpy as np
import osgeo
from osgeo import ogr, osr, gdal, gdalconst
from pystream.mesh.hexagon.create_hexagon_mesh import create_hexagon_mesh
from pystream.mesh.square.create_lat_lon_mesh import create_lat_lon_mesh
from pystream.mesh.square.create_square_mesh import create_square_mesh

from pyearth.gis.gdal.gdal_function import obtain_raster_metadata
from pyearth.gis.gdal.gdal_function import reproject_coordinates
from pyearth.gis.projection.degree_to_meter import degree_to_meter

dResolution=0.5

#we can use the dem extent to setup 
sFilename_geotiff = '/qfs/people/liao313/data/hexwatershed/columbia_river_basin/raster/dem/crbdem.tif'
dPixelWidth, dOriginX, dOriginY, nrow, ncolumn, pSpatialRef, pProjection, pGeotransform = obtain_raster_metadata(sFilename_geotiff)

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



sWorkspace_out = '/compyfs/liao313/04model/pyhexwatershed/columbia_river_basin'


sFilename_shapefile = '/qfs/people/liao313/data/hexwatershed/columbia_river_basin/vector/mesh_id/crb_flowline_remove_small_line_split.shp'


sFilename_output_latlon = os.path.join(sWorkspace_out, 'lat_lon.json')
sFilename_output_square = os.path.join(sWorkspace_out, 'square.json')
sFilename_output_hexagon = os.path.join(sWorkspace_out, 'hexagon.json')
#sFilename_output_mpas = os.path.join(sWorkspace_out, 'mpas.json')
#sFilename_output_tin = os.path.join(sWorkspace_out, 'tin.json')

ncolumn= int( (dLongitude_right - dLongitude_left) / dResolution )
nrow= int( (dLatitude_top - dLatitude_bot) / dResolution )
create_lat_lon_mesh(dLongitude_left, dLatitude_bot, dResolution, ncolumn, nrow, sFilename_output_latlon)

ncolumn= int( (dX_right - dX_left) / dResolution_meter )
nrow= int( (dY_top - dY_bot) / dResolution_meter )
create_square_mesh(dX_left, dY_bot, dResolution_meter, ncolumn, nrow, sFilename_output_square, sFilename_shapefile)
#hexagon edge
dLength_edge = np.sqrt(  2.0 * dArea / (3.0* np.sqrt(3.0))  )
dLength_shift = 0.5 * dLength_edge * np.sqrt(3.0)
dX_spacing = dLength_edge * 1.5
dY_spacing = dLength_edge * np.sqrt(3.0)

ncolumn= int( (dX_right - dX_left) / dX_spacing )
nrow= int( (dY_top - dY_bot) / dY_spacing )
create_hexagon_mesh(dX_left, dY_bot, dResolution_meter, ncolumn, nrow, sFilename_output_hexagon, sFilename_shapefile)