import os, sys
import numpy as np
import osgeo
from osgeo import ogr, osr, gdal, gdalconst
from pystream.mesh.jigsaw.create_mpas_mesh import create_mpas_mesh


sFilename_mesh_in = '/qfs/people/liao313/data/hexcoastal/delaware/vector/mesh/invert_mesh.nc'


sWorkspace_out = '/compyfs/liao313/04model/pyhexwatershed/columbia_river_basin'
sWorkspace_out = '/compyfs/liao313/04model/pyhexwatershed/susquehanna'
sWorkspace_out = '/people/liao313/tmp/susquehanna'

dLongitude_left = -79.44374
dLatitude_bot = 39.00 #,1399152.687,1978258.386

dLongitude_right = -74.24774 
dLatitude_top = 43.00334 # ,1748363.409,2424316.881

sFilename_output_mpas = os.path.join(sWorkspace_out, 'mpas.json')

sFilename_shapefile = '/qfs/people/liao313/data/hexwatershed/susquehanna/vector/hydrology/stream_order78.shp'

aMpas = create_mpas_mesh(sFilename_mesh_in, dLatitude_top, dLatitude_bot, dLongitude_left, dLongitude_right,sFilename_output_mpas)

print('finished')