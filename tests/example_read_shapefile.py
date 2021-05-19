import os, sys
from pyearth.system.define_global_variables import *
from pystream.format.read_flowline_shapefile import read_flowline_shapefile

from pystream.format.export_to_json import export_to_json

sFilename_shapefile_in = '/qfs/people/liao313/data/hexwatershed/columbia_river_basin/vector/hydrology/crb_flowline.shp'
sFilename_mesh = 'hexagon.json'
sWorkspace_out = '/compyfs/liao313/04model/pyhexwatershed/columbia_river_basin'
sFilename_mesh_in = os.path.join(sWorkspace_out, sFilename_mesh)
sFilename_json_out = 'flowline.json'
sWorkspace_out = '/compyfs/liao313/04model/pyhexwatershed/columbia_river_basin'
sFilename_json_out = os.path.join(sWorkspace_out, sFilename_json_out)


aFlowline, pSpatialRef =  read_flowline_shapefile(sFilename_shapefile_in)
sFilename_json_out = sWorkspace_out + slash + 'flowline.json'
export_to_json(aFlowline, pSpatialRef, sFilename_json_out)