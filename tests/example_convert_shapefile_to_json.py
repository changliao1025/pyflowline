import os
from pystream.format.convert_shapefile_to_json import convert_shapefile_to_json
sFilename_shapefile_in = '/qfs/people/liao313/data/hexwatershed/columbia_river_basin/vector/mesh_id/crb_flowline_remove_small_line_split.shp'
sFilename_json_out = 'flowline.json'
sWorkspace_out = '/compyfs/liao313/04model/pyhexwatershed/columbia_river_basin'
sFilename_json_out = os.path.join(sWorkspace_out, sFilename_json_out)

convert_shapefile_to_json(sFilename_shapefile_in, sFilename_json_out)