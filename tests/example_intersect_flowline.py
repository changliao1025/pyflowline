import os
from pystream.intersect.intersect_flowline_with_mesh import intersect_flowline_with_mesh

sWorkspace_out = '/compyfs/liao313/04model/pyhexwatershed/columbia_river_basin'
sFilename_output_latlon = os.path.join(sWorkspace_out, 'lat_lon.json')
sFilename_output_square = os.path.join(sWorkspace_out, 'square.json')
sFilename_output_hexagon = os.path.join(sWorkspace_out, 'hexagon.json')

sFilename_out = 'flowline_merge2.json'
sFilename_flowline = os.path.join(sWorkspace_out, sFilename_out)

#lat-lon
sFilename_mesh=sFilename_output_latlon
sFilename_output= os.path.join(sWorkspace_out, 'lat_lon_intersect.json')
intersect_flowline_with_mesh(sFilename_mesh, sFilename_flowline, sFilename_output)

sFilename_mesh=sFilename_output_square
sFilename_output= os.path.join(sWorkspace_out, 'square_intersect.json')
intersect_flowline_with_mesh(sFilename_mesh, sFilename_flowline, sFilename_output)

sFilename_mesh=sFilename_output_hexagon
sFilename_output= os.path.join(sWorkspace_out, 'hexagon_intersect.json')
intersect_flowline_with_mesh(sFilename_mesh, sFilename_flowline, sFilename_output)