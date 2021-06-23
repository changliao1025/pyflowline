import os
from pystream.shared.vertex import pyvertex
from pystream.intersect.intersect_flowline_with_mesh import intersect_flowline_with_mesh

from pystream.simplification.remove_returning_flowline import remove_returning_flowline

sWorkspace_out = '/compyfs/liao313/04model/pyhexwatershed/columbia_river_basin'
sFilename_output_latlon = os.path.join(sWorkspace_out, 'lat_lon.json')
sFilename_output_square = os.path.join(sWorkspace_out, 'square.json')
sFilename_output_hexagon = os.path.join(sWorkspace_out, 'hexagon.json')

sFilename_out = 'flowline_merge2.json'
sFilename_flowline = os.path.join(sWorkspace_out, sFilename_out)

#lat-lon
sFilename_mesh=sFilename_output_latlon
sFilename_output= os.path.join(sWorkspace_out, 'lat_lon_intersect.json')
#intersect_flowline_with_mesh(sFilename_mesh, sFilename_flowline, sFilename_output)

sFilename_mesh=sFilename_output_square
sFilename_output= os.path.join(sWorkspace_out, 'square_intersect.json')
#intersect_flowline_with_mesh(sFilename_mesh, sFilename_flowline, sFilename_output)

sFilename_mesh=sFilename_output_hexagon
sFilename_output= os.path.join(sWorkspace_out, 'hexagon_intersect.json')
aHexagon, aHexagon_intersect = intersect_flowline_with_mesh(sFilename_mesh, sFilename_flowline, sFilename_output)

#simplify flowline
point= dict()
point['x'] = -2136506.345
point['y'] = 2901799.219
pVertex_outlet=pyvertex(point)

aHexagon_out = remove_returning_flowline(aHexagon_intersect, pVertex_outlet)


