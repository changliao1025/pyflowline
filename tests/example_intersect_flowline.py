import os
from pystream.shared.vertex import pyvertex
from pystream.intersect.intersect_flowline_with_mesh import intersect_flowline_with_mesh

from pystream.simplification.remove_returning_flowline import remove_returning_flowline

from pystream.format.read_flowline_geojson import read_flowline_geojson
from pystream.format.export_flowline_to_json import export_flowline_to_json
from pystream.correct_flowline_direction import correct_flowline_direction
from pystream.loop.remove_flowline_loop import remove_flowline_loop
from pystream.split.find_flowline_vertex import find_flowline_vertex
from pystream.split.split_flowline import split_flowline
from pystream.format.export_vertex_to_json import export_vertex_to_json

sWorkspace_out = '/compyfs/liao313/04model/pyhexwatershed/columbia_river_basin'
sFilename_output_latlon = os.path.join(sWorkspace_out, 'lat_lon.json')
sFilename_output_square = os.path.join(sWorkspace_out, 'square.json')
sFilename_output_hexagon = os.path.join(sWorkspace_out, 'hexagon.json')

sFilename_out = 'flowline_merge2.json'
sFilename_out = 'flowline_segment_order.json'
sFilename_flowline = os.path.join(sWorkspace_out, sFilename_out)
aFlowline0, pSpatialRef = read_flowline_geojson(sFilename_flowline)

sFilename_mesh=sFilename_output_hexagon
sFilename_output= os.path.join(sWorkspace_out, 'hexagon_intersect.json')
aCell, aCell_intersect, aFlowline_intersect = intersect_flowline_with_mesh(1, sFilename_mesh, sFilename_flowline, sFilename_output)

sFilename_mesh=sFilename_output_square
sFilename_output= os.path.join(sWorkspace_out, 'square_intersect.json')
#aCell, aCell_intersect, aFlowline_intersect = intersect_flowline_with_mesh(2, sFilename_mesh, sFilename_flowline, sFilename_output)

#lat-lon
sFilename_mesh=sFilename_output_latlon
sFilename_output= os.path.join(sWorkspace_out, 'lat_lon_intersect.json')
#aCell, aCell_intersect, aFlowline_intersect = intersect_flowline_with_mesh(3, sFilename_mesh, sFilename_flowline, sFilename_output)


#simplify flowline


point= dict()
point['x'] = -2136506.345
point['y'] = 2901799.219
pVertex_outlet=pyvertex(point)

aCell, aFlowline, aFlowline_no_parallel = remove_returning_flowline(aCell, aCell_intersect, pVertex_outlet)
sFilename_out = 'flowline_simplified.json'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
export_flowline_to_json( aFlowline, pSpatialRef, sFilename_out)
#sFilename_out = 'flowline_simplified_no_parallel.json'
#sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
#export_flowline_to_json( aFlowline_no_parallel, pSpatialRef, sFilename_out)

pVertex_outlet=aFlowline[0].pVertex_end
aVertex = find_flowline_vertex(aFlowline)
sFilename_out = 'flowline_vertex_without_confluence2.json'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
export_vertex_to_json( aVertex,pSpatialRef, sFilename_out)

aFlowline = split_flowline(aFlowline, aVertex)
sFilename_out = 'flowline_split_by_point2.json'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
export_flowline_to_json( aFlowline,pSpatialRef, sFilename_out)

aFlowline= correct_flowline_direction(aFlowline,  pVertex_outlet )

aFlowline = remove_flowline_loop(  aFlowline)    
sFilename_out = 'flowline_loop_no_parallel.json'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
export_flowline_to_json( aFlowline,pSpatialRef, sFilename_out)


