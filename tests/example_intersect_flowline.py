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
from pystream.pystream_read_model_configuration_file import pystream_read_model_configuration_file
from pystream.case import streamcase


sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pystream/pystream/case.xml'

    #step 1
aParameter = pystream_read_model_configuration_file(sFilename_configuration_in)

       # iCase_index_in=iCase_index_in, sJob_in=sJob_in, iFlag_mode_in=iFlag_mode_in)

aParameter['sFilename_model_configuration'] = sFilename_configuration_in
oModel = streamcase(aParameter)

#sWorkspace_simulation_case = '/compyfs/liao313/04model/pyhexwatershed/columbia_river_basin'
#sWorkspace_simulation_case = '/people/liao313/tmp/susquehanna'
#sFilename_output_latlon = os.path.join(sWorkspace_simulation_case, 'lat_lon.json')
#sFilename_output_square = os.path.join(sWorkspace_simulation_case, 'square.json')
#sFilename_output_hexagon = os.path.join(sWorkspace_simulation_case, 'hexagon.json')
#sFilename_output_mpas = os.path.join(sWorkspace_simulation_case, 'mpas.json')




#sFilename_out = 'flowline_merge2.json'
sFilename_out = 'flowline_segment_order.shp'
sWorkspace_simulation_case = oModel.sWorkspace_simulation_case
#sFilename_out = 



sFilename_flowline = os.path.join(sWorkspace_simulation_case, sFilename_out)
aFlowline0, pSpatialRef = read_flowline_geojson(sFilename_flowline)

sFilename_mesh=oModel.sFilename_mesh
sFilename_intersect = oModel.sFilename_intersect
#sFilename_output= os.path.join(sWorkspace_simulation_case, 'hexagon_intersect.json')
#aCell, aCell_intersect, aFlowline_intersect = intersect_flowline_with_mesh(1, sFilename_mesh, sFilename_flowline, sFilename_output)

#sFilename_mesh=sFilename_output_square
#sFilename_output= os.path.join(sWorkspace_simulation_case, 'square_intersect.json')
#aCell, aCell_intersect, aFlowline_intersect = intersect_flowline_with_mesh(2, sFilename_mesh, sFilename_flowline, sFilename_output)

#lat-lon
#sFilename_mesh=sFilename_output_latlon
#sFilename_output= os.path.join(sWorkspace_simulation_case, 'lat_lon_intersect.json')
#aCell, aCell_intersect, aFlowline_intersect = intersect_flowline_with_mesh(3, sFilename_mesh, sFilename_flowline, sFilename_output)

#mpas
#sFilename_output= os.path.join(sWorkspace_simulation_case, 'mpas_intersect.shp')
aCell, aCell_intersect, aFlowline_intersect = intersect_flowline_with_mesh(4, sFilename_mesh, sFilename_flowline, sFilename_intersect)


#simplify flowline


point= dict()
point['x'] = -2136506.345
point['y'] = 2901799.219
#susquehanna -76.08247,39.55600,1683240.965,2014539.724
point['x'] = 1683240.965
point['y'] = 2014539.724


pVertex_outlet=pyvertex(point)

aCell, aFlowline, aFlowline_no_parallel = remove_returning_flowline(aCell, aCell_intersect, pVertex_outlet)
sFilename_out = 'flowline_simplified.shp'
sFilename_out = os.path.join(sWorkspace_simulation_case, sFilename_out)

sFilename_flowline_simplified = oModel.sFilename_flowline_simplified
export_flowline_to_json( aFlowline, pSpatialRef, sFilename_flowline_simplified)


pVertex_outlet=aFlowline[0].pVertex_end
aVertex = find_flowline_vertex(aFlowline)
#sFilename_out = 'flowline_vertex_without_confluence2.shp'
#sFilename_out = os.path.join(sWorkspace_simulation_case, sFilename_out)
sFilename_vertex_without_confluence_after_intersection  = oModel.sFilename_vertex_without_confluence_after_intersection
export_vertex_to_json( aVertex,pSpatialRef, sFilename_vertex_without_confluence_after_intersection)

aFlowline = split_flowline(aFlowline, aVertex)
sFilename_out = 'flowline_split_by_point2.shp'
sFilename_out = os.path.join(sWorkspace_simulation_case, sFilename_out)
export_flowline_to_json( aFlowline,pSpatialRef, sFilename_out)

aFlowline= correct_flowline_direction(aFlowline,  pVertex_outlet )

aFlowline = remove_flowline_loop(  aFlowline)    
sFilename_out = 'flowline_loop_no_parallel.shp'
sFilename_out = os.path.join(sWorkspace_simulation_case, sFilename_out)
export_flowline_to_json( aFlowline,pSpatialRef, sFilename_out)


