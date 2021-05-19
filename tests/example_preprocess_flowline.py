import os, sys

import numpy as np
from pystream.shared.vertex import pyvertex
from pyearth.system.define_global_variables import *

from pystream.format.read_flowline_shapefile import read_flowline_shapefile

from pystream.format.export_flowline_to_json import export_flowline_to_json
from pystream.format.export_vertex_to_json import export_vertex_to_json

from pystream.connect.connect_disconnect_flowline import connect_disconnect_flowline
from pystream.correct_flowline_direction import correct_flowline_direction
#
#from hexwatershed.preprocess.stream.merge.merge_flowline import merge_flowline
#from hexwatershed.preprocess.stream.merge.merge_flowline2 import merge_flowline2
#
from pystream.split.split_flowline import split_flowline
#from hexwatershed.preprocess.stream.split.split_flowline2 import split_flowline2
#from hexwatershed.preprocess.stream.split.split_flowline3 import split_flowline3
#
from pystream.split.find_flowline_vertex import find_flowline_vertex
#from hexwatershed.preprocess.stream.split.find_flowline_vertex2 import find_flowline_vertex2
#
#from hexwatershed.preprocess.stream.simplification.remove_flowline_loop import remove_flowline_loop
#from hexwatershed.preprocess.stream.simplification.remove_flowline_loop2 import remove_flowline_loop2
#from hexwatershed.preprocess.stream.simplification.remove_flowline_loop3 import remove_flowline_loop3
#
#from hexwatershed.preprocess.stream.simplification.remove_small_river import remove_small_river
#
#from hexwatershed.preprocess.stream.define_stream_order import define_stream_order
#from hexwatershed.preprocess.mesh.intersect_flowline_with_mesh import intersect_flowline_with_mesh


"""
prepare the flowline using multiple step approach
"""
sFilename_shapefile_in = '/qfs/people/liao313/data/hexwatershed/columbia_river_basin/vector/hydrology/crb_flowline.shp'

sFilename_mesh = 'hexagon.json'

sWorkspace_out = '/compyfs/liao313/04model/pyhexwatershed/columbia_river_basin'

sFilename_mesh_in = os.path.join(sWorkspace_out, sFilename_mesh)

#read shapefile and store information in the list
aFlowline, pSpatialRef = read_flowline_shapefile(sFilename_shapefile_in)
#we also need to save the spatial reference information for the output purpose

sFilename_out = 'flowline.json'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
export_flowline_to_json( aFlowline,pSpatialRef, sFilename_out)

aVertex=list()
point= dict()
point['x'] = -1589612.188
point['y'] = 3068975.112
pVertex=pyvertex(point)
aVertex.append(pVertex)
point['x'] =  -1568732.491
point['y'] = 3064177.639
pVertex=pyvertex(point)
aVertex.append(pVertex)


aThreshold = np.full(2, 300.0, dtype=float)
aFlowline = connect_disconnect_flowline(aFlowline, aVertex, aThreshold)
sFilename_out = 'flowline_connect.json'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
export_flowline_to_json( aFlowline,pSpatialRef, sFilename_out)


aVertex = find_flowline_vertex(aFlowline)
sFilename_out = 'flowline_vertex_without_confluence.json'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
export_vertex_to_json( aVertex,pSpatialRef, sFilename_out)

aFlowline = split_flowline(aFlowline, aVertex)
sFilename_out = 'flowline_split_by_point.json'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
export_flowline_to_json( aFlowline,pSpatialRef, sFilename_out)

#ues location to find outlet
point['x'] = -2136506.345
point['y'] = 2901799.219
pVertex=pyvertex(point)
aFlowline= correct_flowline_direction(aFlowline,  pVertex )
sFilename_in = sFilename_out
sFilename_out = sWorkspace_out + slash + 'flowline_direction.json'
export_flowline_to_json( aFlowline,pSpatialRef, sFilename_out)

#step 4: remove loops

sFilename_in = sFilename_out    
sFilename_out = sWorkspace_out + slash + 'flowline_loop.json'
remove_flowline_loop3(sFilename_in,  sFilename_out)    
sFilename_in = sFilename_out    
sFilename_out = sWorkspace_out + slash + 'flowline_large.json'
#sFilename_in = sWorkspace_out + slash + 'flowline_large.shp'
#remove_small_river(sFilename_in, sFilename_out, 1.0E4)



sFilename_in = sFilename_out
sFilename_out = sWorkspace_out + slash + 'flowline_vertex_with_confluence.json'
#find_flowline_vertex2(sFilename_in,  sFilename_out)

sFilename_in = sWorkspace_out + slash + 'flowline_large.json'

sFilename_in2 = sFilename_out
sFilename_out = sWorkspace_out + slash + 'flowline_merge.json'
#merge_flowline2(sFilename_in, sFilename_in2,  sFilename_out )    
sFilename_in = sFilename_out    
sFilename_out = sWorkspace_out + slash + 'flowline_large2.json'
#remove_small_river(sFilename_in, sFilename_out, 1.0E4)
sFilename_in = sFilename_out
sFilename_out = sWorkspace_out + slash + 'flowline_vertex_with_confluence2.json'
#find_flowline_vertex2(sFilename_in,  sFilename_out)

sFilename_in = sWorkspace_out + slash + 'flowline_large2.json'

sFilename_in2 = sFilename_out
sFilename_out = sWorkspace_out + slash + 'flowline_merge2.json'
#merge_flowline2(sFilename_in, sFilename_in2,  sFilename_out )    

#step 2: merge all as one single feature     
#build stream order 
sFilename_in = sFilename_out
sFilename_out = sWorkspace_out + slash + 'flowline_order.json'
define_stream_order(sFilename_in, sFilename_out)

#step 5: remove small headwater segment

#step 6: intersect with mesh and simplify
sFilename_mesh = '/compyfs/liao313/04model/pyhexwatershed/columbia_river_basin/hexagon.json'
sFilename_flowline= sFilename_out
sFilename_out = sWorkspace_out + slash + 'flowline_intersect.json'
#intersect_flowline_with_mesh(sFilename_mesh, sFilename_flowline, sFilename_out)
#step 7: rebuild index and order
#step 8: calculate properties
print('Finished')




    
