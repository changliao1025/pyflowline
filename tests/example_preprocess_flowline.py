import os, sys
from pystream.define_stream_segment_index import define_stream_segment_index
from pystream.split.find_flowline_confluence import find_flowline_confluence

import numpy as np
from pystream.shared.vertex import pyvertex
from pyearth.system.define_global_variables import *

from pystream.format.read_flowline_shapefile import read_flowline_shapefile

from pystream.format.export_flowline_to_json import export_flowline_to_json
from pystream.format.export_vertex_to_json import export_vertex_to_json

from pystream.connect.connect_disconnect_flowline import connect_disconnect_flowline
from pystream.correct_flowline_direction import correct_flowline_direction
#
from pystream.merge.merge_flowline import merge_flowline

#
from pystream.split.split_flowline import split_flowline
#from hexwatershed.preprocess.stream.split.split_flowline2 import split_flowline2
#from hexwatershed.preprocess.stream.split.split_flowline3 import split_flowline3
#
from pystream.split.find_flowline_vertex import find_flowline_vertex
from pystream.split.find_flowline_confluence import find_flowline_confluence

from pystream.loop.remove_flowline_loop import remove_flowline_loop
#
from pystream.simplification.remove_small_river import remove_small_river
from pystream.define_stream_order import define_stream_order



"""
prepare the flowline using multiple step approach
"""
sFilename_shapefile_in = '/qfs/people/liao313/data/hexwatershed/columbia_river_basin/vector/hydrology/crb_flowline.shp'

sFilename_mesh = 'hexagon.json'

sWorkspace_out = '/compyfs/liao313/04model/pyhexwatershed/columbia_river_basin'
#sWorkspace_out = '/qfs/people/liao313/tmp/columbia_river_basin'


sFilename_mesh_in = os.path.join(sWorkspace_out, sFilename_mesh)

#read shapefile and store information in the list
aFlowline, pSpatialRef = read_flowline_shapefile(sFilename_shapefile_in)
#we also need to save the spatial reference information for the output purpose

sFilename_out = 'flowline.json'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
#export_flowline_to_json( aFlowline,pSpatialRef, sFilename_out)

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
#export_flowline_to_json( aFlowline,pSpatialRef, sFilename_out)


aVertex = find_flowline_vertex(aFlowline)
sFilename_out = 'flowline_vertex_without_confluence.json'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
#export_vertex_to_json( aVertex,pSpatialRef, sFilename_out)

aFlowline = split_flowline(aFlowline, aVertex)
sFilename_out = 'flowline_split_by_point.json'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
#export_flowline_to_json( aFlowline,pSpatialRef, sFilename_out)

#ues location to find outlet
point['x'] = -2136506.345
point['y'] = 2901799.219
pVertex_outlet=pyvertex(point)
aFlowline= correct_flowline_direction(aFlowline,  pVertex_outlet )

sFilename_out = 'flowline_direction.json'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
#export_flowline_to_json( aFlowline,pSpatialRef, sFilename_out)

#step 4: remove loops

aFlowline = remove_flowline_loop(  aFlowline)    
sFilename_out = 'flowline_loop.json'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
export_flowline_to_json( aFlowline,pSpatialRef, sFilename_out)

dThreshold = 1.0E4 
aFlowline = remove_small_river(aFlowline, dThreshold)
sFilename_out = 'flowline_large.json'
sFilename_out =os.path.join(sWorkspace_out, sFilename_out)
export_flowline_to_json( aFlowline, pSpatialRef, sFilename_out)

aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity = find_flowline_confluence(aFlowline,  pVertex_outlet)
sFilename_out = 'flowline_vertex_with_confluence.json'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
export_vertex_to_json( aVertex, pSpatialRef, sFilename_out, aAttribute_data=aConnectivity)

aFlowline = merge_flowline( aFlowline,aVertex, pVertex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence  )  
sFilename_out = 'flowline_merge.json'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
export_flowline_to_json( aFlowline, pSpatialRef, sFilename_out)

aFlowline = remove_small_river(aFlowline, dThreshold)
sFilename_out = 'flowline_large2.json'
sFilename_out =os.path.join(sWorkspace_out, sFilename_out)
export_flowline_to_json( aFlowline, pSpatialRef, sFilename_out)


aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity = find_flowline_confluence(aFlowline,  pVertex_outlet)
sFilename_out = 'flowline_vertex_with_confluence2.json'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
export_vertex_to_json( aVertex, pSpatialRef, sFilename_out, aAttribute_data=aConnectivity)

aFlowline = merge_flowline( aFlowline,aVertex, pVertex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence  )  
sFilename_out = 'flowline_merge2.json'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
export_flowline_to_json( aFlowline, pSpatialRef, sFilename_out)

#build segment index
aFlowline, aStream_segment = define_stream_segment_index(aFlowline)
sFilename_out = 'flowline_segment.json'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
export_flowline_to_json( aFlowline, pSpatialRef, sFilename_out, \
    aAttribute_data=[aStream_segment], aAttribute_field=['iseg'], aAttribute_dtype=['int'])
#build stream order 

aFlowline, aStream_order = define_stream_order(aFlowline)
sFilename_out = 'flowline_order.json'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
export_flowline_to_json( aFlowline, pSpatialRef, sFilename_out, \
    aAttribute_data=[aStream_order], aAttribute_field=['iord'], aAttribute_dtype=['int'])

#simplify flowline after intersect

print('Finished')