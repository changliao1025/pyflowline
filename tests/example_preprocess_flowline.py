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
#sFilename_shapefile_in = '/qfs/people/liao313/data/hexwatershed/columbia_river_basin/vector/hydrology/crb_flowline.shp'

sFilename_shapefile_in = '/qfs/people/liao313/data/hexwatershed/susquehanna/vector/hydrology/flowline.shp'

sFilename_mesh = 'hexagon.json'

sWorkspace_out = '/compyfs/liao313/04model/pyhexwatershed/columbia_river_basin'
sWorkspace_out = '/people/liao313/tmp/susquehanna'

#sWorkspace_out = '/qfs/people/liao313/tmp/columbia_river_basin'



#read shapefile and store information in the list
aFlowline, pSpatialRef = read_flowline_shapefile(sFilename_shapefile_in)
#we also need to save the spatial reference information for the output purpose

sFilename_out = 'flowline.shp'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
#export_flowline_to_json( aFlowline, pSpatialRef, sFilename_out)

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
#aFlowline = connect_disconnect_flowline(aFlowline, aVertex, aThreshold)
sFilename_out = 'flowline_connect.json'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
#export_flowline_to_json( aFlowline,pSpatialRef, sFilename_out)


aVertex = find_flowline_vertex(aFlowline)
sFilename_out = 'flowline_vertex_without_confluence.shp'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
#export_vertex_to_json( aVertex,pSpatialRef, sFilename_out)

aFlowline = split_flowline(aFlowline, aVertex)
sFilename_out = 'flowline_split_by_point.shp'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
export_flowline_to_json( aFlowline,pSpatialRef, sFilename_out)

#ues location to find outlet
#crb
point['x'] = -2136506.345
point['y'] = 2901799.219

#susquehanna -76.00756,39.46286,1691654.168,2005654.678
point['x'] = 1691654.168
point['y'] = 2005654.678

pVertex_outlet=pyvertex(point)
aFlowline= correct_flowline_direction(aFlowline,  pVertex_outlet )

pVertex_outlet = aFlowline[0].pVertex_end

sFilename_out = 'flowline_direction.shp'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
export_flowline_to_json( aFlowline,pSpatialRef, sFilename_out)

#step 4: remove loops

aFlowline = remove_flowline_loop(  aFlowline)    
sFilename_out = 'flowline_loop.shp'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
export_flowline_to_json( aFlowline,pSpatialRef, sFilename_out)

dThreshold = 3.0E4 
aFlowline = remove_small_river(aFlowline, dThreshold)
sFilename_out = 'flowline_large.shp'
sFilename_out =os.path.join(sWorkspace_out, sFilename_out)
export_flowline_to_json( aFlowline, pSpatialRef, sFilename_out)

aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity = find_flowline_confluence(aFlowline,  pVertex_outlet)
sFilename_out = 'flowline_vertex_with_confluence.shp'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
export_vertex_to_json( aVertex, pSpatialRef, sFilename_out, aAttribute_data=aConnectivity)

aFlowline = merge_flowline( aFlowline,aVertex, pVertex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence  )  
sFilename_out = 'flowline_merge.shp'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
export_flowline_to_json( aFlowline, pSpatialRef, sFilename_out)

aFlowline = remove_small_river(aFlowline, dThreshold)
sFilename_out = 'flowline_large2.shp'
sFilename_out =os.path.join(sWorkspace_out, sFilename_out)
export_flowline_to_json( aFlowline, pSpatialRef, sFilename_out)


aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity = find_flowline_confluence(aFlowline,  pVertex_outlet)
sFilename_out = 'flowline_vertex_with_confluence2.shp'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
export_vertex_to_json( aVertex, pSpatialRef, sFilename_out, aAttribute_data=aConnectivity)

aFlowline = merge_flowline( aFlowline,aVertex, pVertex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence  )  
sFilename_out = 'flowline_merge2.shp'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
export_flowline_to_json( aFlowline, pSpatialRef, sFilename_out)

#build segment index
aFlowline, aStream_segment = define_stream_segment_index(aFlowline)
sFilename_out = 'flowline_segment.shp'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
export_flowline_to_json( aFlowline, pSpatialRef, sFilename_out, \
    aAttribute_data=[aStream_segment], aAttribute_field=['iseg'], aAttribute_dtype=['int'])
#build stream order 

aFlowline, aStream_order = define_stream_order(aFlowline)
sFilename_out = 'flowline_segment_order.shp'
sFilename_out = os.path.join(sWorkspace_out, sFilename_out)
export_flowline_to_json( aFlowline, pSpatialRef, sFilename_out, \
    aAttribute_data=[aStream_segment, aStream_order], aAttribute_field=['iseg','iord'], aAttribute_dtype=['int','int'])

#simplify flowline after intersect

print('Finished')