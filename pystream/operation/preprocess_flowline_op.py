import os, sys

import numpy as np

from pystream.shared.vertex import pyvertex
from pyearth.system.define_global_variables import *

from pystream.format.read_flowline_shapefile import read_flowline_shapefile
from pystream.format.export_flowline_to_shapefile import export_flowline_to_shapefile
from pystream.format.export_vertex_to_shapefile import export_vertex_to_shapefile

from pystream.algorithm.connect.connect_disconnect_flowline import connect_disconnect_flowline
from pystream.algorithm.direction.correct_flowline_direction import correct_flowline_direction

#merge
from pystream.algorithm.merge.merge_flowline import merge_flowline

#split
from pystream.algorithm.split.split_flowline import split_flowline
from pystream.algorithm.split.find_flowline_confluence import find_flowline_confluence
from pystream.algorithm.split.find_flowline_vertex import find_flowline_vertex


from pystream.algorithm.loop.remove_flowline_loop import remove_flowline_loop
#
from pystream.algorithm.simplification.remove_small_river import remove_small_river

from pystream.algorithm.index.define_stream_order import define_stream_order
from pystream.algorithm.index.define_stream_segment_index import define_stream_segment_index


"""
prepare the flowline using multiple step approach
"""

def preprocess_flowline_op(oModel_in):
    
    #read shapefile and store information in the list
    iFlag_disconnected = oModel_in.iFlag_disconnected
    dThreshold = oModel_in.dThreshold_small_river

    sFilename_flowlinw_raw = oModel_in.sFilename_flowlinw_raw


    sWorkspace_simulation_flowline = oModel_in.sWorkspace_simulation_flowline
    aFlowline, pSpatialRef = read_flowline_shapefile(sFilename_flowlinw_raw)
    #we also need to save the spatial reference information for the output purpose

    #the flowline should not be in GCS because it cannot be used for distance directly,

    sFilename_out = 'flowline_before_intersect.shp'
    sFilename_out = os.path.join(sWorkspace_simulation_flowline, sFilename_out)
    export_flowline_to_shapefile( aFlowline, pSpatialRef, sFilename_out)

    if iFlag_disconnected ==1:
        #need a better way to include this capability
        
        #aVertex=list()
        #point= dict()
        #point['x'] = -1589612.188
        #point['y'] = 3068975.112
        #pVertex=pyvertex(point)
        #aVertex.append(pVertex)
        #point['x'] =  -1568732.491
        #point['y'] = 3064177.639
        #pVertex=pyvertex(point)
        #aVertex.append(pVertex)


        aThreshold = np.full(2, 300.0, dtype=float)
        #aFlowline = connect_disconnect_flowline(aFlowline, aVertex, aThreshold)
        sFilename_out = 'flowline_connect.json'
        sFilename_out = os.path.join(sWorkspace_simulation_flowline, sFilename_out)
        export_flowline_to_shapefile( aFlowline,pSpatialRef, sFilename_out)
    else:
        pass


    aVertex = find_flowline_vertex(aFlowline)
    sFilename_out = 'flowline_vertex_without_confluence_before_intersect.shp'
    sFilename_out = os.path.join(sWorkspace_simulation_flowline, sFilename_out)
    export_vertex_to_shapefile( aVertex,pSpatialRef, sFilename_out)

    aFlowline = split_flowline(aFlowline, aVertex)
    sFilename_out = 'flowline_split_by_point_before_intersect.shp'
    sFilename_out = os.path.join(sWorkspace_simulation_flowline, sFilename_out)
    export_flowline_to_shapefile( aFlowline,pSpatialRef, sFilename_out)

    #ues location to find outlet
  
    point= dict()
    point['x'] = oModel_in.dx_outlet
    point['y'] = oModel_in.dy_outlet
    pVertex_outlet=pyvertex(point)

    aFlowline= correct_flowline_direction(aFlowline,  pVertex_outlet )

    pVertex_outlet = aFlowline[0].pVertex_end

    sFilename_out = 'flowline_direction_before_intersect.shp'
    sFilename_out = os.path.join(sWorkspace_simulation_flowline, sFilename_out)
    export_flowline_to_shapefile( aFlowline,pSpatialRef, sFilename_out)

    #step 4: remove loops

    aFlowline = remove_flowline_loop(aFlowline)    
    sFilename_out = 'flowline_loop_before_intersect.shp'
    sFilename_out = os.path.join(sWorkspace_simulation_flowline, sFilename_out)
    export_flowline_to_shapefile( aFlowline,pSpatialRef, sFilename_out)

    
    aFlowline = remove_small_river(aFlowline, dThreshold)
    sFilename_out = 'flowline_large_step1_before_intersect.shp'
    sFilename_out =os.path.join(sWorkspace_simulation_flowline, sFilename_out)
    export_flowline_to_shapefile( aFlowline, pSpatialRef, sFilename_out)

    aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity = find_flowline_confluence(aFlowline,  pVertex_outlet)
    sFilename_out = 'flowline_vertex_with_confluence_step1_before_intersect.shp'
    sFilename_out = os.path.join(sWorkspace_simulation_flowline, sFilename_out)
    export_vertex_to_shapefile( aVertex, pSpatialRef, sFilename_out, aAttribute_data=aConnectivity)

    aFlowline = merge_flowline( aFlowline,aVertex, pVertex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence  )  
    sFilename_out = 'flowline_merge_step1_before_intersect.shp'
    sFilename_out = os.path.join(sWorkspace_simulation_flowline, sFilename_out)
    export_flowline_to_shapefile( aFlowline, pSpatialRef, sFilename_out)

    aFlowline = remove_small_river(aFlowline, dThreshold)
    sFilename_out = 'flowline_large_step2_before_intersect.shp'
    sFilename_out =os.path.join(sWorkspace_simulation_flowline, sFilename_out)
    export_flowline_to_shapefile( aFlowline, pSpatialRef, sFilename_out)


    aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity = find_flowline_confluence(aFlowline,  pVertex_outlet)
    sFilename_out = 'flowline_vertex_with_confluence_step2_before_intersect.shp'
    sFilename_out = os.path.join(sWorkspace_simulation_flowline, sFilename_out)
    export_vertex_to_shapefile( aVertex, pSpatialRef, sFilename_out, aAttribute_data=aConnectivity)

    aFlowline = merge_flowline( aFlowline,aVertex, pVertex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence  )  
    sFilename_out = 'flowline_merge_step2_before_intersect.shp'
    sFilename_out = os.path.join(sWorkspace_simulation_flowline, sFilename_out)
    export_flowline_to_shapefile( aFlowline, pSpatialRef, sFilename_out)

    #build segment index
    aFlowline, aStream_segment = define_stream_segment_index(aFlowline)
    sFilename_out = oModel_in.sFilename_flowline_segment_index_before_intersect
    sFilename_out = os.path.join(sWorkspace_simulation_flowline, sFilename_out)
    export_flowline_to_shapefile( aFlowline, pSpatialRef, sFilename_out, \
        aAttribute_data=[aStream_segment], aAttribute_field=['iseg'], aAttribute_dtype=['int'])

    #build stream order 
    aFlowline, aStream_order = define_stream_order(aFlowline)
    sFilename_out = oModel_in.sFilename_flowline_segment_order_before_intersect
    
    export_flowline_to_shapefile( aFlowline, pSpatialRef, sFilename_out, \
        aAttribute_data=[aStream_segment, aStream_order], aAttribute_field=['iseg','iord'], aAttribute_dtype=['int','int'])

    

    print('Finished')