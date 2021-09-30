import os
from pyflowline.shared.vertex import pyvertex
from pyearth.gis.gdal.gdal_function import reproject_coordinates
from pyflowline.format.read_flowline_shapefile import read_flowline_shapefile
from pyflowline.format.read_mesh_shapefile import read_mesh_shapefile
from pyflowline.format.read_flowline_geojson import read_flowline_geojson

from pyflowline.format.export_flowline_to_shapefile import export_flowline_to_shapefile

from pyflowline.algorithm.intersect.intersect_flowline_with_mesh import intersect_flowline_with_mesh

from pyflowline.algorithm.simplification.remove_returning_flowline import remove_returning_flowline
from pyflowline.algorithm.simplification.remove_duplicate_flowline import remove_duplicate_flowline
from pyflowline.algorithm.simplification.remove_duplicate_edge import remove_duplicate_edge
from pyflowline.algorithm.direction.correct_flowline_direction import correct_flowline_direction
from pyflowline.algorithm.loop.remove_flowline_loop import remove_flowline_loop
from pyflowline.algorithm.split.find_flowline_vertex import find_flowline_vertex
from pyflowline.algorithm.split.find_flowline_confluence import find_flowline_confluence
from pyflowline.algorithm.split.split_flowline import split_flowline
from pyflowline.algorithm.split.split_flowline_to_edge import split_flowline_to_edge
from pyflowline.format.export_vertex_to_shapefile import export_vertex_to_shapefile
from pyflowline.algorithm.merge.merge_flowline import merge_flowline

from pyflowline.algorithm.index.define_stream_order import define_stream_order
from pyflowline.algorithm.index.define_stream_segment_index import define_stream_segment_index

def intersect_flowline_with_mesh_with_postprocess_op(opyflowline_in):

    #important
    

    iMesh_type = opyflowline_in.iMesh_type
    
    iFlag_projected = 0
    


    sWorkspace_output = opyflowline_in.sWorkspace_output  

    sFilename_flowline_filter = opyflowline_in.sFilename_flowline_filter


    aFlowline, pSpatialRef_flowline = read_flowline_shapefile(sFilename_flowline_filter)
    

    sFilename_flowline = opyflowline_in.sFilename_flowline_segment_order_before_intersect

    sFilename_mesh=opyflowline_in.sFilename_mesh
    aMesh, pSpatialRef_mesh = read_mesh_shapefile(sFilename_mesh)
    sFilename_flowline_intersect = opyflowline_in.sFilename_flowline_intersect

    
    aCell, aCell_intersect, aFlowline_intersect_all = intersect_flowline_with_mesh(\
        iMesh_type, sFilename_mesh, sFilename_flowline, sFilename_flowline_intersect)


    point= dict()
    
    point['lon'] = opyflowline_in.dLon_outlet
    point['lat'] = opyflowline_in.dLat_outlet
    pVertex_outlet=pyvertex(point)
    
    aFlowline, aFlowline_no_parallel, lCellID_outlet = remove_returning_flowline(iMesh_type, aCell_intersect, pVertex_outlet)
    sFilename_out = 'flowline_simplified_after_intersect.shp'
    sFilename_out = os.path.join(sWorkspace_output, sFilename_out)  
    
    pSpatialRef=  pSpatialRef_mesh
       
    export_flowline_to_shapefile(iFlag_projected, aFlowline, pSpatialRef, sFilename_out)

    #added start
    aFlowline, aEdge = split_flowline_to_edge(aFlowline)
    
    aFlowline = remove_duplicate_flowline(aFlowline)

    sFilename_out = 'flowline_debug.shp'
    sFilename_out = os.path.join(sWorkspace_output, sFilename_out)
    export_flowline_to_shapefile(iFlag_projected, aFlowline, pSpatialRef, sFilename_out)

    aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity\
        = find_flowline_confluence(aFlowline,  pVertex_outlet)

    aFlowline = merge_flowline( aFlowline,aVertex, pVertex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence  )  

    aFlowline = remove_flowline_loop(  aFlowline )    

    aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity\
        = find_flowline_confluence(aFlowline,  pVertex_outlet)

    aFlowline = merge_flowline( aFlowline,aVertex, pVertex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence  ) 
    #added end


    
    
    #pVertex_outlet=aFlowline[0].pVertex_end
    #aVertex = find_flowline_vertex(aFlowline)
    #
    #sFilename_out = 'flowline_vertex_without_confluence_after_intersect.shp'
    #sFilename_out = os.path.join(sWorkspace_output, sFilename_out)
    #export_vertex_to_shapefile(iFlag_projected, aVertex, pSpatialRef, sFilename_out)
    #
    #aFlowline = split_flowline(aFlowline, aVertex)
    #sFilename_out = 'flowline_split_by_point_after_intersect.shp'
    #sFilename_out = os.path.join(sWorkspace_output, sFilename_out)
    #export_flowline_to_shapefile(iFlag_projected, aFlowline, pSpatialRef, sFilename_out)
    #aFlowline= correct_flowline_direction(aFlowline,  pVertex_outlet )
#
#
    #
    #aFlowline = remove_flowline_loop(  aFlowline )    
    #sFilename_out = 'flowline_remove_loop_after_intersect.shp'
    #sFilename_out = os.path.join(sWorkspace_output, sFilename_out)
    #export_flowline_to_shapefile(iFlag_projected, aFlowline, pSpatialRef, sFilename_out)
#
#
    #aFlowline, aEdge = split_flowline_to_edge(aFlowline)
    ##aEdge = remove_duplicate_edge(aEdge)
    #aFlowline = remove_duplicate_flowline(aFlowline)
#
    #aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity\
    #    = find_flowline_confluence(aFlowline,  pVertex_outlet)
#
    #aFlowline = merge_flowline( aFlowline,aVertex, pVertex_outlet, aIndex_headwater,#aIndex_middle, aIndex_confluence  )  
    aFlowline, aStream_segment = define_stream_segment_index(aFlowline)
    aFlowline, aStream_order = define_stream_order(aFlowline)
    
    sFilename_out = 'flowline_final.shp'
    sFilename_out = os.path.join(sWorkspace_output, sFilename_out)
    export_flowline_to_shapefile(iFlag_projected, aFlowline, pSpatialRef, sFilename_out)

    return aCell, aCell_intersect, aFlowline, lCellID_outlet








