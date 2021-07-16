import os
from pystream.shared.vertex import pyvertex

from pystream.format.read_flowline_shapefile import read_flowline_shapefile
from pystream.format.read_flowline_geojson import read_flowline_geojson
from pystream.format.export_flowline_to_shapefile import export_flowline_to_shapefile

from pystream.algorithm.intersect.intersect_flowline_with_mesh import intersect_flowline_with_mesh

from pystream.algorithm.simplification.remove_returning_flowline import remove_returning_flowline


from pystream.algorithm.direction.correct_flowline_direction import correct_flowline_direction
from pystream.algorithm.loop.remove_flowline_loop import remove_flowline_loop
from pystream.algorithm.split.find_flowline_vertex import find_flowline_vertex
from pystream.algorithm.split.split_flowline import split_flowline
from pystream.format.export_vertex_to_shapefile import export_vertex_to_shapefile


def intersect_flowline_with_mesh_with_postprocess_op(oModel_in):

    iMesh_type = oModel_in.iMesh_type

    sWorkspace_simulation_case = oModel_in.sWorkspace_simulation_case  

    sFilename_flowlinw_raw = oModel_in.sFilename_flowlinw_raw


    sWorkspace_simulation_case = oModel_in.sWorkspace_simulation_case
    aFlowline, pSpatialRef = read_flowline_shapefile(sFilename_flowlinw_raw)

    sFilename_flowline = oModel_in.sFilename_flowline_segment_order_before_intersect

    sFilename_mesh=oModel_in.sFilename_mesh
    sFilename_flowline_intersect = oModel_in.sFilename_flowline_intersect

    
    aCell, aCell_intersect, aFlowline_intersect = intersect_flowline_with_mesh(iMesh_type, sFilename_mesh, sFilename_flowline, sFilename_flowline_intersect)


    point= dict()
    point['x'] = oModel_in.dx_outlet
    point['y'] = oModel_in.dy_outlet
    pVertex_outlet=pyvertex(point)
    
    aCell, aFlowline, aFlowline_no_parallel = remove_returning_flowline(iMesh_type, aCell, aCell_intersect, pVertex_outlet)
    sFilename_out = 'flowline_simplified_after_intersect.shp'
    sFilename_out = os.path.join(sWorkspace_simulation_case, sFilename_out)    
    export_flowline_to_shapefile( aFlowline, pSpatialRef, sFilename_out)
    
    
    pVertex_outlet=aFlowline[0].pVertex_end
    aVertex = find_flowline_vertex(aFlowline)
    
    sFilename_out = 'flowline_vertex_without_confluence_after_intersect.shp'
    sFilename_out = os.path.join(sWorkspace_simulation_case, sFilename_out)
    export_vertex_to_shapefile( aVertex,pSpatialRef, sFilename_out)
    
    aFlowline = split_flowline(aFlowline, aVertex)
    sFilename_out = 'flowline_split_by_point_after_intersect.shp'
    sFilename_out = os.path.join(sWorkspace_simulation_case, sFilename_out)
    export_flowline_to_shapefile( aFlowline,pSpatialRef, sFilename_out)
    
    aFlowline= correct_flowline_direction(aFlowline,  pVertex_outlet )
    
    aFlowline = remove_flowline_loop(  aFlowline )    
    sFilename_out = 'flowline_after_intersect.shp'
    sFilename_out = os.path.join(sWorkspace_simulation_case, sFilename_out)
    export_flowline_to_shapefile( aFlowline,pSpatialRef, sFilename_out)

    return aFlowline








