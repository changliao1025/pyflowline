import os
from pystream.shared.vertex import pyvertex

from pystream.format.read_flowline_geojson import read_flowline_geojson
from pystream.format.export_flowline_to_shapefile import export_flowline_to_shapefile

from pystream.algorithm.intersect.intersect_flowline_with_mesh import intersect_flowline_with_mesh



def intersect_flowline_with_mesh_op(oModel_in):

    
    sWorkspace_simulation_case = oModel_in.sWorkspace_simulation_case  



    sFilename_flowline = oModel_in.sFilename_flowline_segment_order_before_intersect

    sFilename_mesh=oModel_in.sFilename_mesh
    sFilename_flowline_intersect = oModel_in.sFilename_flowline_intersect

    
    aCell, aCell_intersect, aFlowline_intersect = intersect_flowline_with_mesh(4, sFilename_mesh, sFilename_flowline, sFilename_flowline_intersect)


    

    return







