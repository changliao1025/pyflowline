import os
from pyflowline.shared.vertex import pyvertex

from pyflowline.format.read_flowline_geojson import read_flowline_geojson
from pyflowline.format.export_flowline_to_shapefile import export_flowline_to_shapefile

from pyflowline.algorithm.intersect.intersect_flowline_with_mesh import intersect_flowline_with_mesh



def intersect_flowline_with_mesh_op(oPyflowline_in):

    iMesh_type = oPyflowline_in.iMesh_type
    sWorkspace_output = oPyflowline_in.sWorkspace_output
    nOutlet = oPyflowline_in.nOutlet

    sFilename_mesh=oPyflowline_in.sFilename_mesh
    

    for i in range(nOutlet):
        sBasin =  "{:03d}".format(i+1)              
        pBasin = oPyflowline_in.aBasin[i]
        sFilename_flowline = pBasin.sFilename_flowline_segment_order_before_intersect
        sFilename_flowline_in = os.path.join(sWorkspace_output, sFilename_flowline)
        sFilename_flowline_intersect = pBasin.sFilename_flowline_intersect
        sFilename_flowline_intersect_out = os.path.join(sWorkspace_output, sFilename_flowline_intersect)

        aCell, aCell_intersect, aFlowline_intersect = intersect_flowline_with_mesh(iMesh_type, sFilename_mesh, \
            sFilename_flowline_in, sFilename_flowline_intersect_out)


    

    return







