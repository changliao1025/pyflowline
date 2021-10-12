import os, sys
from pyearth.system.define_global_variables import *
from pyflowline.shared.vertex import pyvertex

from pyflowline.format.read_flowline_geojson import read_flowline_geojson
from pyflowline.format.export_flowline_to_shapefile import export_flowline_to_shapefile

from pyflowline.algorithm.intersect.intersect_flowline_with_mesh import intersect_flowline_with_mesh



def intersect_flowline_with_mesh_op(oPyflowline_in):

    iMesh_type = oPyflowline_in.iMesh_type
    
    nOutlet = oPyflowline_in.nOutlet
    sFilename_mesh=oPyflowline_in.sFilename_mesh
    sWorkspace_output = oPyflowline_in.sWorkspace_output

    for i in range(nOutlet):
        sBasin =  "{:03d}".format(i+1)    
        sWorkspace_output_basin = oPyflowline_in.sWorkspace_output + slash + sBasin
        Path(sWorkspace_output_basin).mkdir(parents=True, exist_ok=True)          
        pBasin = oPyflowline_in.aBasin[i]
        sFilename_flowline = pBasin.sFilename_flowline_segment_order_before_intersect
        sFilename_flowline_in = os.path.join(sWorkspace_output_basin, sFilename_flowline)
        sFilename_flowline_intersect = pBasin.sFilename_flowline_intersect
        sFilename_flowline_intersect_out = os.path.join(sWorkspace_output_basin, sFilename_flowline_intersect)

        aCell, aCell_intersect, aFlowline_intersect = intersect_flowline_with_mesh(iMesh_type, sFilename_mesh, \
            sFilename_flowline_in, sFilename_flowline_intersect_out)


    

    return







