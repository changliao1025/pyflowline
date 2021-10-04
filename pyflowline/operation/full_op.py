

from pyflowline.operation.preprocess_flowline_op import preprocess_flowline_op
from pyflowline.operation.create_mesh_op import create_mesh_op

from pyflowline.operation.intersect_flowline_with_mesh_with_postprocess_op import intersect_flowline_with_mesh_with_postprocess_op

def full_op(oPyflowline_in):

    iFlag_simplification = oPyflowline_in.iFlag_simplification
    iFlag_create_mesh = oPyflowline_in.iFlag_create_mesh
    iFlag_intersect = oPyflowline_in.iFlag_intersect
    
    if iFlag_simplification == 1:
        preprocess_flowline_op(oPyflowline_in)

    if iFlag_create_mesh==1:
        create_mesh_op(oPyflowline_in)

    if iFlag_intersect ==1:
        
        intersect_flowline_with_mesh_with_postprocess_op(oPyflowline_in)


    