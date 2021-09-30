import os

from pyflowline.operation.intersect_flowline_with_mesh_op import intersect_flowline_with_mesh_op
from pyflowline.format.export_vertex_to_shapefile import export_vertex_to_shapefile
from pyflowline.case.pyflowline_read_model_configuration_file import pyflowline_read_model_configuration_file
from pyflowline.case.pycase import flowlinecase


sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pyflowline/pyflowline/config/case_susquehanna_mpas.xml'
aParameter = pyflowline_read_model_configuration_file(sFilename_configuration_in)
aParameter['sFilename_model_configuration'] = sFilename_configuration_in
oModel = flowlinecase(aParameter)





intersect_flowline_with_mesh_op(oModel)



