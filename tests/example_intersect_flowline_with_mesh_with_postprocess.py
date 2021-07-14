import os

from pystream.case.pystream_read_model_configuration_file import pystream_read_model_configuration_file
from pystream.case.pycase import streamcase


from pystream.operation.intersect_flowline_with_mesh_with_postprocess_op import intersect_flowline_with_mesh_with_postprocess_op

sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pystream/pystream/config/pystream_susquehanna_mpas.xml'
aParameter = pystream_read_model_configuration_file(sFilename_configuration_in)
aParameter['sFilename_model_configuration'] = sFilename_configuration_in
oModel = streamcase(aParameter)


intersect_flowline_with_mesh_with_postprocess_op(oModel)



