import os, sys
from pystream.format.export_mesh_info_to_json import export_mesh_info_to_json
import numpy as np


from pystream.case.pycase import streamcase
from pystream.case.pystream_read_model_configuration_file import pystream_read_model_configuration_file

from pystream.operation.create_mesh_op import create_mesh_op

##############################################################################################
#you only need to change the xml configuration file, which contains all the required information
##############################################################################################
sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pystream/pystream/config/pystream_susquehanna_mpas.xml'
aParameter = pystream_read_model_configuration_file(sFilename_configuration_in)
aParameter['sFilename_model_configuration'] = sFilename_configuration_in
oModel = streamcase(aParameter)
aCell = create_mesh_op(oModel)

#export_mesh_info_to_json(aCell, sFilename_json_out=oModel.sFilename_mesh_info)
#
#sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pystream/pystream/config/case_susquehanna_square.xml'
#aParameter = pystream_read_model_configuration_file(sFilename_configuration_in)
#aParameter['sFilename_model_configuration'] = sFilename_configuration_in
#oModel = streamcase(aParameter)
#aCell = create_mesh_op(oModel)
#
#sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pystream/pystream/config/case_susquehanna_latlon.xml'
#aParameter = pystream_read_model_configuration_file(sFilename_configuration_in)
#aParameter['sFilename_model_configuration'] = sFilename_configuration_in
#oModel = streamcase(aParameter)
#aCell = create_mesh_op(oModel)

#sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pystream/pystream/config/pystream_susquehanna_mpas.xml'
#aParameter = pystream_read_model_configuration_file(sFilename_configuration_in)
#aParameter['sFilename_model_configuration'] = sFilename_configuration_in
#oModel = streamcase(aParameter)
#aCell = create_mesh_op(oModel)

#sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pystream/pystream/config/case_susquehanna_tin.xml'
#aParameter = pystream_read_model_configuration_file(sFilename_configuration_in)
#aParameter['sFilename_model_configuration'] = sFilename_configuration_in
#oModel = streamcase(aParameter)
#aCell = create_mesh_op(oModel)
