import os, sys
from pystream.format.export_mesh_info_to_json import export_mesh_info_to_json
import numpy as np


from pystream.case.pycase import streamcase
from pystream.case.pystream_read_model_configuration_file import pystream_read_model_configuration_file

from pystream.operation.create_mesh_op import create_mesh_op

##############################################################################################
#you only need to change the json configuration file, which contains all the required information
##############################################################################################
#sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pystream/pystream/config/pystream_susquehanna_hexagon.json'

#sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pystream/pystream/config/case_susquehanna_square.json'

sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pystream/pystream/config/pystream_susquehanna_mpas.json'

#sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pystream/pystream/config/case_susquehanna_tin.json'

oModel = pystream_read_model_configuration_file(sFilename_configuration_in)

aCell = create_mesh_op(oModel)
