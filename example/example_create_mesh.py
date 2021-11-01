import os, sys
from pyflowline.format.export_mesh_info_to_json import export_mesh_info_to_json
import numpy as np


from pyflowline.case.pycase import flowlinecase
from pyflowline.case.pyflowline_read_model_configuration_file import pyflowline_read_model_configuration_file

from pyflowline.operation.create_mesh_op import create_mesh_op

##############################################################################################
#you only need to change the json configuration file, which contains all the required information
##############################################################################################
#sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pyflowline/pyflowline/config/pyflowline_susquehanna_hexagon.json'

#sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pyflowline/pyflowline/config/case_susquehanna_square.json'

sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pyflowline/pyflowline/config/pyflowline_susquehanna_mpas.json'

#sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pyflowline/pyflowline/config/case_susquehanna_tin.json'

oModel = pyflowline_read_model_configuration_file(sFilename_configuration_in)

aCell = create_mesh_op(oModel)
