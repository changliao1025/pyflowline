import os, sys
import numpy as np

from pyflowline.case.pyflowline_read_model_configuration_file import pyflowline_read_model_configuration_file

from pyflowline.operation.preprocess_flowline_op import preprocess_flowline_op

sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pyflowline/pyflowline/config/pyflowline_susquehanna_mpas.json'
#sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pyflowline/pyflowline/config/pyflowline_icom_mpas.json'
#sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pyflowline/pyflowline/config/pyflowline_global_mpas.json'
oPyflowline = pyflowline_read_model_configuration_file(sFilename_configuration_in)

preprocess_flowline_op(oPyflowline)

print('Finished')