#this is an exmaple to run all the stepsimport os, sys
import numpy as np

from pyflowline.classes.pycase import flowlinecase
from pyflowline.classes.pyflowline_read_model_configuration_file import pyflowline_read_model_configuration_file
from pyflowline.operation.full_op import full_op

sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pyflowline/pyflowline/config/pyflowline_susquehanna_mpas.json'
#sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pyflowline/pyflowline/config/pyflowline_icom_mpas.json'
#sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pyflowline/pyflowline/config/pyflowline_global_mpas.json'

oPyflowline = pyflowline_read_model_configuration_file(sFilename_configuration_in )

aCell = full_op(oPyflowline)