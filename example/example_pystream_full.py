#this is an exmaple to run all the stepsimport os, sys
import numpy as np

from pyflowline.case.pycase import flowlinecase
from pyflowline.case.pyflowline_read_model_configuration_file import pyflowline_read_model_configuration_file
from pyflowline.operation.full_op import full_op
sDate_in = '20210909'

iCase_index = 10

sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pyflowline/pyflowline/config/pyflowline_susquehanna_mpas.json'
opyflowline = pyflowline_read_model_configuration_file(sFilename_configuration_in, \
        sDate_in = sDate_in,\
        iCase_index_in=iCase_index )

aCell = full_op(opyflowline)