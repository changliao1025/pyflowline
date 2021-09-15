#this is an exmaple to run all the stepsimport os, sys
import numpy as np

from pystream.case.pycase import streamcase
from pystream.case.pystream_read_model_configuration_file import pystream_read_model_configuration_file
from pystream.operation.full_op import full_op
sDate_in = '20210909'

iCase_index = 10

sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pystream/pystream/config/pystream_susquehanna_mpas.json'
oPystream = pystream_read_model_configuration_file(sFilename_configuration_in, \
        sDate_in = sDate_in,\
        iCase_index_in=iCase_index )

aCell = full_op(oPystream)