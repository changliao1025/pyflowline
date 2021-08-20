#this is an exmaple to run all the stepsimport os, sys
import numpy as np

from pystream.case.pycase import streamcase
from pystream.case.pystream_read_model_configuration_file import pystream_read_model_configuration_file
from pystream.operation.full_op import full_op


sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pystream/pystream/config/pystream_susquehanna_latlon.json'
oPystream = pystream_read_model_configuration_file(sFilename_configuration_in)

aCell = full_op(oPystream)