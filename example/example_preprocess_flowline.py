import os, sys
import numpy as np

from pystream.case.pystream_read_model_configuration_file import pystream_read_model_configuration_file

from pystream.operation.preprocess_flowline_op import preprocess_flowline_op


sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pystream/pystream/config/pystream_NHDPLUS_H_0204_HU4_3_mpas.json'
sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pystream/pystream/config/pystream_susquehanna_mpas.json'
oPystream = pystream_read_model_configuration_file(sFilename_configuration_in)

preprocess_flowline_op(oPystream)


print('Finished')