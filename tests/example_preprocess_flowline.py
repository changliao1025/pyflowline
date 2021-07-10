import os, sys
import numpy as np


from pystream.case.pycase import streamcase
from pystream.case.pystream_read_model_configuration_file import pystream_read_model_configuration_file
from pystream.operation.preprocess_flowline import preprocess_flowline




sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pystream/pystream/config/case_susquehanna_mpas.xml'
aParameter = pystream_read_model_configuration_file(sFilename_configuration_in)
aParameter['sFilename_model_configuration'] = sFilename_configuration_in
oModel = streamcase(aParameter)
preprocess_flowline(oModel)


print('Finished')