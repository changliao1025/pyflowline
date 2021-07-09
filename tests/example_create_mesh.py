import os, sys
import numpy as np


from pystream.case.pycase import streamcase
from pystream.case.pystream_read_model_configuration_file import pystream_read_model_configuration_file
from pystream.operation.create_mesh import create_mesh

##############################################################################################
#you only need to change the xml configuration file, which contains all the required information
##############################################################################################
sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pystream/pystream/config/case_susquehanna_hexagon.xml'
#sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pystream/pystream/config/case_susquehanna_square.xml'
#sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pystream/pystream/config/case_susquehanna_latlon.xml'
#sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pystream/pystream/config/case_susquehanna_mpas.xml'
#sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pystream/pystream/config/case_susquehanna_tin.xml'



aParameter = pystream_read_model_configuration_file(sFilename_configuration_in)
aParameter['sFilename_model_configuration'] = sFilename_configuration_in
oModel = streamcase(aParameter)
aCell = create_mesh(oModel)
