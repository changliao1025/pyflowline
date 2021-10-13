import os, sys
import numpy as np

from pyflowline.case.pyflowline_read_model_configuration_file import pyflowline_read_model_configuration_file

from pyflowline.operation.preprocess_flowline_op import preprocess_flowline_op

sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pyflowline/pyflowline/config/pyflowline_icom_mpas.json'
oPyflowline = pyflowline_read_model_configuration_file(sFilename_configuration_in)

preprocess_flowline_op(oPyflowline)

exit
#sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pyflowline/pyflowline/config/pyflowline_NHDPLUS_H_0204_HU4_mpas.json'
#
#
#sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pyflowline/pyflowline/config/pyflowline_NHDPLUS_H_0206_HU4_2_mpas.json'
#oPyflowline = pyflowline_read_model_configuration_file(sFilename_configuration_in)
##preprocess_flowline_op(oPyflowline)
#
#sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pyflowline/pyflowline/config/pyflowline_NHDPLUS_H_0206_HU4_3_mpas.json'
#oPyflowline = pyflowline_read_model_configuration_file(sFilename_configuration_in)
##preprocess_flowline_op(oPyflowline)
#
#sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pyflowline/pyflowline/config/pyflowline_NHDPLUS_H_0206_HU4_4_mpas.json'
#oPyflowline = pyflowline_read_model_configuration_file(sFilename_configuration_in)
##preprocess_flowline_op(oPyflowline)
#
#sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pyflowline/pyflowline/config/pyflowline_NHDPLUS_H_0206_HU4_5_mpas.json'
#oPyflowline = pyflowline_read_model_configuration_file(sFilename_configuration_in)
##preprocess_flowline_op(oPyflowline)
#
#sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pyflowline/pyflowline/config/pyflowline_NHDPLUS_H_0207_HU4_mpas.json'
#oPyflowline = pyflowline_read_model_configuration_file(sFilename_configuration_in)
##preprocess_flowline_op(oPyflowline)
#
#sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pyflowline/pyflowline/config/pyflowline_NHDPLUS_H_0208_HU4_1_mpas.json'
#oPyflowline = pyflowline_read_model_configuration_file(sFilename_configuration_in)
#preprocess_flowline_op(oPyflowline)
#
#sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pyflowline/pyflowline/config/pyflowline_NHDPLUS_H_0208_HU4_2_mpas.json'
#oPyflowline = pyflowline_read_model_configuration_file(sFilename_configuration_in)
#preprocess_flowline_op(oPyflowline)
#
#sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pyflowline/pyflowline/config/pyflowline_NHDPLUS_H_0208_HU4_3_mpas.json'
#oPyflowline = pyflowline_read_model_configuration_file(sFilename_configuration_in)
#preprocess_flowline_op(oPyflowline)

print('Finished')