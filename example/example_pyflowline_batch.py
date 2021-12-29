#this is an exmaple to run all the stepsimport os, sys
import numpy as np

from pyflowline.classes.pycase import flowlinecase
from pyflowline.classes.pyflowline_read_model_configuration_file import pyflowline_read_model_configuration_file
from pyflowline.operation.full_op import full_op




aResolution= [1000, 10000, 50000]
sDate_in = '20210909'
for i in range(0,3):

    #hexagon
    iCase_index = i * 3 + 1
    dReresolution_meter = aResolution[i]

    sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pyflowline/pyflowline/config/pyflowline_susquehanna_hexagon.json'
    oPyflowline = pyflowline_read_model_configuration_file(sFilename_configuration_in, \
        sDate_in = sDate_in,\
        iCase_index_in=iCase_index, dResolution_meter_in =dReresolution_meter )

    aCell = full_op(oPyflowline)


    #square
    iCase_index = iCase_index + 1
    sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pyflowline/pyflowline/config/pyflowline_susquehanna_square.json'
    oPyflowline = pyflowline_read_model_configuration_file(sFilename_configuration_in, \
        sDate_in = sDate_in,\
        iCase_index_in=iCase_index, dResolution_meter_in =dReresolution_meter )

    aCell = full_op(oPyflowline)

    #latlon
    iCase_index = iCase_index + 1
    sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pyflowline/pyflowline/config/pyflowline_susquehanna_latlon.json'
    oPyflowline = pyflowline_read_model_configuration_file(sFilename_configuration_in, \
        sDate_in = sDate_in,\
        iCase_index_in=iCase_index, dResolution_meter_in =dReresolution_meter )

    aCell = full_op(oPyflowline)


