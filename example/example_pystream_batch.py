#this is an exmaple to run all the stepsimport os, sys
import numpy as np

from pystream.case.pycase import streamcase
from pystream.case.pystream_read_model_configuration_file import pystream_read_model_configuration_file
from pystream.operation.full_op import full_op




aResolution= [1000, 10000, 50000]
sDate_in = '20210909'
for i in range(0,3):

    #hexagon
    iCase_index = i * 3 + 1
    dReresolution_meter = aResolution[i]

    sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pystream/pystream/config/pystream_susquehanna_hexagon.json'
    oPystream = pystream_read_model_configuration_file(sFilename_configuration_in, \
        sDate_in = sDate_in,\
        iCase_index_in=iCase_index, dResolution_meter_in =dReresolution_meter )

    aCell = full_op(oPystream)


    #square
    iCase_index = iCase_index + 1
    sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pystream/pystream/config/pystream_susquehanna_square.json'
    oPystream = pystream_read_model_configuration_file(sFilename_configuration_in, \
        sDate_in = sDate_in,\
        iCase_index_in=iCase_index, dResolution_meter_in =dReresolution_meter )

    aCell = full_op(oPystream)

    #latlon
    iCase_index = iCase_index + 1
    sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pystream/pystream/config/pystream_susquehanna_latlon.json'
    oPystream = pystream_read_model_configuration_file(sFilename_configuration_in, \
        sDate_in = sDate_in,\
        iCase_index_in=iCase_index, dResolution_meter_in =dReresolution_meter )

    aCell = full_op(oPystream)


