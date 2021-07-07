import os, sys
import numpy as np
import osgeo
from pystream.case.pycase import streamcase
from pystream.case.pystream_read_model_configuration_file import pystream_read_model_configuration_file
from pystream.mesh.jigsaw.create_mpas_mesh import create_mpas_mesh

sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pystream/pystream/config/case_susquehanna_mpas.xml'
aParameter = pystream_read_model_configuration_file(sFilename_configuration_in)
aParameter['sFilename_model_configuration'] = sFilename_configuration_in
oModel = streamcase(aParameter)




sFilename_spatial_reference = oModel.sFilename_spatial_reference 
#'/qfs/people/liao313/data/hexwatershed/susquehanna/vector/hydrology/stream_order78.shp'

aMpas = create_mpas_mesh(oModel)

print('finished')