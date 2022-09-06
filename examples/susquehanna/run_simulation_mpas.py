import os, sys
from pathlib import Path
from os.path import realpath
import argparse
import logging
for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)

logging.basicConfig(format='%(asctime)s %(message)s')
logging.warning('is the time Pyflowline simulation started.')

from pyflowline.classes.pycase import flowlinecase
from pyflowline.pyflowline_read_model_configuration_file import pyflowline_read_model_configuration_file

sDate='20220630'

dataPath = str(Path(__file__).parents[2]) # data is located two dir's up

sWorkspace_data = realpath( dataPath +  '/data/susquehanna' )
sWorkspace_input =  str(Path(sWorkspace_data)  /  'input')
sWorkspace_output=  str(Path(sWorkspace_data)  /  'output')

#an examples configuration file is provided with the repository, but you need to update this file based on your own case study
sMesh_type = 'mpas'
iCase_index = 1
sPath = str( Path().resolve() )

sFilename_configuration_in = realpath( sPath +  '/examples/susquehanna/pyflowline_susquehanna_mpas.json' )
   
oPyflowline = pyflowline_read_model_configuration_file(sFilename_configuration_in, \
    iCase_index_in=iCase_index, sDate_in=sDate)
oPyflowline.aBasin[0].dLatitude_outlet_degree=39.462000
oPyflowline.aBasin[0].dLongitude_outlet_degree=-76.009300
oPyflowline.setup()
oPyflowline.flowline_simplification()
aCell = oPyflowline.mesh_generation()
oPyflowline.reconstruct_topological_relationship(aCell)
oPyflowline.analyze()
oPyflowline.evaluate()
oPyflowline.export()

logging.basicConfig(format='%(asctime)s %(message)s')
logging.warning('is the time Pyflowline simulation finished.')
