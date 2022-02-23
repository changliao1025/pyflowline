import sys
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
from pyflowline.pyflowline_generate_template_configuration_json_file import pyflowline_generate_template_configuration_json_file

parser = argparse.ArgumentParser()
parser.add_argument("--sMesh_type", help = "sMesh_type",  type = str)
parser.add_argument("--iCase_index", help = "iCase_index",  type = int)
parser.add_argument("--dResolution_meter", help = "dResolution_meter",  type = float)
parser.add_argument("--sDate", help = "sDate",  type = str)

#example
#python notebook.py --sMesh_type hexagon --iCase_index 1 --dResolution_meter 50000 --sDate 20220201
pArgs = parser.parse_args()
if len(sys.argv) == 1:
    sMesh_type = 'mpas'
    iCase_index = 10
    dResolution_meter=5000
    sDate='20220201'
else:
    if len(sys.argv)> 1:
        sMesh_type = pArgs.sMesh_type
        iCase_index = pArgs.iCase_index
        dResolution_meter=pArgs.dResolution_meter
        sDate = pArgs.sDate
        print(sMesh_type, iCase_index, dResolution_meter, sDate)
    else:
        print(len(sys.argv), 'Missing arguaments')
        pass

iFlag_option = 1
if iFlag_option ==1:
    sPath = str(Path(__file__).parent.resolve())
    sFilename_configuration_in = sPath +  '/../configurations/template.json' 
    oPyflowline = pyflowline_generate_template_configuration_json_file(sFilename_configuration_in)
    print(oPyflowline.tojson())
    #now you can customize the model object again
    oPyflowline.iCase_index = 1
    print(oPyflowline.tojson())
else: 
    if iFlag_option == 2:
        #an example configuration file is provided with the repository, but you need to update this file based on your own case study
        #linux
        sPath = str(Path(__file__).parent.resolve())
        if sMesh_type=='hexagon':
            sFilename_configuration_in = sPath +  '/../tests/configurations/pyflowline_susquehanna_hexagon.json' 
        else:
            if sMesh_type=='square':
                sFilename_configuration_in = sPath +  '/../tests/configurations/pyflowline_susquehanna_square.json' 
            else:
                if sMesh_type=='latlon':
                    sFilename_configuration_in = sPath +  '/../tests/configurations/pyflowline_susquehanna_latlon.json' 
                else:
                    sFilename_configuration_in = sPath +  '/../tests/configurations/pyflowline_susquehanna_mpas.json' 
        
        
        oPyflowline = pyflowline_read_model_configuration_file(sFilename_configuration_in, \
            iCase_index_in=iCase_index, dResolution_meter_in=dResolution_meter, sDate_in=sDate)
        #print the case information in details
        print(oPyflowline.tojson())

#pyflowline can process multiple basins within one singel run
#the total number of basin is controlled by the nOutlet variable
#convert the raw flowline into geojson in WGS84 system        
oPyflowline.convert_flowline_to_json()

oPyflowline.flowline_simplification()

oPyflowline.mesh_generation()

oPyflowline.reconstruct_topological_relationship()

oPyflowline.analyze()

oPyflowline.export()

print('Finished')

logging.basicConfig(format='%(asctime)s %(message)s')
logging.warning('is the time Pyflowline simulation finished.')
