import sys
from pathlib import Path
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
pArgs = parser.parse_args()
if len(sys.argv) == 1:
    sMesh_type = 'mpas'
    iCase_index = 1
    dResolution_meter=50000
else:
    if len(sys.argv)==4:
        sMesh_type = pArgs.sMesh_type
        iCase_index = pArgs.iCase_index
        dResolution_meter=pArgs.dResolution_meter
    else:
        pass

iFlag_option = 2
if iFlag_option ==1:

    sPath = str(Path(__file__).parent.resolve())
    sFilename_configuration_in = sPath +  '/../tests/configurations/template.json' 
    
    oPyflowline = pyflowline_generate_template_configuration_json_file(sFilename_configuration_in)
    print(oPyflowline.tojson())
    #now you can customize the model object
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
        
        print(sFilename_configuration_in)
        oPyflowline = pyflowline_read_model_configuration_file(sFilename_configuration_in, \
            iCase_index_in=iCase_index, dResolution_meter_in=dResolution_meter)
        #print the case information in details
        print(oPyflowline.tojson())

#pyflowline can process multiple basins within one singel run
#the total number of basin is controlled by the nOutlet variable
#convert the raw flowline into geojson in WGS84 system        
#oPyflowline.convert_flowline_to_json()
#oPyflowline.plot(sVariable_in = 'flowline_filter_json')
#oPyflowline.flowline_simplification()
#oPyflowline.plot(sVariable_in = 'flowline_simplified')
#oPyflowline.plot_study_area()
#exit()
#oPyflowline.mesh_generation()
#oPyflowline.plot(sVariable_in = 'mesh')
#exit()
#aExtent_full = [-78.5,-75.5, 39.2,42.5]
#aExtent_zoom = [-76.0,-76.5, 39.5,40.0] #outlet
#aExtent_zoom = [-76.5,-76.2, 41.6,41.9] #meander
#aExtent_zoom = [-77.3,-76.5, 40.2,41.0] #braided
#aExtent_zoom = [-77.3,-76.5, 40.2,41.0] #confluence

oPyflowline.reconstruct_topological_relationship()
#oPyflowline.plot(sVariable_in = 'final')
#oPyflowline.plot(sVariable_in = 'overlap',aExtent_in=aExtent_full )

#replace conceptual flowline with real flowline length
oPyflowline.analyze()
#oPyflowline.evaluate()
oPyflowline.export()

#sFilename_dem_flowline ='/qfs/people/liao313/data/hexwatershed/susquehanna/vector/swat/swat10k.shp'
#oPyflowline.compare_with_raster_dem_method(sFilename_dem_flowline,aExtent_in=aExtent_zoom )



print('Finished')

logging.basicConfig(format='%(asctime)s %(message)s')
logging.warning('is the time Pyflowline simulation finished.')
