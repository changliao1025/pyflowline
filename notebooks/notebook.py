import os
from pathlib import Path

import logging
for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)

logging.basicConfig(format='%(asctime)s %(message)s')
logging.warning('is the Pyflowline simulation started.')

from pyflowline.classes.pycase import flowlinecase
from pyflowline.pyflowline_read_model_configuration_file import pyflowline_read_model_configuration_file

iFlag_option = 2
if iFlag_option ==1:
    oPyflowline=flowlinecase()
    oPyflowline.iCase_index = 1
else: 
    if iFlag_option == 2:
        #an example configuration file is provided with the repository, but you need to update this file based on your own case study
        #linux
  
        sPath = str(Path(__file__).parent.resolve())
        sFilename_configuration_in = sPath +  '/../configurations/pyflowline_susquehanna_hexagon.json' 
        sFilename_configuration_in = sPath +  '/../tests/configurations/pyflowline_susquehanna_mpas.json' 
        #sFilename_configuration_in = str(Path.cwd()) +  '/configurations/pyflowline_susquehanna_latlon.json' 
        #sFilename_configuration_in = str(Path.cwd()) +  '/configurations/pyflowline_susquehanna_square.json' 
        #mac
        #sFilename_configuration_in = '/Users/liao313/workspace/python/pyflowline/configurations/pyflowline_susquehanna_hexagon_mac.json'
        print(sFilename_configuration_in)
        oPyflowline = pyflowline_read_model_configuration_file(sFilename_configuration_in)
        #print the case information in details
        print(oPyflowline.tojson())

#pyflowline can process multiple basins within one singel run
#the total number of basin is controlled by the nOutlet variable
#convert the raw flowline into geojson in WGS84 system        
#oPyflowline.convert_flowline_to_json()
#oPyflowline.plot(sVariable_in = 'flowline_filter_json')
oPyflowline.flowline_simplification()
#oPyflowline.plot(sVariable_in = 'flowline_simplified')


#oPyflowline.mesh_generation()
#oPyflowline.plot(sVariable_in = 'mesh')
#exit()
aExtent_full = [-78.5,-75.5, 39.2,42.5]
#aExtent_zoom = [-76.0,-76.5, 39.5,40.0] #outlet
#aExtent_zoom = [-76.5,-76.2, 41.6,41.9] #meander
#aExtent_zoom = [-77.3,-76.5, 40.2,41.0] #braided
aExtent_zoom = [-77.3,-76.5, 40.2,41.0] #confluence

oPyflowline.reconstruct_topological_relationship()
#oPyflowline.plot(sVariable_in = 'final')
#oPyflowline.plot(sVariable_in = 'overlap',aExtent_in=aExtent_full )

oPyflowline.export()

#sFilename_dem_flowline ='/qfs/people/liao313/data/hexwatershed/susquehanna/vector/swat/swat10k.shp'
#oPyflowline.compare_with_raster_dem_method(sFilename_dem_flowline,aExtent_in=aExtent_zoom )

#oPyflowline.evaluate()

print('Finished')

logging.basicConfig(format='%(asctime)s %(message)s')
logging.warning('is the time Pyflowline simulation finished.')
