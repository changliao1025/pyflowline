import os
from pathlib import Path
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
        sFilename_configuration_in = str(Path.cwd()) +  '/configurations/pyflowline_susquehanna_hexagon.json' 
        sFilename_configuration_in = str(Path.cwd()) +  '/configurations/pyflowline_susquehanna_mpas.json' 
        sFilename_configuration_in = str(Path.cwd()) +  '/configurations/pyflowline_susquehanna_latlon.json' 
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
#oPyflowline.preprocess_flowline()
#oPyflowline.plot(sVariable_in = 'flowline_simplified')


#oPyflowline.create_mesh()
oPyflowline.plot(sVariable_in = 'mesh')
#exit()
aExtent_in = [-76.0,-76.5, 39.5,40.0]
#oPyflowline.plot(sVariable_in = 'overlap_filter' )
#oPyflowline.plot(sVariable_in = 'overlap_filter',aExtent_in=aExtent_in )

oPyflowline.intersect_flowline_with_mesh()
oPyflowline.plot(sVariable_in = 'final')
oPyflowline.plot(sVariable_in = 'overlap')

print('Finished')
