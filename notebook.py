import os
from pathlib import Path
from pyflowline.case.pycase import flowlinecase
from pyflowline.case.pyflowline_read_model_configuration_file import pyflowline_read_model_configuration_file

iFlag_option = 2

if iFlag_option ==1:

    oPyflowline=flowlinecase()
    oPyflowline.iCase_index = 1


else: 
    if iFlag_option == 2:

        #an example configuration file is provided with the repository, but you need to update this file based on your own case study
        
        #sFilename_configuration_in = str(Path.cwd()) +  '/pyflowline/config/pyflowline_susquehanna_hexagon.json' 
        sFilename_configuration_in = '/Users/liao313/workspace/python/pyflowline/pyflowline/config/pyflowline_susquehanna_hexagon_mac.json'
        print(sFilename_configuration_in)
        oPyflowline = pyflowline_read_model_configuration_file(sFilename_configuration_in)

        #print the case information in details
        print(oPyflowline)

#pyflowline can process multiple basins within one singel run
#the total number of basin is controlled by the nOutlet variable
#convert the raw flowline into geojson in WGS84 system        
from pyflowline.format.convert_shapefile_to_json import convert_shapefile_to_json
nOutlet = oPyflowline.nOutlet
for i in range(nOutlet):
    sBasin =  "{:03d}".format(i+1)    
    sWorkspace_output_basin = Path(oPyflowline.sWorkspace_output) / sBasin
    Path(sWorkspace_output_basin).mkdir(parents=True, exist_ok=True)                      
    pBasin = oPyflowline.aBasin[i]
    #the original flowline in shapefile format
    sFilename_raw = pBasin.sFilename_flowline_filter
    #the new flowine in geojson format in WGS84
    
    sFilename_out = pBasin.sFilename_flowline_filter_json
    #convert_shapefile_to_json(sFilename_raw, sFilename_out)
    pass


#now we can visualize the flowline
from pyflowline.plot.pyflowline_plot_flowline import pyflowline_plot_flowline
for i in range(nOutlet):
    pBasin = oPyflowline.aBasin[i]
    #pyflowline_plot_flowline(pBasin, sVariable_in = 'flowline_filter_json') 
    pass

from pyflowline.operation.preprocess_flowline_op import preprocess_flowline_op
#preprocess_flowline_op(oPyflowline)

for i in range(nOutlet):
    pBasin = oPyflowline.aBasin[i]
    #pyflowline_plot_flowline(pBasin, sVariable_in = 'flowline_simplified')
    pass

from pyflowline.operation.create_mesh_op import create_mesh_op
aCell = create_mesh_op(oPyflowline)

from pyflowline.plot.pyflowline_plot_mesh import pyflowline_plot_mesh
pyflowline_plot_mesh(oPyflowline)



from pyflowline.operation.intersect_flowline_with_mesh_with_postprocess_op import intersect_flowline_with_mesh_with_postprocess_op
intersect_flowline_with_mesh_with_postprocess_op(oPyflowline)

for i in range(nOutlet):
    pBasin = oPyflowline.aBasin[i]
    pyflowline_plot_flowline(pBasin, sVariable_in = 'flowline_final')
    pass

