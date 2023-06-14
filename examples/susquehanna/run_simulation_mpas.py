import os, sys
from pathlib import Path
from os.path import realpath
#===================================
#set up workspace path
#===================================
sPath_parent = str(Path(__file__).parents[2]) # data is located two dir's up

sys.path.append(sPath_parent)
from pyflowline.pyflowline_read_model_configuration_file import pyflowline_read_model_configuration_file

sPath_data = realpath( sPath_parent +  '/data/susquehanna' )
sWorkspace_input =  str(Path(sPath_data)  /  'input')
sWorkspace_output=  str(Path(sPath_data)  /  'output')


#===================================
#you need to update this file based on your own case study
#===================================
sFilename_configuration_in = realpath( sPath_parent +  '/examples/susquehanna/pyflowline_susquehanna_mpas.json' )
if os.path.isfile(sFilename_configuration_in):
    pass
else:
    print('This configuration does not exist: ', sFilename_configuration_in )

#===================================
#setup case information
#===================================
iCase_index = 1
iFlag_visualization = 1
sMesh = 'mpas'
sDate='20230101'

oPyflowline = pyflowline_read_model_configuration_file(sFilename_configuration_in, 
    iCase_index_in=iCase_index, sDate_in=sDate)


oPyflowline.aBasin[0].dLatitude_outlet_degree=39.462000
oPyflowline.aBasin[0].dLongitude_outlet_degree=-76.009300

oPyflowline.setup()
if iFlag_visualization ==1:
    oPyflowline.plot(sFilename_output_in = 'filter_flowline.png', sVariable_in = 'flowline_filter' )
    pass

oPyflowline.flowline_simplification()

if iFlag_visualization == 1:
     
    aExtent_meander = [-76.5,-76.2, 41.6,41.9] 
    oPyflowline.plot( iFlag_title=1 ,sVariable_in='flowline_simplified' , sFilename_output_in = 'flowline_simplified.png', ) 

    oPyflowline.plot( iFlag_title=1 ,sVariable_in='flowline_simplified' , aExtent_in =aExtent_meander, sFilename_output_in = 'flowline_simplified_zoom.png', ) 

    pass
aCell = oPyflowline.mesh_generation()

if iFlag_visualization == 1:
    oPyflowline.plot( iFlag_title=1 ,sVariable_in='mesh', sFilename_output_in = 'mesh.png' ) 
    pass

oPyflowline.reconstruct_topological_relationship(aCell)

if iFlag_visualization == 1:
    oPyflowline.plot( iFlag_title=1, sVariable_in='overlap', sFilename_output_in = 'mesh_w_flowline.png',)
    pass
oPyflowline.export()

print('Finished')
