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
sFilename_configuration_in = realpath( sPath_parent +  '/data/susquehanna/input/pyflowline_susquehanna_mpas.json' )
if os.path.isfile(sFilename_configuration_in):
    pass
else:
    print('This configuration does not exist: ', sFilename_configuration_in )

#===================================
#setup case information
#===================================
iCase_index = 1
iFlag_simulation = 0
iFlag_visualization = 1
sMesh = 'mpas'
sDate='20230701'

oPyflowline = pyflowline_read_model_configuration_file(sFilename_configuration_in, 
    iCase_index_in=iCase_index, sDate_in=sDate)

#take a look at the model parameters

oPyflowline.print()

#now we can change the following model parameters
#there are two ways to change the model parameters
#use a function or assign a value directly

oPyflowline.change_model_parameter(sVariable_in='sWorkspace_output', sValue_in=sWorkspace_output)

#if you need to change a parameter for a basin instead of the whole model domain, use the iFlag_basin_in option, this will change all the basins
oPyflowline.change_model_parameter(sVariable_in='dLatitude_outlet_degree', sValue_in=39.462000, iFlag_basin_in=1)
oPyflowline.change_model_parameter(sVariable_in='dLongitude_outlet_degree', sValue_in=-76.009300, iFlag_basin_in=1)


#the second way is to assign a value directly
oPyflowline.aBasin[0].dLatitude_outlet_degree=39.462000
oPyflowline.aBasin[0].dLongitude_outlet_degree=-76.009300

#oPyflowline.setup()
if iFlag_visualization ==1:
    #oPyflowline.plot(sVariable_in = 'flowline_filter', sFilename_output_in = 'filter_flowline.png'  )
    pass

if iFlag_simulation == 1:
    #oPyflowline.flowline_simplification()
    pass

if iFlag_visualization == 1:
     
    aExtent_meander = [-76.5,-76.2, 41.6,41.9] 
    #oPyflowline.plot( sVariable_in='flowline_simplified' , sFilename_output_in = 'flowline_simplified.png' ) 
    #oPyflowline.plot( sVariable_in='flowline_simplified' , sFilename_output_in = 'flowline_simplified_zoom.png', aExtent_in =aExtent_meander ) 

    pass

if iFlag_simulation == 1:
    aCell = oPyflowline.mesh_generation()

if iFlag_visualization == 1:
    oPyflowline.plot( sVariable_in='mesh', sFilename_output_in = 'mesh.png' ) 
    pass

if iFlag_simulation == 1:
    oPyflowline.reconstruct_topological_relationship(aCell)

if iFlag_visualization == 1:
    oPyflowline.plot(  sVariable_in='overlap', sFilename_output_in = 'mesh_w_flowline.png',)
    pass

oPyflowline.export()

print('Finished')
