from pathlib import Path
from pyflowline.case.pycase import flowlinecase
from pyflowline.case.pyflowline_read_model_configuration_file import pyflowline_read_model_configuration_file

iFlag_option = 2

if iFlag_option ==1:

    oPyflowline=flowlinecase()
    oPyflowline.iCase_index = 1


else: 
    if iFlag_option ==2:

        #an example configuration file is provided with the repository, but you need to update this file based on your own case study
        
        #sFilename_configuration_in = str(Path.cwd()) +  '/pyflowline/config/pyflowline_susquehanna_hexagon.json' 
        sFilename_configuration_in = '/qfs/people/liao313/workspace/python/pyflowline/pyflowline/config/pyflowline_susquehanna_hexagon.json'
        print(sFilename_configuration_in)
        oPyflowline = pyflowline_read_model_configuration_file(sFilename_configuration_in)

        #print the case information in details
        print(oPyflowline)
        

from pyflowline.plot.pyflowline_plot_flowline import pyflowline_plot_flowline
pyflowline_plot_flowline(oPyflowline, sVariable_in = 'flowline_raw') 

from pyflowline.operation.preprocess_flowline_op import preprocess_flowline_op
preprocess_flowline_op(oPyflowline)

pyflowline_plot_flowline(oPyflowline, sVariable_in = 'flowline_simplified')

from pyflowline.operation.create_mesh_op import create_mesh_op
aCell = create_mesh_op(oPyflowline)

from pyflowline.plot.pyflowline_plot_mesh import pyflowline_plot_mesh
pyflowline_plot_mesh(oPyflowline)

from pyflowline.operation.intersect_flowline_with_mesh_with_postprocess_op import intersect_flowline_with_mesh_with_postprocess_op
intersect_flowline_with_mesh_with_postprocess_op(oPyflowline)

pyflowline_plot_flowline(oPyflowline, sVariable_in = 'flowline_final')

