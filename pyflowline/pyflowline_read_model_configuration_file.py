import os 
import datetime
import json
from pyflowline.classes.pycase import flowlinecase

pDate = datetime.datetime.today()
sDate_default = "{:04d}".format(pDate.year) + "{:02d}".format(pDate.month) + "{:02d}".format(pDate.day)

def pyflowline_read_model_configuration_file(sFilename_configuration_in,\
    iFlag_standalone_in= None,\
        iFlag_use_mesh_dem_in=None,\
     iCase_index_in=None, \
         dResolution_degree_in = None,\
         dResolution_meter_in = None,\
         sMesh_type_in = None, \
             sModel_in = None,\
                 sDate_in = None,\
                     sWorkspace_output_in = None):

    if not os.path.isfile(sFilename_configuration_in):
        print(sFilename_configuration_in + ' does not exist')
        return
    
    # Opening JSON file
    with open(sFilename_configuration_in) as json_file:
        data = json.load(json_file)   
    
    if iCase_index_in is not None:        
        iCase_index = iCase_index_in
    else:       
        iCase_index = int( data['iCase_index'])
    
    if iFlag_standalone_in is not None:        
        iFlag_standalone = iFlag_standalone_in
    else:       
        iFlag_standalone = int( data['iFlag_standalone'])

    if iFlag_use_mesh_dem_in is not None:        
        iFlag_use_mesh_dem = iFlag_use_mesh_dem_in
    else:       
        iFlag_use_mesh_dem = int(data['iFlag_use_mesh_dem'])

     if sMesh_type_in is not None:
        sMesh_type = sMesh_type_in
    else:
        sMesh_type = data['sModel']
        pass
        
    if sModel_in is not None:
        sModel = sModel_in
    else:
        sModel = data['sModel']
        pass

    if sDate_in is not None:
        sDate = sDate_in
    else:
        sDate = data['sDate']
        pass

    if dResolution_degree_in is not None:
        dResolution_degree = dResolution_degree_in
    else:
        dResolution_degree = data['dResolution_degree']
        pass

    if dResolution_meter_in is not None:
        dResolution_meter = dResolution_meter_in
    else:
        dResolution_meter = data['dResolution_meter']
        pass

    if sWorkspace_output_in is not None:
        sWorkspace_output = sWorkspace_output_in
    else:
        sWorkspace_output = data['sWorkspace_output']
        pass

    
    data['iCase_index'] = iCase_index
    data['iFlag_standalone'] = iFlag_standalone
    data['iFlag_use_mesh_dem'] = iFlag_use_mesh_dem
    data['dResolution_degree'] = dResolution_degree
    data['dResolution_meter'] = dResolution_meter

    data['sDate'] = sDate
    data['sModel'] = sModel
    data['sWorkspace_output'] = sWorkspace_output
   
    

    #based on global variable, a few variables are calculate once
    #calculate the modflow simulation period
    #https://docs.python.org/3/library/datetime.html#datetime-objects
   
   
    
    #data
    
    #simulation
    
    oPyflowline = flowlinecase(data)
   
    
    return oPyflowline