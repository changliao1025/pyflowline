import os 
import datetime
import json
from pyflowline.classes.pycase import flowlinecase

pDate = datetime.datetime.today()
sDate_default = "{:04d}".format(pDate.year) + "{:02d}".format(pDate.month) + "{:02d}".format(pDate.day)

def pyflowline_read_model_configuration_file(sFilename_configuration_in,   \
    iFlag_standalone_in= None, \
        iFlag_use_mesh_dem_in=None, \
        iCase_index_in=None,   \
            dResolution_degree_in = None,  \
            dResolution_meter_in = None,  \
                sMesh_type_in = None,  \
                sModel_in = None, \
                    sDate_in = None,\
                    sWorkspace_output_in = None):
    """read a model configuration

    Args:
        sFilename_configuration_in (str): _description_
        iFlag_standalone_in (int, optional): _description_. Defaults to None.
        iFlag_use_mesh_dem_in (int, optional): _description_. Defaults to None.
        iCase_index_in (int, optional): _description_. Defaults to None.
        dResolution_degree_in (float, optional): _description_. Defaults to None.
        dResolution_meter_in (float, optional): _description_. Defaults to None.
        sMesh_type_in (str, optional): _description_. Defaults to None.
        sModel_in (str, optional): _description_. Defaults to None.
        sDate_in (str, optional): _description_. Defaults to None.
        sWorkspace_output_in (str, optional): _description_. Defaults to None.

    Returns:
        _type_: _description_
    """

    if not os.path.isfile(sFilename_configuration_in):
        print(sFilename_configuration_in + ' does not exist')
        return
    
    # Opening JSON file
    with open(sFilename_configuration_in) as json_file:
        aConfig = json.load(json_file)   
    
    if iCase_index_in is not None:        
        iCase_index = iCase_index_in
    else:       
        iCase_index = int( aConfig['iCase_index'])
    
    if iFlag_standalone_in is not None:        
        iFlag_standalone = iFlag_standalone_in
    else:       
        iFlag_standalone = int( aConfig['iFlag_standalone'])

    if iFlag_use_mesh_dem_in is not None:        
        iFlag_use_mesh_dem = iFlag_use_mesh_dem_in
    else:       
        iFlag_use_mesh_dem = int(aConfig['iFlag_use_mesh_dem'])

    if sMesh_type_in is not None:
        sMesh_type = sMesh_type_in
    else:
        sMesh_type = aConfig['sModel']
        pass
        
    if sModel_in is not None:
        sModel = sModel_in
    else:
        sModel = aConfig['sModel']
        pass

    if sDate_in is not None:
        sDate = sDate_in
    else:
        sDate = aConfig['sDate']
        pass

    if dResolution_degree_in is not None:
        dResolution_degree = dResolution_degree_in
    else:
        dResolution_degree = aConfig['dResolution_degree']
        pass

    if dResolution_meter_in is not None:
        dResolution_meter = dResolution_meter_in
    else:
        dResolution_meter = aConfig['dResolution_meter']
        pass

    if sWorkspace_output_in is not None:
        sWorkspace_output = sWorkspace_output_in
    else:
        sWorkspace_output = aConfig['sWorkspace_output']
        pass

    
    aConfig['iCase_index'] = iCase_index
    aConfig['iFlag_standalone'] = iFlag_standalone
    aConfig['iFlag_use_mesh_dem'] = iFlag_use_mesh_dem
    aConfig['dResolution_degree'] = dResolution_degree
    aConfig['dResolution_meter'] = dResolution_meter

    aConfig['sDate'] = sDate
    aConfig['sModel'] = sModel
    aConfig['sWorkspace_output'] = sWorkspace_output
   
    

    #based on global variable, a few variables are calculate once
    #calculate the modflow simulation period
    #https://docs.python.org/3/library/datetime.html#datetime-objects
   
   
    
    #aConfig
    
    #simulation
    
    oPyflowline = flowlinecase(aConfig)
   
    
    return oPyflowline