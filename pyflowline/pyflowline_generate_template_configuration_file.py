import os, sys
from pathlib import Path
from os.path import realpath
#use this function to generate an initial json file for hexwatershed
import json
#once it's generated, you can modify it and use it for different simulations
from pyflowline.classes.pycase import flowlinecase
from pyflowline.classes.basin import pybasin

def pyflowline_generate_basin_template_configuration_file(sFilename_basins_json, nBasin, sWorkspace_input_in, sWorkspace_output_in):
    """generate baisn configureation

    Args:
        sFilename_basins_json (str): the filename
        nBasin (int): the total number of basin
        sWorkspace_input_in (str): the input data path
        sWorkspace_output (str): the output path

    Returns:
        basin: a basin object
    """

    aBasin_out = list()
    for i in range(nBasin):
        sBasin =  "{:03d}".format(i+1)   
        aConfig_basin = {}
        aConfig_basin['iFlag_dam'] = 0
        aConfig_basin['iFlag_disconnected'] = 0
        aConfig_basin['lBasinID'] = i + 1
        aConfig_basin['dLatitude_outlet_degree'] = -180
        aConfig_basin['dLongitude_outlet_degree'] = 180
        aConfig_basin['dAccumulation_threshold'] = -90
        aConfig_basin['dThreshold_small_river'] = 90        
        aConfig_basin['sFilename_dam'] = str(Path(sWorkspace_input_in)  /  'ICoM_dams.csv')        
        aConfig_basin['sFilename_flowline_filter'] = str(Path(sWorkspace_input_in)  /  'streamord7above.shp')        
        aConfig_basin['sFilename_flowline_raw'] = str(Path(sWorkspace_input_in)  /  'allflowline.shp')        
        aConfig_basin['sFilename_flowline_topo'] = str(Path(sWorkspace_input_in)  /  'flowline.csv')        
        aConfig_basin['sWorkspace_output_basin'] = str(Path(sWorkspace_output_in) / sBasin )
        pBasin = pybasin(aConfig_basin)    
        aBasin_out.append(pBasin)        
        pass
        
    #export basin config to a file    
    sFilename_basins_json = sFilename_basins_json 
    with open(sFilename_basins_json, 'w', encoding='utf-8') as f:
        sJson = json.dumps([json.loads(ob.tojson()) for ob in aBasin_out], indent = 4)        
        f.write(sJson)    
        f.close()

    return aBasin_out

def pyflowline_generate_template_configuration_file(sFilename_json, \
    sWorkspace_input, sWorkspace_output, iFlag_standalone_in=None, \
        iFlag_use_mesh_dem_in=None, iFlag_use_shapefile_extent_in=None,  iCase_index_in=None,  dResolution_degree_in = None, dResolution_meter_in = None,   sDate_in = None,  sMesh_type_in = None,    sModel_in = None):
    """generate pyflowline config template file

    Args:
        sFilename_json (str): _description_
        sWorkspace_input (str): _description_
        iFlag_standalone_in (int, optional): _description_. Defaults to None.
        iFlag_use_mesh_dem_in (int, optional): _description_. Defaults to None.
        iFlag_use_shapefile_extent_in (int, optional): _description_. Defaults to None.
        iCase_index_in (int, optional): _description_. Defaults to None.
        dResolution_degree_in (float, optional): _description_. Defaults to None.
        dResolution_meter_in (float, optional): _description_. Defaults to None.
        sDate_in (str, optional): _description_. Defaults to None.
        sMesh_type_in (str, optional): _description_. Defaults to None.
        sModel_in (str, optional): _description_. Defaults to None.
        sWorkspace_output_in (str, optional): _description_. Defaults to None.

    Returns:
        _type_: _description_
    """
    if os.path.exists(sFilename_json):         
        os.remove(sFilename_json)

    if iCase_index_in is not None:        
        iCase_index = iCase_index_in
    else:       
        iCase_index = 1
    
    if iFlag_standalone_in is not None:        
        iFlag_standalone = iFlag_standalone_in
    else:       
        iFlag_standalone = 1

    if iFlag_use_mesh_dem_in is not None:        
        iFlag_use_mesh_dem = iFlag_use_mesh_dem_in
    else:       
        iFlag_use_mesh_dem = 0

    if iFlag_use_shapefile_extent_in is not None:        
        iFlag_use_shapefile_extent = iFlag_use_shapefile_extent_in
    else:       
        iFlag_use_shapefile_extent = 0

    if sMesh_type_in is not None:
        sMesh_type = sMesh_type_in
    else:
        sMesh_type = 'hexagon'
        pass
    if sDate_in is not None:
        sDate = sDate_in
    else:
        sDate = '20220202'
        pass
    
    
    nBasin = 1

 
    #use a dict to initialize the class
    aConfig = {}
    
    
    aConfig['iFlag_use_shapefile_extent'] = iFlag_use_shapefile_extent 
    aConfig['iFlag_use_mesh_dem'] = iFlag_use_mesh_dem
    aConfig['iFlag_save_mesh'] = 1
    aConfig['iFlag_flowline'] = 1
    aConfig['iFlag_simplification']=1
    aConfig['iFlag_create_mesh']=1
    aConfig['iFlag_intersect']=1
    aConfig['iFlag_resample_method']=1
    aConfig['iFlag_global']=0
    aConfig['iFlag_multiple_outlet']=0
    aConfig['iFlag_elevation_profile']=1
    aConfig['iFlag_rotation']=0
    aConfig['iFlag_stream_burning_topology']=1
    aConfig['iFlag_save_elevation']=1     
    aConfig['nOutlet'] = nBasin
    aConfig['dResolution_degree'] = 0.5
    aConfig['dResolution_meter'] = 50000
    aConfig['dLongitude_left'] = -180
    aConfig['dLongitude_right'] = 180
    aConfig['dLatitude_bot'] = -90
    aConfig['dLatitude_top'] = 90
    aConfig['sFilename_model_configuration']  = sFilename_json 
    aConfig['sWorkspace_input'] = sWorkspace_input    
    aConfig['sWorkspace_output'] = sWorkspace_output        
    aConfig['sRegion'] = 'susquehanna'
    aConfig['sModel'] = 'pyflowline'
    aConfig['iCase_index'] = iCase_index   
    aConfig['sMesh_type'] = sMesh_type
    aConfig['sJob'] = 'pyflowline'
    aConfig['sDate']= sDate
    #full path
    aConfig['sFilename_mesh_netcdf'] = str(Path(sWorkspace_input)  /  'lnd_cull_mesh.nc')
    aConfig['sFilename_spatial_reference'] =  str(Path(sWorkspace_input)  /  'boundary_proj.shp')   
    #relative path
    aConfig['flowline_info'] = 'flowline_info.json'
    aConfig['sFilename_mesh_info'] = 'mesh_info.json'
    aConfig['sFilename_elevation'] = 'elevation.json' 
    oModel = flowlinecase(aConfig)
    #generate basin
    sDirname = os.path.dirname(sFilename_json)
    sFilename =  Path(sFilename_json).stem + '_basins.json'
    sFilename_basins_json = os.path.join(sDirname, sFilename)   
    aBasin = pyflowline_generate_basin_template_configuration_file(sFilename_basins_json, nBasin, sWorkspace_input,  oModel.sWorkspace_output)
    oModel.aBasin = aBasin
    oModel.sFilename_basins = sFilename_basins_json
    oModel.export_config_to_json(sFilename_json)    

    return oModel

