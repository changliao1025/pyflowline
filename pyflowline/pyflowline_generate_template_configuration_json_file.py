import os, sys
from pathlib import Path
from os.path import realpath
#use this function to generate an initial json file for hexwatershed
import json
#once it's generated, you can modify it and use it for different simulations
from pyflowline.classes.pycase import flowlinecase
from pyflowline.classes.basin import pybasin

sPath = str(Path(__file__).parent.resolve())
sPath_data = realpath(sPath +  '/../data/susquehanna/' )

def pyflowline_generate_basin_template_configuration_json_file(sFilename_basins_json, nBasin, sWorkspace_output):
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
        
        aConfig_basin['sFilename_dam'] = str(Path(sPath_data)  /  'ICoM_dams.csv')
        #"/qfs/people/liao313/data/hexwatershed/susquehanna/auxiliary/ICoM_dams.csv"
        aConfig_basin['sFilename_flowline_filter'] = str(Path(sPath_data)  /  'streamord7above.shp')
        #"/qfs/people/liao313/data/hexwatershed/susquehanna/vector/hydrology/streamord7above.shp"
        aConfig_basin['sFilename_flowline_raw'] = str(Path(sPath_data)  /  'allflowline.shp')
        #"/qfs/people/liao313/data/hexwatershed/susquehanna/vector/hydrology/allflowline.shp"
        aConfig_basin['sFilename_flowline_topo'] = str(Path(sPath_data)  /  'flowline.csv')
        #"/qfs/people/liao313/data/hexwatershed/susquehanna/auxiliary/flowline.csv"
        aConfig_basin['sWorkspace_output_basin'] = str(Path(sWorkspace_output) / sBasin )
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
def pyflowline_generate_template_configuration_json_file(sFilename_json):
   
    if os.path.exists(sFilename_json):         
        os.remove(sFilename_json)
    
    

    nBasin = 1

 
    #use a dict to initialize the class
    aConfig = {}
    aConfig['iFlag_flowline'] = 1
    aConfig['iFlag_merge_reach'] = 1
    
    aConfig['iFlag_resample_method'] = 2   
    aConfig['iFlag_use_shapefile_extent'] = 1 
     
    aConfig['nOutlet'] = nBasin
    aConfig['dResolution_degree'] = 0.5
    aConfig['dResolution_meter'] = 50000
    aConfig['dLongitude_left'] = -180
    aConfig['dLongitude_right'] = 180
    aConfig['dLatitude_bot'] = -90
    aConfig['dLatitude_top'] = 90
    aConfig['sFilename_model_configuration']  = sFilename_json 
    aConfig['sWorkspace_data'] = sPath_data
    #'/people/liao313/data'    
    aConfig['sWorkspace_project'] = 'pyflowline'
    #aConfig['sWorkspace_bin'] = '/people/liao313/bin'
    aConfig['sWorkspace_output'] = str(Path(sPath_data)  /  'output')
    #"/compyfs/liao313/04model/pyflowline/susquehanna"
    aConfig['sRegion'] = 'susquehanna'
    aConfig['sModel'] = 'pyflowline'
    aConfig['iCase_index'] = 1
   
    aConfig['sMesh_type'] = 'hexagon'
    aConfig['sJob'] = 'pyflowline'
    aConfig['sDate']= '20220110'
    aConfig['sFilename_mesh'] = str(Path(sPath_data)  /  'lnd_cull_mesh.nc')
    #"/qfs/people/liao313/data/icom/mesh/delaware_lnd_60_30_5_2_v2/lnd_cull_mesh.nc"    
    aConfig['flowline_info'] = 'flowline_info.json'
    aConfig['sFilename_mesh_info'] = 'mesh_info.json'
    aConfig['sFilename_elevation'] = 'elevation.json'
    #aConfig['sFilename_dem']  = '/qfs/people/liao313/data/hexwatershed/susquehanna/raster/dem/dem_ext.tif'    
    aConfig['sFilename_spatial_reference'] =  str(Path(sPath_data)  /  'boundary_proj.shp')
    #'/qfs/people/liao313/data/hexwatershed/susquehanna/vector/hydrology/boundary_proj.shp'
    

    oModel = flowlinecase(aConfig)

    #generate basin
    sDirname = os.path.dirname(sFilename_json)
    sFilename =  Path(sFilename_json).stem + '_basins.json'
    sFilename_basins_json = os.path.join(sDirname, sFilename)

    aBasin = pyflowline_generate_basin_template_configuration_json_file(sFilename_basins_json, nBasin, oModel.sWorkspace_output)

    oModel.aBasin = aBasin
    oModel.sFilename_basins = sFilename_basins_json
    oModel.export_config_to_json(sFilename_json)

    

    return oModel

