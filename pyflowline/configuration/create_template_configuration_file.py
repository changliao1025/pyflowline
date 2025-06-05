import os
from pathlib import Path
#use this function to generate an initial json file for hexwatershed
import json
#once it's generated, you can modify it and use it for different simulations
from pyflowline.classes.pycase import flowlinecase
from pyflowline.classes.basin import pybasin

def create_template_basin_configuration_file(
        sFilename_basins_json,
        nBasin,
        sWorkspace_input_in,
        sWorkspace_output_in):
    """generate basin configuration

    Args:
        sFilename_basins_json (str or Path): the filename
        nBasin (int): the total number of basin
        sWorkspace_input_in (str or Path): the input data path
        sWorkspace_output (str or Path): the output path

    Returns:
        basin: a basin object
    """

    # Ensure input pathnames are strings
    sFilename_basins_json = str(sFilename_basins_json)
    sWorkspace_input_in = str(sWorkspace_input_in)
    sWorkspace_output_in = str(sWorkspace_output_in)

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
    with open(sFilename_basins_json, 'w', encoding='utf-8') as f:
        sJson = json.dumps([json.loads(ob.tojson()) for ob in aBasin_out], indent = 4)
        f.write(sJson)
        f.close()

    return aBasin_out

def create_template_configuration_file(
        sFilename_json,
        sWorkspace_input=None,
        sWorkspace_output=None,
        iFlag_standalone_in=None,
        iFlag_use_mesh_dem_in=None,
        iCase_index_in=None,
        dResolution_degree_in=None,
        dResolution_meter_in=None,
        sDate_in=None,
        sMesh_type_in=None,
        sModel_in=None,
        iFlag_save_mesh_in=None,
        iFlag_simplification_in=None,
        iFlag_create_mesh_in=None,
        iFlag_intersect_in=None,
        iFlag_resample_method_in=None,
        iFlag_global_in=None,
        iFlag_multiple_outlet_in=None,
        iFlag_elevation_profile_in=None,
        iFlag_rotation_in=None,
        iFlag_stream_burning_topology_in=None,
        iFlag_save_elevation_in=None,
        nOutlet_in=None,
        dLongitude_left_in=None,
        dLongitude_right_in=None,
        dLatitude_bot_in=None,
        dLatitude_top_in=None,
        sRegion_in=None,
        sJob_in=None,
        flowline_info_in=None,
        sFilename_mesh_info_in=None,
        sFilename_elevation_in=None):
    """generate pyflowline config template file

    Args:
        sFilename_json (str or Path): Path to output configuration JSON file.
        sWorkspace_input (str or Path, optional): Input workspace path. Defaults to current directory.
        sWorkspace_output (str or Path, optional): Output workspace path. Defaults to current directory.
        iFlag_standalone_in (int, optional): Flag for standalone mode. Defaults to 1.
        iFlag_use_mesh_dem_in (int, optional): Flag to use mesh DEM. Defaults to 0.
        iCase_index_in (int, optional): Case index. Defaults to 1.
        dResolution_degree_in (float, optional): Resolution in degrees. Defaults to 0.5.
        dResolution_meter_in (float, optional): Resolution in meters. Defaults to 50000.
        sDate_in (str, optional): Date string. Defaults to '20220202'.
        sMesh_type_in (str, optional): Mesh type. Defaults to 'hexagon'.
        sModel_in (str, optional): Model name. Defaults to 'pyflowline'.
        iFlag_save_mesh_in (int, optional): Flag to save mesh. Defaults to 1.
        iFlag_simplification_in (int, optional): Flag for simplification. Defaults to 1.
        iFlag_create_mesh_in (int, optional): Flag to create mesh. Defaults to 1.
        iFlag_intersect_in (int, optional): Flag for intersection. Defaults to 1.
        iFlag_resample_method_in (int, optional): Resampling method flag. Defaults to 1.
        iFlag_global_in (int, optional): Global flag. Defaults to 0.
        iFlag_multiple_outlet_in (int, optional): Multiple outlet flag. Defaults to 0.
        iFlag_elevation_profile_in (int, optional): Elevation profile flag. Defaults to 1.
        iFlag_rotation_in (int, optional): Rotation flag. Defaults to 0.
        iFlag_stream_burning_topology_in (int, optional): Stream burning topology flag. Defaults to 1.
        iFlag_save_elevation_in (int, optional): Save elevation flag. Defaults to 1.
        nOutlet_in (int, optional): Number of outlets. Defaults to 1.
        dLongitude_left_in (float, optional): Left longitude boundary. Defaults to -180.
        dLongitude_right_in (float, optional): Right longitude boundary. Defaults to 180.
        dLatitude_bot_in (float, optional): Bottom latitude boundary. Defaults to -90.
        dLatitude_top_in (float, optional): Top latitude boundary. Defaults to 90.
        sRegion_in (str, optional): Region name. Defaults to 'susquehanna'.
        sJob_in (str, optional): Job name. Defaults to 'pyflowline'.
        flowline_info_in (str, optional): Flowline info filename. Defaults to 'flowline_info.json'.
        sFilename_mesh_info_in (str, optional): Mesh info filename. Defaults to 'mesh_info.json'.
        sFilename_elevation_in (str, optional): Elevation filename. Defaults to 'elevation.json'.

    Returns:
        flowlinecase: A configured flowlinecase object
    """

    # Ensure input pathnames are strings
    sFilename_json = str(sFilename_json)
    #use the current working directory as the input and output workspace
    if sWorkspace_input is None:
        sWorkspace_input = str(Path.cwd())
    else:
        sWorkspace_input = str(sWorkspace_input)

    if sWorkspace_output is None:
        sWorkspace_output = str(Path.cwd())
    else:
        sWorkspace_output = str(sWorkspace_output)

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

    if sModel_in is not None:
        sModel = sModel_in
    else:
        sModel = 'pyflowline'
        pass

    if dResolution_degree_in is not None:
        dResolution_degree = dResolution_degree_in
    else:
        dResolution_degree = 0.5
        pass

    if dResolution_meter_in is not None:
        dResolution_meter = dResolution_meter_in
    else:
        dResolution_meter = 50000
        pass

    nBasin = 1 if nOutlet_in is None else nOutlet_in

    # Use a dict to initialize the class
    aConfig = {}

    # Set core configuration values
    aConfig['iFlag_use_mesh_dem'] = iFlag_use_mesh_dem
    aConfig['iFlag_save_mesh'] = 1 if iFlag_save_mesh_in is None else iFlag_save_mesh_in
    aConfig['iFlag_simplification'] = 1 if iFlag_simplification_in is None else iFlag_simplification_in
    aConfig['iFlag_create_mesh'] = 1 if iFlag_create_mesh_in is None else iFlag_create_mesh_in
    aConfig['iFlag_intersect'] = 1 if iFlag_intersect_in is None else iFlag_intersect_in
    aConfig['iFlag_resample_method'] = 1 if iFlag_resample_method_in is None else iFlag_resample_method_in
    aConfig['iFlag_global'] = 0 if iFlag_global_in is None else iFlag_global_in
    aConfig['iFlag_multiple_outlet'] = 0 if iFlag_multiple_outlet_in is None else iFlag_multiple_outlet_in
    aConfig['iFlag_elevation_profile'] = 1 if iFlag_elevation_profile_in is None else iFlag_elevation_profile_in
    aConfig['iFlag_rotation'] = 0 if iFlag_rotation_in is None else iFlag_rotation_in
    aConfig['iFlag_stream_burning_topology'] = 1 if iFlag_stream_burning_topology_in is None else iFlag_stream_burning_topology_in
    aConfig['iFlag_save_elevation'] = 1 if iFlag_save_elevation_in is None else iFlag_save_elevation_in
    aConfig['nOutlet'] = nBasin
    aConfig['dResolution_degree'] = dResolution_degree
    aConfig['dResolution_meter'] = dResolution_meter
    aConfig['dLongitude_left'] = -180 if dLongitude_left_in is None else dLongitude_left_in
    aConfig['dLongitude_right'] = 180 if dLongitude_right_in is None else dLongitude_right_in
    aConfig['dLatitude_bot'] = -90 if dLatitude_bot_in is None else dLatitude_bot_in
    aConfig['dLatitude_top'] = 90 if dLatitude_top_in is None else dLatitude_top_in
    aConfig['sFilename_model_configuration'] = sFilename_json
    aConfig['sWorkspace_input'] = sWorkspace_input
    aConfig['sWorkspace_output'] = sWorkspace_output
    aConfig['sRegion'] = 'susquehanna' if sRegion_in is None else sRegion_in
    aConfig['sModel'] = sModel
    aConfig['iCase_index'] = iCase_index
    aConfig['sMesh_type'] = sMesh_type
    aConfig['sJob'] = 'pyflowline' if sJob_in is None else sJob_in
    aConfig['sDate'] = sDate

    # File path configurations
    aConfig['flowline_info'] = 'flowline_info.json' if flowline_info_in is None else flowline_info_in
    aConfig['sFilename_mesh_info'] = 'mesh_info.json' if sFilename_mesh_info_in is None else sFilename_mesh_info_in
    aConfig['sFilename_elevation'] = 'elevation.json' if sFilename_elevation_in is None else sFilename_elevation_in

    oModel = flowlinecase(aConfig)

    # Generate basin
    sDirname = os.path.dirname(sFilename_json)
    sFilename = Path(sFilename_json).stem + '_basins.json'
    sFilename_basins_json = os.path.join(sDirname, sFilename)

    aBasin = create_template_basin_configuration_file(
        sFilename_basins_json,
        nBasin,
        sWorkspace_input,
        oModel.sWorkspace_output)

    oModel.aBasin = aBasin
    oModel.sFilename_basins = sFilename_basins_json
    oModel.export_config_to_json(sFilename_json)

    return oModel
