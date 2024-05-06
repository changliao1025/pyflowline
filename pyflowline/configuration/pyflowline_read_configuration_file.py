import os
from pathlib import Path
import datetime
import json
from pyflowline.classes.pycase import flowlinecase
pDate = datetime.datetime.today()
sDate_default = "{:04d}".format(
    pDate.year) + "{:02d}".format(pDate.month) + "{:02d}".format(pDate.day)


def pyflowline_read_configuration_file(sFilename_configuration_in,
                                             iFlag_standalone_in=None,
                                             iFlag_use_mesh_dem_in=None,
                                             iCase_index_in=None,
                                             iResolution_index_in = None,
                                             dResolution_degree_in=None,
                                             dResolution_meter_in=None,
                                             sMesh_type_in=None,
                                             sModel_in=None,
                                             sDate_in=None,
                                              sDggrid_type_in = None,
                                             sWorkspace_output_in=None):
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
        iCase_index = int(aConfig['iCase_index'])

    if iResolution_index_in is not None:
        iResolution_index = iResolution_index_in
    else:
        if "iResolution_index" in aConfig:
            iResolution_index =  int( aConfig['iResolution_index'])
        else:
            iResolution_index = 10

        pass

    if iFlag_standalone_in is not None:
        iFlag_standalone = iFlag_standalone_in
    else:
        iFlag_standalone = int(aConfig['iFlag_standalone'])

    if iFlag_use_mesh_dem_in is not None:
        iFlag_use_mesh_dem = iFlag_use_mesh_dem_in
    else:
        iFlag_use_mesh_dem = int(aConfig['iFlag_use_mesh_dem'])

    if sMesh_type_in is not None:
        sMesh_type = sMesh_type_in
    else:
        sMesh_type = aConfig['sMesh_type']
        pass

    if sDggrid_type_in is not None:
        sDggrid_type = sDggrid_type_in
    else:
        if "sDggrid_type" in aConfig:
            sDggrid_type = aConfig["sDggrid_type"]
        else:
            sDggrid_type = 'ISEA3H'

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
        # try to create this output folder first using

    try:
        print(sWorkspace_output)
        Path(sWorkspace_output).mkdir(parents=True, exist_ok=True)
    except ValueError:
        print("The specified output workspace cannot be created!")
        exit

    aConfig['iCase_index'] = iCase_index
    aConfig['iFlag_standalone'] = iFlag_standalone
    aConfig['iFlag_use_mesh_dem'] = iFlag_use_mesh_dem
    aConfig["iResolution_index"] = iResolution_index
    aConfig['dResolution_degree'] = dResolution_degree
    aConfig['dResolution_meter'] = dResolution_meter

    aConfig['sDate'] = sDate
    aConfig['sModel'] = sModel
    aConfig['sMesh_type'] = sMesh_type
    aConfig["sDggrid_type"] = sDggrid_type
    aConfig['sWorkspace_output'] = sWorkspace_output

    aConfig["sFilename_model_configuration"] = sFilename_configuration_in

    # based on global variable, a few variables are calculate once
    # calculate the modflow simulation period
    # https://docs.python.org/3/library/datetime.html#datetime-objects

    # aConfig

    # simulation

    oPyflowline = flowlinecase(aConfig)

    return oPyflowline
