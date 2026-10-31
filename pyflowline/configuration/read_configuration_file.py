from pathlib import Path

from pyflowline.classes.pycase import flowlinecase
from pyflowline.configuration.config_manager import _resolve_config_path, load_pyflowline_config


def pyflowline_read_configuration_file(
    sFilename_configuration_in,
    iFlag_standalone_in=1,
    iFlag_create_directory_in=None,
    iFlag_use_mesh_dem_in=None,
    iCase_index_in=None,
    iResolution_index_in=None,
    dResolution_degree_in=None,
    dResolution_meter_in=None,
    sMesh_type_in=None,
    sModel_in="pyflowline",
    sDate_in=None,
    sDggrid_type_in=None,
    sWorkspace_output_in=None,
):
    """Read a PyFlowline model configuration from YAML or JSON.

    Supports:
        - YAML (.yaml / .yml) with comments (#)
        - JSON (.json) for backward compatibility
        - Extensionless file paths and directory paths
        - Two-tier sectioned specification or legacy flat dictionary
        - Runtime parameter overrides

    Args:
        sFilename_configuration_in (str or Path): Path to configuration file or directory
        iFlag_standalone_in (int, optional): Standalone mode flag. Defaults to 1.
        iFlag_create_directory_in (int, optional): Directory creation flag.
        iFlag_use_mesh_dem_in (int, optional): Use mesh DEM flag.
        iCase_index_in (int, optional): Case index override.
        iResolution_index_in (int, optional): Resolution index override.
        dResolution_degree_in (float, optional): Resolution in degrees override.
        dResolution_meter_in (float, optional): Resolution in meters override.
        sMesh_type_in (str, optional): Mesh type override.
        sModel_in (str, optional): Model name override.
        sDate_in (str, optional): Simulation date override.
        sDggrid_type_in (str, optional): DGGRID type override.
        sWorkspace_output_in (str or Path, optional): Output workspace override.

    Returns:
        flowlinecase: Configured flowlinecase object ready for simulation.
    """
    resolved_path = _resolve_config_path(sFilename_configuration_in)

    if not resolved_path.is_file():
        print(f"{sFilename_configuration_in} does not exist")
        return None

    # Load configuration (handles YAML comments, relative paths, platform paths)
    aConfig = load_pyflowline_config(resolved_path)

    # Apply runtime argument overrides
    if iCase_index_in is not None:
        aConfig["iCase_index"] = int(iCase_index_in)

    if iResolution_index_in is not None:
        aConfig["iResolution_index"] = int(iResolution_index_in)

    if iFlag_standalone_in is not None:
        aConfig["iFlag_standalone"] = int(iFlag_standalone_in)

    if iFlag_use_mesh_dem_in is not None:
        aConfig["iFlag_use_mesh_dem"] = int(iFlag_use_mesh_dem_in)

    if sMesh_type_in is not None:
        aConfig["sMesh_type"] = str(sMesh_type_in).lower()

    if sDggrid_type_in is not None:
        aConfig["sDggrid_type"] = str(sDggrid_type_in)

    if sModel_in is not None:
        aConfig["sModel"] = str(sModel_in)

    if sDate_in is not None:
        aConfig["sDate"] = str(sDate_in)

    if dResolution_degree_in is not None:
        aConfig["dResolution_degree"] = float(dResolution_degree_in)

    if dResolution_meter_in is not None:
        aConfig["dResolution_meter"] = float(dResolution_meter_in)

    if sWorkspace_output_in is not None:
        aConfig["sWorkspace_output"] = str(sWorkspace_output_in)

    sWorkspace_output = aConfig.get("sWorkspace_output", "./output")
    try:
        Path(sWorkspace_output).mkdir(parents=True, exist_ok=True)
    except Exception as exc:
        print(f"Warning: Could not create output workspace {sWorkspace_output}: {exc}")

    aConfig["sFilename_model_configuration"] = str(resolved_path)

    oPyflowline = flowlinecase(aConfig, iFlag_create_directory_in=iFlag_create_directory_in)

    return oPyflowline
