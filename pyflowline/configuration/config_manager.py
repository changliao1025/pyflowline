from __future__ import annotations

import datetime
import importlib.resources as importlib_resources
import json
import os
import platform
from pathlib import Path
from typing import Any

try:
    import yaml

    HAVE_YAML = True
except ImportError:
    HAVE_YAML = False

# ==============================================================================
# Project & Package Path Resolution (formerly path_manager)
# ==============================================================================


def pyflowline_package_root(package_name: str = "pyflowline") -> Path:
    """Get the installation root directory of the package."""
    try:
        return Path(str(importlib_resources.files(package_name))).resolve()
    except Exception:
        return Path(__file__).resolve().parents[1]


def root_path_from_setup_file() -> Path:
    """Return the top-level (root) path of the pyflowline project.

    Navigates upwards from this file's path until finding a directory
    containing setup.py or pyproject.toml.
    """
    this_path = Path(__file__).resolve()
    for parent in this_path.parents:
        if (parent / "setup.py").exists() or (parent / "pyproject.toml").exists():
            return parent
    raise FileNotFoundError("Project root marker (setup.py or pyproject.toml) not found.")


def root_path_from_pyflowline_package_root() -> Path:
    """Get the parent directory of the pyflowline package root."""
    return pyflowline_package_root().parent


def root_path_from_path_manager_location() -> Path:
    """Return the top-level root path assuming fixed directory layout:
    root_path/pyflowline/configuration/config_manager.py
    """
    return Path(__file__).resolve().parents[2]


def pyflowline_project_root() -> Path:
    """Find the root path of the project.

    Searches upward for a project marker (setup.py, pyproject.toml),
    falling back to package root parent.
    """
    try:
        return root_path_from_setup_file()
    except FileNotFoundError:
        return root_path_from_pyflowline_package_root()
    except Exception as exc:
        raise RuntimeError(
            f"An unexpected error occurred while locating the project root: {exc}"
        ) from exc


def join_project_path(*args: str) -> Path:
    """Join components of a path relative to the project root."""
    return pyflowline_project_root().joinpath(*args)


class _PathManagerNamespace:
    """Backwards-compatible path_manager namespace object."""

    pyflowline_project_root = staticmethod(pyflowline_project_root)
    pyflowline_package_root = staticmethod(pyflowline_package_root)
    join_project_path = staticmethod(join_project_path)
    root_path_from_setup_file = staticmethod(root_path_from_setup_file)
    root_path_from_pyflowline_package_root = staticmethod(root_path_from_pyflowline_package_root)
    root_path_from_path_manager_location = staticmethod(root_path_from_path_manager_location)


path_manager = _PathManagerNamespace()

# ==============================================================================
# PyFlowline Workflow Configuration Manager
#
# Primary specification format is YAML (.yaml / .yml) to allow inline and block
# comments (#), making workflow specifications transparent, readable, and
# user-friendly for both case and basin configurations.
#
# JSON (.json) is fully supported for programmatic workflows and backward
# compatibility.
# ==============================================================================

DEFAULT_BASIN_CONFIG: dict[str, Any] = {
    "lBasinID": 1,
    "dLatitude_outlet_degree": -180.0,
    "dLongitude_outlet_degree": 180.0,
    "dAccumulation_threshold": -90.0,
    "dThreshold_small_river": 90.0,
    "iFlag_dam": 0,
    "iFlag_debug": 0,
    "iFlag_disconnected": 0,
    "sFilename_dam": "ICoM_dams.csv",
    "sFilename_flowline_filter": "streamord7above.shp",
    "sFilename_flowline_raw": "allflowline.shp",
    "sFilename_flowline_topo": "flowline.csv",
    "sWorkspace_output_basin": "",
}

DEFAULT_FLOWLINE_CONFIG: dict[str, Any] = {
    "iFlag_use_mesh_dem": 0,
    "iFlag_save_mesh": 1,
    "iFlag_simplification": 1,
    "iFlag_create_mesh": 1,
    "iFlag_intersect": 1,
    "iFlag_resample_method": 1,
    "iFlag_global": 0,
    "iFlag_multiple_outlet": 0,
    "iFlag_elevation_profile": 1,
    "iFlag_rotation": 0,
    "iFlag_stream_burning_topology": 1,
    "iFlag_save_elevation": 1,
    "iFlag_standalone": 1,
    "iFlag_flowline": 1,
    "iFlag_run_jigsaw": 0,
    "iFlag_run_dggrid": 0,
    "iFlag_fill_hole": 0,
    "iFlag_break_by_distance": 0,
    "iFlag_user_provided_binary": 0,
    "iResolution_index": 10,
    "sDggrid_type": "ISEA3H",
    "nOutlet": 1,
    "dResolution_degree": 0.5,
    "dResolution_meter": 50000.0,
    "dLongitude_left": -180.0,
    "dLongitude_right": 180.0,
    "dLatitude_bot": -90.0,
    "dLatitude_top": 90.0,
    "dLongitude_outlet_degree": -180.0,
    "dLatitude_outlet_degree": 90.0,
    "dAccumulation_threshold": 100000.0,
    "dThreshold_small_river": 10000.0,
    "dThreshold_break_by_distance": 5000.0,
    "sRegion": "susquehanna",
    "sModel": "pyflowline",
    "iCase_index": 1,
    "sMesh_type": "hexagon",
    "sJob": "pyflowline",
    "sDate": "20220202",
    "flowline_info": "flowline_info.json",
    "sFilename_mesh_info": "mesh_info.json",
    "sFilename_elevation": "elevation.json",
    "sWorkspace_input": "",
    "sWorkspace_output": "./output",
    "sWorkspace_data": "",
    "sWorkspace_project": "",
    "sWorkspace_bin": "",
    "sFilename_model_configuration": "",
    "sFilename_spatial_reference": "",
    "sFilename_dem": "",
    "sFilename_mesh_netcdf": "",
    "sFilename_mpas_mesh_netcdf": "",
    "sFilename_jigsaw_mesh_netcdf": "",
    "sFilename_mesh_boundary": "",
    "sFilename_coastline_boundary": "",
    "sFilename_basins": "",
    "sFilename_jigsaw_configuration": "",
    "sFilename_dggrid": "",
}


def _select_platform_path(value: Any) -> str:
    """Resolve platform-specific dictionary or string path."""
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, str):
        return value
    if not isinstance(value, dict):
        raise ValueError("Path values must be strings or platform-specific mappings")

    current_platform = platform.system()
    selected = value.get(current_platform, value.get("default", ""))
    if not isinstance(selected, str):
        raise ValueError(f"Path for {current_platform} must be a string")
    return selected


def _resolve_config_path(
    filename: str | Path,
    candidates: tuple[str, ...] = (
        "pyflowline_config.yaml",
        "pyflowline_config.yml",
        "workflow_config.yaml",
        "workflow_config.yml",
        "configuration.yaml",
        "configuration.yml",
        "pyflowline_config.json",
        "workflow_config.json",
        "configuration.json",
    ),
) -> Path:
    """Resolve configuration file path handling directories and missing extensions."""
    path = Path(filename).expanduser().resolve()
    if path.is_dir():
        for candidate_name in candidates:
            candidate = path / candidate_name
            if candidate.is_file():
                return candidate
        raise FileNotFoundError(
            f"No configuration file found in directory: {path}. Candidates checked: {candidates}"
        )
    if not path.is_file():
        for ext in (".yaml", ".yml", ".json"):
            candidate = path.with_suffix(ext)
            if candidate.is_file():
                return candidate
    return path


def _parse_config_file(path: Path) -> Any:
    """Parse a YAML or JSON configuration file, seamlessly allowing comments."""
    ext = path.suffix.lower()
    if ext in (".yaml", ".yml"):
        if not HAVE_YAML:
            raise ImportError(
                "PyYAML is required to parse YAML configuration files. "
                "Please install pyyaml with `pip install pyyaml` or `conda install pyyaml`."
            )
        with path.open(encoding="utf-8") as f:
            raw = yaml.safe_load(f)
    elif ext == ".json":
        with path.open(encoding="utf-8") as f:
            raw = json.load(f)
    else:
        # Fallback: attempt YAML first (YAML accepts comments and JSON syntax), then JSON
        try:
            if HAVE_YAML:
                with path.open(encoding="utf-8") as f:
                    raw = yaml.safe_load(f)
            else:
                with path.open(encoding="utf-8") as f:
                    raw = json.load(f)
        except Exception:
            with path.open(encoding="utf-8") as f:
                raw = json.load(f)

    if raw is None:
        raw = {}
    return raw


def _dump_config_file(data: Any, path: Path, fmt: str | None = None) -> None:
    """Dump data to a YAML or JSON file."""
    path = Path(path).expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)

    if fmt is None:
        ext = path.suffix.lower()
        if ext in (".yaml", ".yml"):
            fmt = "yaml"
        elif ext == ".json":
            fmt = "json"
        else:
            fmt = "yaml"
    else:
        fmt = fmt.lower()

    # Recursively convert Path objects to strings and exclude non-serializable fields
    def _clean_for_export(obj: Any) -> Any:
        if isinstance(obj, Path):
            return str(obj)
        if isinstance(obj, dict):
            return {
                k: _clean_for_export(v)
                for k, v in obj.items()
                if k not in ("config_file", "pRTree_mesh", "pRTree_flowline", "pRTree_edge")
            }
        if isinstance(obj, list):
            return [_clean_for_export(v) for v in obj]
        return obj

    cleaned = _clean_for_export(data)

    if fmt in ("yaml", "yml"):
        if not HAVE_YAML:
            raise ImportError(
                "PyYAML is required to write YAML configuration files. "
                "Please install pyyaml with `pip install pyyaml`."
            )
        with path.open("w", encoding="utf-8") as f:
            yaml.safe_dump(cleaned, f, sort_keys=False, default_flow_style=False)
    elif fmt == "json":
        with path.open("w", encoding="utf-8") as f:
            json.dump(cleaned, f, indent=4)
    else:
        raise ValueError(f"Unsupported configuration format {fmt!r}. Expected 'yaml' or 'json'.")


def _normalize_basin_dict(basin: dict[str, Any]) -> dict[str, Any]:
    """Normalize a basin dictionary from either sectioned YAML or flat keys."""
    res = DEFAULT_BASIN_CONFIG.copy()

    # Direct flat keys if provided
    for key, val in basin.items():
        if key in res:
            res[key] = val

    # Sectioned / friendly aliases
    if "id" in basin:
        res["lBasinID"] = int(basin["id"])
    if "basin_id" in basin:
        res["lBasinID"] = int(basin["basin_id"])

    outlet = basin.get("outlet", {})
    if isinstance(outlet, dict):
        if "longitude" in outlet:
            res["dLongitude_outlet_degree"] = float(outlet["longitude"])
        if "latitude" in outlet:
            res["dLatitude_outlet_degree"] = float(outlet["latitude"])

    if "accumulation_threshold" in basin:
        res["dAccumulation_threshold"] = float(basin["accumulation_threshold"])
    if "small_river_threshold" in basin:
        res["dThreshold_small_river"] = float(basin["small_river_threshold"])
    if "small_river" in basin:
        res["dThreshold_small_river"] = float(basin["small_river"])

    flags = basin.get("flags", {})
    if isinstance(flags, dict):
        if "dam" in flags:
            res["iFlag_dam"] = int(flags["dam"])
        if "debug" in flags:
            res["iFlag_debug"] = int(flags["debug"])
        if "disconnected" in flags:
            res["iFlag_disconnected"] = int(flags["disconnected"])

    paths = basin.get("paths", {})
    if isinstance(paths, dict):
        if "dam" in paths:
            res["sFilename_dam"] = _select_platform_path(paths["dam"])
        if "flowline_filter" in paths:
            res["sFilename_flowline_filter"] = _select_platform_path(paths["flowline_filter"])
        if "flowline_raw" in paths:
            res["sFilename_flowline_raw"] = _select_platform_path(paths["flowline_raw"])
        if "flowline_topo" in paths:
            res["sFilename_flowline_topo"] = _select_platform_path(paths["flowline_topo"])
        if "output_basin" in paths:
            res["sWorkspace_output_basin"] = _select_platform_path(paths["output_basin"])

    return res


class WorkflowSpec:
    """Mesh-agnostic front-end specification for a pyflowline run.

    Tier 1 of the Two-Tier Compiler Pattern: captures user intent in a clean,
    human-readable structure without requiring backend-internal parameter names.
    """

    def __init__(
        self,
        mesh_type: str = "hexagon",
        case_index: int = 1,
        date: str | None = None,
        model: str = "pyflowline",
        job: str | None = None,
        region: dict[str, Any] | None = None,
        resolutions: dict[str, Any] | None = None,
        paths: dict[str, Any] | None = None,
        flags: dict[str, Any] | None = None,
        basins: list[dict[str, Any]] | None = None,
        jigsaw: dict[str, Any] | None = None,
    ):
        self.mesh_type = str(mesh_type).lower()
        self.case_index = int(case_index)
        self.date = date
        self.model = model
        self.job = job or self.mesh_type
        self.region = dict(region or {})
        self.resolutions = dict(resolutions or {})
        self.paths = dict(paths or {})
        self.flags = dict(flags or {})
        self.basins = list(basins or [])
        self.jigsaw = dict(jigsaw or {})

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> WorkflowSpec:
        workflow = data.get("workflow", {}) or {}
        return cls(
            mesh_type=workflow.get("mesh_type", data.get("sMesh_type", "hexagon")),
            case_index=workflow.get("case_index", data.get("iCase_index", 1)),
            date=workflow.get("date", data.get("sDate")),
            model=workflow.get("model", data.get("sModel", "pyflowline")),
            job=workflow.get("job", data.get("sJob")),
            region=data.get("region", {}),
            resolutions=data.get("resolutions", {}),
            paths=data.get("paths", {}),
            flags=data.get("flags", {}),
            basins=data.get("basins", []),
            jigsaw=data.get("jigsaw", {}),
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "workflow": {
                "mesh_type": self.mesh_type,
                "case_index": self.case_index,
                "date": self.date,
                "model": self.model,
                "job": self.job,
            },
            "resolutions": self.resolutions,
            "region": self.region,
            "paths": self.paths,
            "flags": self.flags,
            "basins": self.basins,
            "jigsaw": self.jigsaw,
        }


def compile_workflow(spec_or_config: dict[str, Any] | WorkflowSpec) -> dict[str, Any]:
    """Compile Tier 1 workflow specification into backend-native pyflowline configuration.

    Supports both sectioned Tier 1 configurations (with 'workflow', 'resolutions',
    'region', 'paths', 'flags', 'basins', 'jigsaw' sections) and flat native dictionaries.
    """
    if isinstance(spec_or_config, WorkflowSpec):
        raw = spec_or_config.to_dict()
    elif isinstance(spec_or_config, dict):
        raw = spec_or_config
    else:
        raise TypeError("spec_or_config must be a dictionary or WorkflowSpec")

    # Start with default configuration
    native = DEFAULT_FLOWLINE_CONFIG.copy()

    # Pass through any existing top-level native keys (for backward compatibility)
    for k, v in raw.items():
        if k in native:
            native[k] = v

    # Section 1: workflow
    workflow = raw.get("workflow", {})
    if isinstance(workflow, dict) and workflow:
        if "mesh_type" in workflow:
            native["sMesh_type"] = str(workflow["mesh_type"]).lower()
        if "case_index" in workflow:
            native["iCase_index"] = int(workflow["case_index"])
        if "date" in workflow and workflow["date"]:
            native["sDate"] = str(workflow["date"])
        if "model" in workflow:
            native["sModel"] = str(workflow["model"])
        if "job" in workflow:
            native["sJob"] = str(workflow["job"])
        elif "sJob" not in native or not native["sJob"]:
            native["sJob"] = native["sMesh_type"]
        if "region" in workflow:
            native["sRegion"] = str(workflow["region"])
        if "standalone" in workflow:
            native["iFlag_standalone"] = 1 if workflow["standalone"] else 0
        if "dggrid_type" in workflow:
            native["sDggrid_type"] = str(workflow["dggrid_type"])
        if "outlets" in workflow:
            native["nOutlet"] = int(workflow["outlets"])
        elif "nOutlet" in workflow:
            native["nOutlet"] = int(workflow["nOutlet"])

    # Section 2: resolutions
    resolutions = raw.get("resolutions", {})
    if isinstance(resolutions, dict) and resolutions:
        if "degree" in resolutions and resolutions["degree"] is not None:
            native["dResolution_degree"] = float(resolutions["degree"])
        if "meter" in resolutions and resolutions["meter"] is not None:
            native["dResolution_meter"] = float(resolutions["meter"])
        if "resolution_index" in resolutions and resolutions["resolution_index"] is not None:
            native["iResolution_index"] = int(resolutions["resolution_index"])
        if (
            "accumulation_threshold" in resolutions
            and resolutions["accumulation_threshold"] is not None
        ):
            native["dAccumulation_threshold"] = float(resolutions["accumulation_threshold"])
        if (
            "small_river_threshold" in resolutions
            and resolutions["small_river_threshold"] is not None
        ):
            native["dThreshold_small_river"] = float(resolutions["small_river_threshold"])
        elif "small_river" in resolutions and resolutions["small_river"] is not None:
            native["dThreshold_small_river"] = float(resolutions["small_river"])
        if "break_by_distance" in resolutions and resolutions["break_by_distance"] is not None:
            native["dThreshold_break_by_distance"] = float(resolutions["break_by_distance"])

    # Section 3: region
    region = raw.get("region", {})
    if isinstance(region, dict) and region:
        bbox = region.get("bbox", {})
        if isinstance(bbox, dict):
            if "left" in bbox:
                native["dLongitude_left"] = float(bbox["left"])
            if "right" in bbox:
                native["dLongitude_right"] = float(bbox["right"])
            if "bottom" in bbox:
                native["dLatitude_bot"] = float(bbox["bottom"])
            if "top" in bbox:
                native["dLatitude_top"] = float(bbox["top"])

        # Also support direct coordinates in region
        for k_src, k_dst in (
            ("longitude_left", "dLongitude_left"),
            ("longitude_right", "dLongitude_right"),
            ("latitude_bot", "dLatitude_bot"),
            ("latitude_bottom", "dLatitude_bot"),
            ("latitude_top", "dLatitude_top"),
        ):
            if k_src in region and region[k_src] is not None:
                native[k_dst] = float(region[k_src])

        outlet = region.get("outlet", {})
        if isinstance(outlet, dict):
            if "longitude" in outlet:
                native["dLongitude_outlet_degree"] = float(outlet["longitude"])
            if "latitude" in outlet:
                native["dLatitude_outlet_degree"] = float(outlet["latitude"])
        elif "dLongitude_outlet_degree" in region:
            native["dLongitude_outlet_degree"] = float(region["dLongitude_outlet_degree"])
            native["dLatitude_outlet_degree"] = float(region["dLatitude_outlet_degree"])

        if "cellid_outlet" in region:
            native["lCellID_outlet"] = int(region["cellid_outlet"])

    # Section 4: flags
    flags = raw.get("flags", {})
    if isinstance(flags, dict) and flags:
        flag_map = {
            "create_mesh": "iFlag_create_mesh",
            "save_mesh": "iFlag_save_mesh",
            "simplification": "iFlag_simplification",
            "intersect": "iFlag_intersect",
            "flowline": "iFlag_flowline",
            "use_mesh_dem": "iFlag_use_mesh_dem",
            "rotation": "iFlag_rotation",
            "stream_burning_topology": "iFlag_stream_burning_topology",
            "resample_method": "iFlag_resample_method",
            "global": "iFlag_global",
            "multiple_outlet": "iFlag_multiple_outlet",
            "run_jigsaw": "iFlag_run_jigsaw",
            "run_dggrid": "iFlag_run_dggrid",
            "fill_hole": "iFlag_fill_hole",
            "break_by_distance": "iFlag_break_by_distance",
            "user_provided_binary": "iFlag_user_provided_binary",
            "elevation_profile": "iFlag_elevation_profile",
            "save_elevation": "iFlag_save_elevation",
            "disconnected": "iFlag_disconnected",
            "land_ocean_mask": "iFlag_land_ocean_mask",
            "standalone": "iFlag_standalone",
            "mesh_boundary": "iFlag_mesh_boundary",
            "simulation": "iFlag_simulation",
            "visualization": "iFlag_visualization",
        }
        for flag_k, native_k in flag_map.items():
            if flag_k in flags and flags[flag_k] is not None:
                native[native_k] = int(flags[flag_k])

    # Section 5: paths
    paths = raw.get("paths", {})
    if isinstance(paths, dict) and paths:
        path_map = {
            "workspace_input": "sWorkspace_input",
            "workspace_output": "sWorkspace_output",
            "workspace_data": "sWorkspace_data",
            "workspace_project": "sWorkspace_project",
            "workspace_bin": "sWorkspace_bin",
            "mesh_boundary": "sFilename_mesh_boundary",
            "mesh_netcdf": "sFilename_mesh_netcdf",
            "mpas_mesh_netcdf": "sFilename_mpas_mesh_netcdf",
            "jigsaw_mesh_netcdf": "sFilename_jigsaw_mesh_netcdf",
            "coastline_boundary": "sFilename_coastline_boundary",
            "spatial_reference": "sFilename_spatial_reference",
            "dem": "sFilename_dem",
            "basins": "sFilename_basins",
            "jigsaw_configuration": "sFilename_jigsaw_configuration",
            "dggrid": "sFilename_dggrid",
            "flowline_info": "flowline_info",
            "mesh_info": "sFilename_mesh_info",
            "elevation": "sFilename_elevation",
            "model_configuration": "sFilename_model_configuration",
        }
        for path_k, native_k in path_map.items():
            if path_k in paths and paths[path_k] is not None:
                native[native_k] = _select_platform_path(paths[path_k])

        # Also support direct sFilename_ / sWorkspace_ in paths
        for k, v in paths.items():
            if k.startswith("sFilename_") or k.startswith("sWorkspace_"):
                native[k] = _select_platform_path(v)

    # Section 6: basins (inline list or single dict of basin specifications)
    basins = raw.get("basins") or raw.get("aBasin")
    if basins:
        if isinstance(basins, dict):
            basins = [basins]
        if isinstance(basins, list):
            compiled_basins = [_normalize_basin_dict(b) for b in basins if isinstance(b, dict)]
            native["aBasin"] = compiled_basins
            native["nOutlet"] = len(compiled_basins)
            if len(compiled_basins) > 1:
                native["iFlag_multiple_outlet"] = 1

    # Section 7: jigsaw
    jigsaw = raw.get("jigsaw")
    if isinstance(jigsaw, dict) and jigsaw:
        native["jigsaw"] = jigsaw
        native["iFlag_run_jigsaw"] = 1

    # Date normalization
    if not native.get("sDate"):
        pDate = datetime.datetime.today()
        native["sDate"] = f"{pDate.year:04d}" + f"{pDate.month:02d}" + f"{pDate.day:02d}"

    # Carry forward config_file reference if present
    if "config_file" in raw:
        native["config_file"] = raw["config_file"]
        native["sFilename_model_configuration"] = str(raw["config_file"])

    return native


def load_pyflowline_config(filename: str | Path) -> dict[str, Any]:
    """Load a workflow YAML or JSON file, apply defaults, and resolve paths.

    Features:
        - Primary format: YAML (.yaml, .yml) with inline/block comments (#)
        - Full backward compatibility with JSON (.json)
        - Directory discovery (checks for default configuration names)
        - Extensionless stem resolution (.yaml -> .yml -> .json)
        - Relative path resolution against config file location
        - Platform-specific path mapping resolution (Linux, Darwin, Windows, default)
    """
    resolved_path = _resolve_config_path(filename)
    raw_dict = _parse_config_file(resolved_path)

    if not isinstance(raw_dict, dict):
        raise ValueError(f"Configuration in {resolved_path} must be a dictionary / mapping object")

    raw_dict["config_file"] = resolved_path

    # Compile into native dictionary
    config = compile_workflow(raw_dict)

    # Resolve relative paths relative to the configuration file directory
    base_dir = resolved_path.parent
    path_keys = [
        "sWorkspace_input",
        "sWorkspace_output",
        "sWorkspace_data",
        "sWorkspace_project",
        "sWorkspace_bin",
        "sFilename_spatial_reference",
        "sFilename_dem",
        "sFilename_mesh_netcdf",
        "sFilename_mpas_mesh_netcdf",
        "sFilename_jigsaw_mesh_netcdf",
        "sFilename_mesh_boundary",
        "sFilename_coastline_boundary",
        "sFilename_basins",
        "sFilename_jigsaw_configuration",
        "sFilename_dggrid",
    ]

    for key in path_keys:
        val = config.get(key)
        if val and isinstance(val, str) and not os.path.isabs(val):
            # Only resolve if the path doesn't already start with an absolute root
            resolved_file = (base_dir / val).resolve()
            config[key] = str(resolved_file)

    # If inline basins are present, resolve relative paths inside each basin
    if "aBasin" in config and isinstance(config["aBasin"], list):
        basin_path_keys = [
            "sFilename_dam",
            "sFilename_flowline_filter",
            "sFilename_flowline_raw",
            "sFilename_flowline_topo",
            "sWorkspace_output_basin",
        ]
        for b in config["aBasin"]:
            if isinstance(b, dict):
                for bk in basin_path_keys:
                    b_val = b.get(bk)
                    if b_val and isinstance(b_val, str) and not os.path.isabs(b_val):
                        b[bk] = str((base_dir / b_val).resolve())

    # If external basins file is referenced and nOutlet wasn't explicitly configured, inspect basin count
    workflow_raw = (
        raw_dict.get("workflow", {}) if isinstance(raw_dict.get("workflow"), dict) else {}
    )
    if (
        config.get("sFilename_basins")
        and not workflow_raw.get("outlets")
        and not workflow_raw.get("nOutlet")
    ):
        try:
            basin_file = Path(config["sFilename_basins"])
            if basin_file.is_file():
                basin_records = BasinConfigManager.load_config(basin_file)
                if isinstance(basin_records, list) and len(basin_records) > 0:
                    config["nOutlet"] = len(basin_records)
                    if len(basin_records) > 1:
                        config["iFlag_multiple_outlet"] = 1
        except Exception:
            pass

    return config


# Alias for cross-project naming consistency with pyearthmesh
load_workflow_config = load_pyflowline_config


def dump_pyflowline_config(
    config: dict[str, Any],
    filename: str | Path,
    fmt: str | None = None,
) -> None:
    """Dump a configuration dictionary to a YAML or JSON file."""
    _dump_config_file(config, Path(filename), fmt=fmt)


# Alias for cross-project naming consistency with pyearthmesh
dump_workflow_config = dump_pyflowline_config


# ==============================================================================
# Model Loading & Workflow Factory (formerly read_configuration_file)
# ==============================================================================


def pyflowline_read_configuration_file(
    sFilename_configuration_in: str | Path,
    iFlag_standalone_in: int | None = 1,
    iFlag_create_directory_in: int | None = None,
    iFlag_use_mesh_dem_in: int | None = None,
    iCase_index_in: int | None = None,
    iResolution_index_in: int | None = None,
    dResolution_degree_in: float | None = None,
    dResolution_meter_in: float | None = None,
    sMesh_type_in: str | None = None,
    sModel_in: str | None = "pyflowline",
    sDate_in: str | None = None,
    sDggrid_type_in: str | None = None,
    sWorkspace_output_in: str | Path | None = None,
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
        flowlinecase | None: Configured flowlinecase object ready for simulation.
    """
    from pyflowline.classes.pycase import flowlinecase

    resolved_path = _resolve_config_path(sFilename_configuration_in)

    if not resolved_path.is_file():
        print(f"{sFilename_configuration_in} does not exist")
        return None

    # Load configuration (handles YAML comments, relative paths, platform paths)
    aConfig: dict[str, Any] = load_pyflowline_config(resolved_path)

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


# Alias for cross-project consistency with mpas_land_mesh / pyearthmesh
load_configuration = pyflowline_read_configuration_file


def create_workflow(config_or_path: str | Path | dict[str, Any], **kwargs):
    """Factory function: load configuration and return an instantiated flowlinecase."""
    from pyflowline.classes.pycase import flowlinecase

    if isinstance(config_or_path, (str, Path)):
        return load_configuration(config_or_path, **kwargs)
    elif isinstance(config_or_path, dict):
        if "workflow" in config_or_path or "resolutions" in config_or_path:
            config = compile_workflow(config_or_path)
        else:
            config = config_or_path.copy()
    elif isinstance(config_or_path, WorkflowSpec):
        config = compile_workflow(config_or_path)
    else:
        raise TypeError("config_or_path must be a path string, Path, or dictionary")

    # Apply any kwargs overrides
    for k, v in kwargs.items():
        config[k] = v

    iFlag_create_directory = kwargs.get("iFlag_create_directory_in", None)
    return flowlinecase(config, iFlag_create_directory_in=iFlag_create_directory)


class BasinConfigManager:
    """Configuration manager for basin objects supporting YAML and JSON."""

    @staticmethod
    def get_default_config():
        """Returns a dictionary with all default basin configuration values."""
        return DEFAULT_BASIN_CONFIG.copy()

    @staticmethod
    def create_template_config(basin_id=1, custom_values=None):
        """Create a configuration dictionary for a basin with default values."""
        config = BasinConfigManager.get_default_config()
        config["lBasinID"] = basin_id

        sWorkspace_input = Path.cwd()
        config["sFilename_dam"] = str(sWorkspace_input / config["sFilename_dam"])
        config["sFilename_flowline_filter"] = str(
            sWorkspace_input / config["sFilename_flowline_filter"]
        )
        config["sFilename_flowline_raw"] = str(sWorkspace_input / config["sFilename_flowline_raw"])
        config["sFilename_flowline_topo"] = str(
            sWorkspace_input / config["sFilename_flowline_topo"]
        )

        sWorkspace_output = Path.cwd()
        basin_folder = f"{basin_id:09d}"
        config["sWorkspace_output_basin"] = str(sWorkspace_output / basin_folder)

        if custom_values:
            for key, value in custom_values.items():
                if isinstance(value, Path):
                    value = str(value)
                config[key] = value

        return config

    @staticmethod
    def create_multiple_basins(sFilename_basins, nBasin, custom_values=None):
        """Generate configuration for multiple basins and save to YAML or JSON file."""
        from pyflowline.classes.basin import pybasin

        sFilename_basins_path = Path(sFilename_basins).expanduser().resolve()
        sFilename_basins_path.parent.mkdir(parents=True, exist_ok=True)

        aBasin_out = []
        basin_dicts = []
        for i in range(nBasin):
            aConfig_basin = BasinConfigManager.create_template_config(
                basin_id=i + 1, custom_values=custom_values
            )
            pBasin = pybasin(aConfig_basin, iFlag_create_directory_in=0)
            aBasin_out.append(pBasin)
            basin_dicts.append(json.loads(pBasin.tojson()))

        # Determine format from extension
        ext = sFilename_basins_path.suffix.lower()
        fmt = "json" if ext == ".json" else "yaml"
        _dump_config_file(basin_dicts, sFilename_basins_path, fmt=fmt)

        return aBasin_out

    @staticmethod
    def load_config(filename: str | Path) -> list[dict[str, Any]] | dict[str, Any]:
        """Load basin configuration from YAML or JSON file."""
        resolved = _resolve_config_path(
            filename,
            candidates=("basins.yaml", "basins.yml", "basins.json"),
        )
        data = _parse_config_file(resolved)
        if isinstance(data, list):
            return [_normalize_basin_dict(b) for b in data]
        if isinstance(data, dict):
            return _normalize_basin_dict(data)
        return data

    @staticmethod
    def dump_config(config: Any, filename: str | Path, fmt: str | None = None) -> None:
        """Dump basin configuration to YAML or JSON file."""
        _dump_config_file(config, Path(filename), fmt=fmt)


class FlowlineConfigManager:
    """Configuration manager for pyflowline supporting YAML and JSON."""

    @staticmethod
    def get_default_config():
        """Returns a dictionary with all default configuration values."""
        return DEFAULT_FLOWLINE_CONFIG.copy()

    @staticmethod
    def create_template_config(output_filename, custom_values=None):
        """Create a configuration file with default values, optionally customized."""
        config = FlowlineConfigManager.get_default_config()
        config["sWorkspace_input"] = str(Path.cwd())
        config["sWorkspace_output"] = str(Path.cwd())
        config["sFilename_model_configuration"] = str(output_filename)

        if custom_values:
            for key, value in custom_values.items():
                if isinstance(value, Path):
                    value = str(value)
                config[key] = value

        dest = Path(output_filename).expanduser().resolve()
        dest.parent.mkdir(parents=True, exist_ok=True)

        ext = dest.suffix.lower()
        fmt = "json" if ext == ".json" else "yaml"
        _dump_config_file(config, dest, fmt=fmt)

        return config

    @staticmethod
    def load_config(filename: str | Path) -> dict[str, Any]:
        """Load pyflowline configuration from YAML or JSON file."""
        return load_pyflowline_config(filename)

    @staticmethod
    def dump_config(config: dict[str, Any], filename: str | Path, fmt: str | None = None) -> None:
        """Dump pyflowline configuration to YAML or JSON file."""
        dump_pyflowline_config(config, filename, fmt=fmt)


class JigsawConfigManager:
    """Configuration manager for JIGSAW that handles YAML and JSON serialization."""

    @staticmethod
    def get_default_config():
        """Returns a dictionary with all default JIGSAW configuration values."""
        return {
            "ncolumn_space": 360,
            "nrow_space": 180,
            "dSpac_value": 100.0,
            "iFlag_geom": False,
            "iFlag_spac": False,
            "iFlag_init": False,
            "iFlag_opts": False,
            "iFlag_spac_ocean": False,
            "iFlag_spac_land": True,
            "iFlag_geom_dam": False,
            "iFlag_spac_dam": False,
            "iFlag_geom_city": False,
            "iFlag_spac_city": False,
            "iFlag_geom_river_network": False,
            "iFlag_spac_river_network": False,
            "iFlag_geom_coastline": False,
            "iFlag_spac_coastline": False,
            "iFlag_geom_watershed_boundary": False,
            "iFlag_spac_watershed_boundary": False,
            "iFlag_geom_lake_boundary": False,
            "iFlag_spac_lake_boundary": False,
            "dResolution_land": 45.0,
            "dResolution_dam": 4.0,
            "dResolution_city": 4.0,
            "dResolution_river_network": 4.0,
            "dResolution_coastline": 4.0,
            "dResolution_watershed_boundary": 4.0,
            "dResolution_lake_boundary": 4.0,
            "geom_mshID": "ellipsoid-mesh",
            "spac_mshID": "ellipsoid-grid",
            "FULL_SPHERE_RADIUS": 6371.0,
            "dhdx_lim": 0.25,
            "sFilename_dam_vector": None,
            "sFilename_dam_raster": None,
            "sFilename_city_vector": None,
            "sFilename_city_raster": None,
            "sFilename_river_network_vector": None,
            "sFilename_river_network_raster": None,
            "sFilename_coastline_vector": None,
            "sFilename_coastline_raster": None,
            "sFilename_watershed_boundary_vector": None,
            "sFilename_watershed_boundary_raster": None,
            "sFilename_lake_boundary_vector": None,
            "sFilename_lake_boundary_raster": None,
            "hfun_hmax": "inf",
            "hfun_hmin": 0.0,
            "hfun_scal": "absolute",
            "mesh_dims": 2,
            "bisection": -1,
            "optm_qlim": 0.95,
            "optm_iter": 32,
            "optm_qtol": 1.0e-05,
            "mesh_rad2": 1.5,
            "mesh_rad3": 2.0,
            "mesh_eps1": 0.333,
            "mesh_eps2": 0.333,
            "mesh_top": 1,
            "mesh_iter": 3,
            "verbosity": 0,
            "geom_file": None,
            "hfun_file": None,
            "mesh_file": None,
            "mesh_kern": "delfront",
            "optm_kern": "odt+dqdx",
            "geom_feat": True,
            "mesh_type": "euclidean-mesh",
            "output_formats": ["vtk", "gmsh"],
        }

    @staticmethod
    def create_template_config(output_filename, custom_values=None):
        """Create a JIGSAW configuration file in YAML or JSON."""
        config = JigsawConfigManager.get_default_config()
        if custom_values:
            for key, value in custom_values.items():
                if isinstance(value, Path):
                    value = str(value)
                config[key] = value

        dest = Path(output_filename).expanduser().resolve()
        dest.parent.mkdir(parents=True, exist_ok=True)

        ext = dest.suffix.lower()
        fmt = "json" if ext == ".json" else "yaml"
        _dump_config_file(config, dest, fmt=fmt)
        return config

    @staticmethod
    def load_config(filename: str | Path) -> dict[str, Any]:
        """Load a JIGSAW configuration from a YAML or JSON file."""
        resolved = _resolve_config_path(
            filename,
            candidates=("jigsaw.yaml", "jigsaw.yml", "jigsaw.json"),
        )
        return _parse_config_file(resolved)

    @staticmethod
    def dump_config(config: dict[str, Any], filename: str | Path, fmt: str | None = None) -> None:
        """Dump a JIGSAW configuration to a YAML or JSON file."""
        _dump_config_file(config, Path(filename), fmt=fmt)


# Functions for backwards compatibility
def create_basin_template_configuration_file(sFilename_basins, nBasin):
    """Generate basin configuration using the BasinConfigManager."""
    return BasinConfigManager.create_multiple_basins(sFilename_basins, nBasin)


def create_pyflowline_template_configuration_file(sFilename_config, **kwargs):
    """Generate pyflowline config template file using parameter keywords."""
    config = FlowlineConfigManager.create_template_config(sFilename_config, kwargs)

    from pyflowline.classes.pycase import flowlinecase

    oModel = flowlinecase(config, iFlag_create_directory_in=0)

    dest_path = Path(sFilename_config)
    ext = dest_path.suffix.lower()
    basin_ext = ext if ext in (".yaml", ".yml", ".json") else ".yaml"
    sFilename_basins = str(dest_path.parent / (dest_path.stem + "_basins" + basin_ext))

    aBasin = create_basin_template_configuration_file(sFilename_basins, config.get("nOutlet", 1))

    oModel.aBasin = aBasin
    oModel.sFilename_basins = sFilename_basins

    # Export configuration in the target format
    fmt = "json" if ext == ".json" else "yaml"
    dump_pyflowline_config(config, sFilename_config, fmt=fmt)

    return oModel


def create_jigsaw_template_configuration_file(sFilename_config, **kwargs):
    """Generate JIGSAW config template file using parameter keywords."""
    return JigsawConfigManager.create_template_config(sFilename_config, custom_values=kwargs)


# Backwards-compatibility stubs for legacy scripts
def change_config_key_value(*args: Any, **kwargs: Any) -> None:
    """Deprecated legacy stub. Modify configuration in memory instead."""
    pass


change_json_key_value = change_config_key_value
change_config_keys_values = change_config_key_value
change_json_keys_values = change_config_key_value


__all__ = [
    "DEFAULT_BASIN_CONFIG",
    "DEFAULT_FLOWLINE_CONFIG",
    "BasinConfigManager",
    "FlowlineConfigManager",
    "JigsawConfigManager",
    "WorkflowSpec",
    "compile_workflow",
    "create_basin_template_configuration_file",
    "create_jigsaw_template_configuration_file",
    "create_pyflowline_template_configuration_file",
    "create_workflow",
    "dump_pyflowline_config",
    "dump_workflow_config",
    "join_project_path",
    "load_configuration",
    "load_pyflowline_config",
    "load_workflow_config",
    "path_manager",
    "pyflowline_package_root",
    "pyflowline_project_root",
    "pyflowline_read_configuration_file",
    "root_path_from_path_manager_location",
    "root_path_from_pyflowline_package_root",
    "root_path_from_setup_file",
]
