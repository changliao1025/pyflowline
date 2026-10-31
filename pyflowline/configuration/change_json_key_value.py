import shutil
import tempfile
from pathlib import Path
from typing import Any

from pyflowline.configuration.config_manager import (
    _dump_config_file,
    _parse_config_file,
    _resolve_config_path,
)


def change_config_key_value(
    sFilename_config_in: str | Path,
    sKey: str,
    new_value: Any,
    iFlag_basin_in: bool | int | None = None,
    iBasin_index_in: int | None = None,
) -> None:
    """Change the value associated with the specified key in a YAML or JSON file safely.

    Args:
        sFilename_config_in (str or Path): Path to YAML or JSON configuration file.
        sKey (str): Key whose value needs to be changed.
        new_value: New value to set for the specified key.
        iFlag_basin_in (bool or int, optional): Flag indicating if the key is for a basin.
        iBasin_index_in (int, optional): Index of the basin in a list of basins.
    """
    path = _resolve_config_path(sFilename_config_in)

    if isinstance(new_value, Path):
        new_value = str(new_value)

    data = _parse_config_file(path)

    # Determine whether it's a sectioned config or flat
    if iFlag_basin_in is None:
        if isinstance(data, dict):
            # Check if key belongs to a section in a sectioned config
            found = False
            for sec in ("workflow", "resolutions", "region", "paths", "flags", "jigsaw"):
                if sec in data and isinstance(data[sec], dict) and sKey in data[sec]:
                    data[sec][sKey] = new_value
                    found = True
                    break
            if not found:
                data[sKey] = new_value
        else:
            data[sKey] = new_value
    else:
        # Basin configuration
        idx = int(iBasin_index_in) if iBasin_index_in is not None else 0
        if isinstance(data, list):
            if idx < len(data) and isinstance(data[idx], dict):
                data[idx][sKey] = new_value
        elif isinstance(data, dict):
            if "basins" in data and isinstance(data["basins"], list):
                if idx < len(data["basins"]) and isinstance(data["basins"][idx], dict):
                    data["basins"][idx][sKey] = new_value
            else:
                data[sKey] = new_value

    ext = path.suffix.lower()
    suffix = ext if ext in (".yaml", ".yml", ".json") else ".yaml"
    fmt = "json" if suffix == ".json" else "yaml"

    with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=suffix) as temp_file:
        temp_path = Path(temp_file.name)

    _dump_config_file(data, temp_path, fmt=fmt)
    shutil.move(str(temp_path), str(path))


def change_config_keys_values(
    sFilename_config_in: str | Path,
    key_value_pairs: dict[str, Any],
    iFlag_basin_in: bool | int | None = None,
) -> None:
    """Change multiple key-value pairs in a YAML or JSON configuration file safely."""
    for key, value in key_value_pairs.items():
        change_config_key_value(sFilename_config_in, key, value, iFlag_basin_in=iFlag_basin_in)


# Backward-compatibility aliases
change_json_key_value = change_config_key_value
change_json_keys_values = change_config_keys_values
