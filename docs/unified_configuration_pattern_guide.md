# Unified Configuration Architecture Pattern Guide

> **A reusable architectural blueprint for scientific modeling, geospatial, and computational Python packages.**
>
> *Standardized across `pyflowline`, `pyearthmesh`, and `mpas_land_mesh`.*

---

## 1. Overview & Core Philosophy

Scientific simulation frameworks often face two competing demands:
1. **Human Usability:** Domain scientists need readable, structured configuration files where assumptions, data sources, thresholds, and parameters can be documented directly with `# inline comments`.
2. **Backend Execution Rigor:** Computational engines require strict, low-level data structures, numeric flags, resolved absolute paths, and normalized types.

Historically, projects frequently suffer from **configuration module sprawl**:
- Separate `path_manager.py` for root detection.
- Separate `read_configuration_file.py` for loading.
- Separate `change_config_key_value.py` mutating JSON files on disk.
- Fragile coupling between parent domain cases and child features/basins.

The **Unified Configuration Manager Pattern** consolidates all path discovery, format parsing, cross-platform path resolution, compiler transformation, and model instantiation into **one coherent, authoritative module** (`config_manager.py`).

```text
┌────────────────────────────────────────────────────────────────────────┐
│                        User Workflow YAML File                         │
│  (Human-friendly, inline comments #, sectioned: workflow/paths/flags)  │
└───────────────────────────────────┬────────────────────────────────────┘
                                    │
                                    ▼
┌────────────────────────────────────────────────────────────────────────┐
│            Unified Config Manager (config/config_manager.py)           │
│  ├── 1. Package & Project Root Discovery (importlib.resources)         │
│  ├── 2. Auto-Discovery (.yaml -> .yml -> .json / directory defaults)   │
│  ├── 3. OS Path Resolution (Darwin / Linux / Windows / default)        │
│  ├── 4. Relative Path Resolution (relative to config file directory)   │
│  ├── 5. Two-Tier Compiler Boundary (compile_workflow)                  │
│  └── 6. Runtime Overrides (in-memory kwargs without disk mutation)     │
└───────────────────────────────────┬────────────────────────────────────┘
                                    │
                                    ▼
┌────────────────────────────────────────────────────────────────────────┐
│                        Engine Execution Object                         │
│               (flowlinecase / BaseWorkflow / CaseAdapter)              │
│               setup() -> run() -> export() -> plot()                   │
└────────────────────────────────────────────────────────────────────────┘
```

---

## 2. Key Distinction: `load_workflow_config` vs. `load_configuration`

A fundamental distinction in this architecture is the clean separation between **Data Compilation** (`dict`) and **Execution Instantiation** (`object`):

| Feature / Aspect | `load_workflow_config(path)` | `load_configuration(path, **kwargs)` |
| :--- | :--- | :--- |
| **Return Type** | **Compiled Dictionary** (`dict[str, Any]`) | **Instantiated Engine Object** (`flowlinecase`, `Workflow`) |
| **Pipeline Stage** | **Data Layer** (Parses, validates, resolves paths, compiles) | **Execution Layer** (Wraps compiled config into simulation object) |
| **Model State** | Passive data structure (no simulation state) | Active OOP object with simulation state, arrays, and caches |
| **Callable Methods** | None (read/write dictionary keys: `config["sCase"]`) | Executable methods: `.setup()`, `.run()`, `.export()`, `.plot()` |
| **Runtime Overrides** | Returns compiled dictionary from file | Accepts in-memory `**kwargs` overrides (`dResolution_meter_in=...`) |
| **Directory I/O** | Does not create workspace directories | Automatically prepares/creates workspace directories |
| **Primary Use Cases** | Test suites, HPC batch job generators, inspecting metadata | Workflow runner scripts (`run_simulation_*.py`), interactive notebooks |
| **Common Aliases** | `load_pyflowline_config` | `create_workflow` |

### Side-by-Side Usage Comparison

```python
# ==============================================================================
# 1. load_workflow_config: Pure Data Pipeline (Dictionary)
# ==============================================================================
from pyflowline.configuration import load_workflow_config

# Loads YAML/JSON, resolves OS/relative paths, compiles into backend dictionary
config = load_workflow_config("pyflowline_susquehanna_mpas.yaml")

print(type(config))             # <class 'dict'>
print(config["sMesh_type"])     # "mpas"
print(config["dResolution_meter"]) # 5000.0
# Ideal for: inspecting paths, generating SLURM job scripts, unit testing assertions


# ==============================================================================
# 2. load_configuration: Active Simulation Pipeline (Object)
# ==============================================================================
from pyflowline.configuration import load_configuration

# Loads config AND instantiates the flowlinecase engine with optional overrides
case = load_configuration(
    "pyflowline_susquehanna_mpas.yaml",
    dResolution_meter_in=2500.0, # Runtime in-memory override
)

print(type(case))               # <class 'pyflowline.classes.pycase.flowlinecase'>
case.setup()                    # Sets up domain, flowline vectors, coordinates
case.flowline_simplification()  # Simplifies river reaches
aCell = case.mesh_generation()  # Generates MPAS Voronoi mesh
case.reconstruct_topological_relationship(aCell)
case.export()                   # Writes final GeoJSON, CSV, topology files
```

---

## 3. Standard Directory Layout

Instead of scattering multiple helper modules, each project adopts a minimal, identical layout:

```text
my_package/
├── __init__.py                 # Exposes: import my_package.config as config
└── config/                     # (or configuration/)
    ├── __init__.py             # Public API re-exports + legacy shims
    └── config_manager.py       # The SINGLE source of truth
```

### Module Alias in `my_package/__init__.py`
Enable cross-project import parity (`import my_package.config` or `from my_package.config import ...`):

```python
from __future__ import annotations
import sys
from my_package import config
from my_package.config import load_configuration, create_workflow

# Cross-project alias
sys.modules.setdefault("my_package.config", config)
sys.modules.setdefault("my_package.config.config_manager", config.config_manager)
```

---

## 4. Key Architectural Pillars

### 4.1 Pillar 1: Dual-Tier Configuration Specification
Configurations are organized into high-level logical sections rather than a flat dictionary of 50+ obscure keys:

```yaml
# ==============================================================================
# Simulation Case Specification (YAML)
# ==============================================================================
workflow:
  date: "20261101"                              # Run date (YYYYMMDD)
  mesh_type: mpas                               # Target mesh or engine family
  case_index: 1                                 # Run index
  model: pyflowline                             # Simulation engine
  standalone: true                              # Standalone execution
  debug: false                                  # Diagnostic logging

resolutions:
  degree: 0.05                                  # Resolution in degrees
  meter: 3000.0                                 # High-resolution corridor spacing (m)
  accumulation_threshold: 100000.0              # Upstream area threshold

region:
  bbox:
    left: 95.0
    right: 120.0
    bottom: 32.0
    top: 42.0

flags:
  simulation: true                              # Enable simulation execution pipeline
  visualization: false                          # Enable map visualization output
  save_mesh: true                               # Export results to disk

paths:
  workspace_input:
    Darwin: /Users/user/data/project
    Linux: /public/home/user/data/project
    default: ./input
  workspace_output:
    Darwin: /Users/user/scratch/project
    Linux: /data/share/user/project
    default: ./output
  mesh_netcdf: ./culled_mesh.nc
  basins: ./project_basins.yaml                 # Child specification file
```

---

### 4.2 Pillar 2: Cross-Platform Path Mapping
Modelers execute runs locally on macOS laptops and deploy batch simulations on Linux HPC clusters. The loader evaluates `platform.system()` dynamically:

```python
def _select_platform_path(value: Any) -> str:
    """Resolve platform-specific dictionary or string path."""
    if isinstance(value, str):
        return value
    if isinstance(value, Path):
        return str(value)
    if not isinstance(value, dict):
        raise ValueError("Path values must be strings or platform-specific mappings")

    current_platform = platform.system()  # 'Darwin', 'Linux', 'Windows'
    selected = value.get(current_platform, value.get("default", ""))
    return str(selected)
```

---

### 4.3 Pillar 3: Flexible 1-to-N Relationships (Split vs. Inlined)
When a domain simulation includes multiple sub-domains, watersheds, basins, or boundary features:
- **Modular / Split Model (Recommended for Multi-Entity):**
  Parent case points to child specification via `paths.basins: ./basins.yaml`. The child file can be reused across multiple mesh setups without duplicate definitions.
- **Unified / Inlined Model:**
  Child entities are listed directly inside the parent YAML under `basins: [...]`.

The loader transparently supports both, auto-detecting entity counts:

```python
# In compile_workflow:
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
```

---

### 4.4 Pillar 4: In-Memory Overrides (No Disk Mutation)
Never mutate YAML or JSON files on disk before running a case (which strips user comments and creates file race conditions). Instead, pass runtime parameter overrides:

```python
# Pass runtime overrides directly to load_configuration
case = load_configuration(
    "workflow_config.yaml",
    dResolution_meter_in=5000.0,
    sWorkspace_output_in="./custom_output",
)
case.setup()
case.run()
```

---

### 4.5 Pillar 5: First-Class Template Generation (`create_model_configuration`)
Rather than maintaining loose, ad-hoc template generator scripts in `examples/` that drift, suffer from typos, or get left unmaintained, template generation is an intrinsic capability of `config_manager.py`. It produces either human-readable, sectioned YAML with paired basins or flat/sectioned JSON based on the destination file extension or `fmt` argument:

```python
from my_package import create_model_configuration

# 1. Generate YAML template (primary)
create_model_configuration("workflow_config.yaml", mesh_type="hexagon")

# 2. Generate JSON template (programmatic)
create_model_configuration("workflow_config.json", mesh_type="mpas")
```

---

## 5. Reusable Implementation Blueprint

Below is the template for [`config_manager.py`](file:///Users/changliao/workspace/python/pyflowline/pyflowline/configuration/config_manager.py) ready to be dropped into any project:

```python
"""Unified Configuration Manager Blueprint.

Single-source configuration, path resolution, two-tier compiler,
and workflow factory.
"""

from __future__ import annotations

import datetime
import importlib.resources as importlib_resources
import json
import os
import platform
import sys
from pathlib import Path
from typing import Any

try:
    import yaml
    HAVE_YAML = True
except ImportError:
    HAVE_YAML = False


# ==============================================================================
# 1. Package & Project Root Path Discovery
# ==============================================================================

def package_root(package_name: str = "my_package") -> Path:
    """Locate package installation root directory."""
    try:
        return Path(str(importlib_resources.files(package_name))).resolve()
    except Exception:
        return Path(__file__).resolve().parents[1]


def project_root() -> Path:
    """Find repository root by searching upward for project marker files."""
    this_path = Path(__file__).resolve()
    for parent in this_path.parents:
        if (parent / "setup.py").exists() or (parent / "pyproject.toml").exists() or (parent / ".git").exists():
            return parent
    return package_root().parent


def join_project_path(*args: str) -> Path:
    """Join path components relative to the project root."""
    return project_root().joinpath(*args)


# ==============================================================================
# 2. Path Resolution & Discovery Helpers
# ==============================================================================

def _select_platform_path(value: Any) -> str:
    """Resolve platform-specific dictionary or plain string path."""
    if isinstance(value, str):
        return value
    if isinstance(value, Path):
        return str(value)
    if not isinstance(value, dict):
        raise ValueError("Path values must be strings or platform-specific mappings")

    current_platform = platform.system()
    selected = value.get(current_platform, value.get("default", ""))
    return str(selected)


def _resolve_config_path(
    path_like: str | Path,
    candidates: tuple[str, ...] = ("workflow_config.yaml", "workflow_config.yml", "workflow_config.json"),
) -> Path:
    """Resolve file path supporting directory lookup and extensionless stems."""
    path = Path(path_like).expanduser().resolve()

    if path.is_file():
        return path

    if path.is_dir():
        for candidate in candidates:
            cand_path = path / candidate
            if cand_path.is_file():
                return cand_path
        raise FileNotFoundError(f"Directory {path} contains none of: {candidates}")

    # Check extension stems
    for ext in (".yaml", ".yml", ".json"):
        with_ext = path.with_suffix(ext)
        if with_ext.is_file():
            return with_ext

    return path


def _parse_config_file(path: Path) -> dict[str, Any] | list[Any]:
    """Parse YAML or JSON configuration file."""
    if not path.is_file():
        raise FileNotFoundError(f"Configuration file not found: {path}")

    text = path.read_text(encoding="utf-8")
    ext = path.suffix.lower()

    if ext in (".yaml", ".yml"):
        if not HAVE_YAML:
            raise ImportError("PyYAML is required to parse .yaml/.yml files: pip install pyyaml")
        return yaml.safe_load(text) or {}

    return json.loads(text)


def _dump_config_file(data: Any, path: Path, fmt: str | None = None) -> None:
    """Dump configuration data to YAML or JSON."""
    ext = path.suffix.lower()
    use_yaml = (fmt == "yaml") or (fmt is None and ext in (".yaml", ".yml"))

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        if use_yaml:
            if not HAVE_YAML:
                raise ImportError("PyYAML is required to dump .yaml/.yml files: pip install pyyaml")
            yaml.dump(data, f, sort_keys=False, default_flow_style=False)
        else:
            json.dump(data, f, indent=4)


# ==============================================================================
# 3. Two-Tier Compiler Boundary
# ==============================================================================

def compile_workflow(raw_dict: dict[str, Any]) -> dict[str, Any]:
    """Compile human-friendly sectioned dictionary into backend-native dictionary."""
    native = {}

    workflow = raw_dict.get("workflow", {})
    if workflow:
        native["mesh_type"] = str(workflow.get("mesh_type", "default")).lower()
        native["case_index"] = int(workflow.get("case_index", 1))
        native["date"] = str(workflow.get("date", datetime.date.today().strftime("%Y%m%d")))
        native["debug"] = bool(workflow.get("debug", False))

    resolutions = raw_dict.get("resolutions", {})
    if resolutions:
        for k, v in resolutions.items():
            if v is not None:
                native[f"resolution_{k}"] = float(v)

    flags = raw_dict.get("flags", {})
    if flags:
        for k, v in flags.items():
            native[f"flag_{k}"] = int(bool(v)) if isinstance(v, bool) else v

    paths = raw_dict.get("paths", {})
    if paths:
        for k, v in paths.items():
            native[f"path_{k}"] = _select_platform_path(v)

    return native


# ==============================================================================
# 4. Loader & Model Factory
# ==============================================================================

def load_workflow_config(filename: str | Path) -> dict[str, Any]:
    """Load configuration file, resolve paths, and compile into backend representation."""
    resolved_path = _resolve_config_path(filename)
    raw = _parse_config_file(resolved_path)
    if not isinstance(raw, dict):
        raise ValueError(f"Configuration must be a dictionary, got: {type(raw)}")

    compiled = compile_workflow(raw)

    # Resolve relative paths against configuration file directory
    base_dir = resolved_path.parent
    for k, v in compiled.items():
        if k.startswith("path_") and isinstance(v, str) and v and not os.path.isabs(v):
            compiled[k] = str((base_dir / v).resolve())

    return compiled


def load_configuration(filename: str | Path, **kwargs) -> Any:
    """Load workflow configuration, apply overrides, and instantiate execution object."""
    config = load_workflow_config(filename)

    # Apply runtime kwargs overrides
    for k, v in kwargs.items():
        if v is not None:
            config[k] = v

    # Instantiate case engine (import locally to avoid circular dependencies)
    from my_package.classes.case import CaseEngine
    return CaseEngine(config)


def create_model_configuration(
    filename: str | Path,
    mesh_type: str = "default",
    case_index: int = 1,
    date: str | None = None,
    resolution_meter: float = 5000.0,
    resolution_degree: float = 0.5,
    region_name: str = "default",
    workspace_input: str | Path | None = None,
    workspace_output: str | Path | None = None,
    n_outlets: int = 1,
    fmt: str | None = None,
    sectioned: bool = True,
    **kwargs: Any,
) -> dict[str, Any]:
    """Generate a template workflow configuration file in YAML or JSON."""
    dest_path = Path(filename).expanduser().resolve()
    dest_path.parent.mkdir(parents=True, exist_ok=True)

    if fmt is None:
        ext = dest_path.suffix.lower()
        fmt = "json" if ext == ".json" else "yaml"
    else:
        fmt = fmt.lower()

    date_str = date or datetime.date.today().strftime("%Y%m%d")
    ws_in = str(workspace_input) if workspace_input else str(dest_path.parent / "input")
    ws_out = str(workspace_output) if workspace_output else str(dest_path.parent / "output")

    if fmt in ("yaml", "yml") and sectioned:
        config_data = {
            "workflow": {
                "mesh_type": mesh_type.lower(),
                "case_index": int(case_index),
                "date": date_str,
            },
            "resolutions": {
                "degree": float(resolution_degree),
                "meter": float(resolution_meter),
            },
            "region": {"name": region_name},
            "flags": {"simulation": True, "visualization": False},
            "paths": {
                "workspace_input": ws_in,
                "workspace_output": ws_out,
            },
        }
    else:
        config_data = {
            "sMesh_type": mesh_type,
            "iCase_index": case_index,
            "sDate": date_str,
            "dResolution_meter": resolution_meter,
            "dResolution_degree": resolution_degree,
            "sWorkspace_input": ws_in,
            "sWorkspace_output": ws_out,
        }
        config_data.update(kwargs)

    _dump_config_file(config_data, dest_path, fmt=fmt)
    return config_data


# Unified cross-project aliases
create_workflow = load_configuration
load_pyflowline_config = load_workflow_config
dump_workflow_config = _dump_config_file
```

---

## 6. Companion `__init__.py` Blueprint

```python
"""Configuration package interface."""
from __future__ import annotations

import sys
from my_package.config import config_manager
from my_package.config.config_manager import (
    compile_workflow,
    create_model_configuration,
    create_workflow,
    dump_workflow_config,
    join_project_path,
    load_configuration,
    load_workflow_config,
    package_root,
    project_root,
)

# Backwards-compatibility module shims
sys.modules.setdefault("my_package.config.config_manager", config_manager)
sys.modules.setdefault("my_package.config.path_manager", config_manager)
sys.modules.setdefault("my_package.config.read_configuration_file", config_manager)

__all__ = [
    "compile_workflow",
    "create_model_configuration",
    "create_workflow",
    "dump_workflow_config",
    "load_configuration",
    "load_workflow_config",
    "package_root",
    "project_root",
    "join_project_path",
    "config_manager",
]
```

---

## 7. Adoption Checklist for New / Existing Projects

When migrating a new project (e.g. `pyearthbuffer`, `pyearthriver`, `hexwatershed_utility`):

1. [ ] Create `my_package/config/` with `config_manager.py` and `__init__.py`.
2. [ ] Define the project-specific sections in `DEFAULT_CONFIG` (`workflow`, `resolutions`, `region`, `flags`, `paths`).
3. [ ] Define `_select_platform_path` for macOS / Linux / Windows compatibility.
4. [ ] Define `compile_workflow` mapping user-friendly names to backend-native keys.
5. [ ] Wire `load_configuration` and `create_workflow` to instantiate the project's primary case or pipeline object.
6. [ ] Delete legacy `path_manager.py`, `read_*.py`, and `change_*_key_value.py` files.
7. [ ] Add `my_package.config` alias in `my_package/__init__.py`.
8. [ ] Replace on-disk configuration mutations with clean `load_configuration(path, override_arg=val)`.
9. [ ] Provide `create_model_configuration` supporting dual YAML/JSON template generation.
