# Two-Tier Compiler Pattern for PyFlowline Workflow Architecture

> **Status: ACTIVE SPECIFICATION (CORE ARCHITECTURE IMPLEMENTED & YAML-EXTENDED)**
> - **Core Compiler Boundary:** Implemented via [`compile_workflow()`](../pyflowline/configuration/config_manager.py) and [`WorkflowSpec`](../pyflowline/configuration/config_manager.py).
> - **Specification Formats:** **YAML (primary, comment-enabled, user-friendly)** and **JSON (backward-compatible, programmatic)** via [`load_pyflowline_config()`](../pyflowline/configuration/config_manager.py) and [`dump_pyflowline_config()`](../pyflowline/configuration/config_manager.py).
> - **Case & Basin Dual Support:** Both parent domain cases and child basin specifications natively support YAML with `# comments`.
> - **Supported Mesh Families:** Hexagon, Square, LatLon, MPAS, DGGRID (ISEA3H), and TIN.
> - **Test Suite:** Passing in [`tests/test_yaml_configuration.py`](../tests/test_yaml_configuration.py).

---

## 1. Purpose & Philosophy

Hydrologic river network mesh generation requires reconciling two competing demands:
1. **User Friendliness & Maintainability:** Domain scientists and hydrologists need readable configuration files where choices (e.g. *why* drainage accumulation threshold is 100,000 $m^2$, *why* stream burning topology is enabled, which SLURM queue is targeted, or how river reaches are pruned) can be documented directly inline via `# comments`.
2. **Execution Rigor:** PyFlowline's internal case execution engine (`flowlinecase` and `pybasin`) requires dozens of strict, low-level parameters, matrix flags, and file paths.

The Two-Tier Compiler Pattern bridges this divide cleanly:

```text
User Intent (Human-Friendly YAML with Comments)
        │
        ├── Parent Case: workflow, resolutions, region, flags, paths
        └── Child Basins: outlet coordinates, dams, flowline vectors (split or inline)
        │
        ▼  [Tier 1: load_pyflowline_config / WorkflowSpec]
Domain Specification
        │
        ▼  [Compiler: compile_workflow()]
Backend-Native Configuration
        │
        ▼  [Tier 2: flowlinecase & pybasin]
Execution & Simulation (setup -> simplification -> mesh_generation -> export)
```

---

## 2. Why YAML is the Primary Tier 1 Specification Format

While standard JSON is effective for machine-to-machine exchange, it poses severe usability barriers for human modelers:
- **No Comments Allowed:** Standard JSON strictly prohibits comments (`#` or `//`). Modelers cannot explain parameter justifications, cite data sources, or temporarily annotate settings.
- **Syntax Clutter:** Quoted strings on every key, bracket balancing, and syntax errors caused by trailing commas make managing configuration files frustrating.
- **Separation of Files:** In JSON, parent case files and child basin files were forced to be separate files connected via brittle absolute paths.

By adopting **YAML** as the primary human-facing format:
- **Inline and Block Comments (`# ...`):** Users document domain boundaries, outlet locations, and thresholds directly where they are defined.
- **Clean Indentation-Based Syntax:** Eliminates bracket noise and trailing comma syntax errors.
- **Inlined or Modular Basin Definitions:** Basins can either be referenced in external files (`paths.basins: ./basins.yaml`) or embedded directly within the parent case YAML under `basins:`.
- **Full Backward Compatibility on Load:** JSON is a formal syntactic subset of YAML; PyFlowline accepts `.yaml`, `.yml`, and legacy `.json` files interchangeably, while the export API provides a clean, unified interface via `pyflowline_export_config()` and `pyflowline_export_basin_config()` (defaulting to YAML).

---

## 3. Key Challenges and Architectural Solutions

### 3.1 Challenge 1: Hierarchical Relational Structure (Case vs. Basin Coupling)
- **Challenge:** PyFlowline models have a 1-to-N parent-child relationship between the case (domain bounding box, mesh generator, output directory) and basins (individual outlet coordinates, stream order filters, dam locations). Managing multiple decoupled JSON files with hardcoded file paths caused frequent broken link errors.
- **Solution:**
  - The YAML compiler supports **both split and unified models**:
    - *Unified / Inlined:* Users can define `basins:` directly as a list inside the case YAML file.
    - *Modular / Split:* Users can reference external YAML or JSON basin files via `paths.basins`.
  - Relative paths inside both case and basin files are automatically resolved relative to the configuration file's directory.

### 3.2 Challenge 2: Portable Paths Across Different Computing Environments
- **Challenge:** Modelers run simulations locally on macOS/Windows laptops and deploy production runs on Linux HPC clusters (e.g. SLURM/Compy/Constance). Hardcoded file paths in legacy JSON configs caused cross-platform failures.
- **Solution:**
  - Support platform-specific dictionaries:
    ```yaml
    paths:
      mesh_netcdf:
        Linux: /compyfs/user/data/lnd_cull_mesh.nc
        Darwin: /Users/user/data/lnd_cull_mesh.nc
        default: ./lnd_cull_mesh.nc
    ```
  - Automatically evaluate the host operating system (`platform.system()`) at runtime.

### 3.3 Challenge 3: YAML Type Coercion Nuances
- **Challenge:** YAML parsers may implicitly convert numeric dates (e.g. `20220202`) into integers, case identifiers (`00000001`) into octal/integers, or boolean flags into unexpected types.
- **Solution:**
  - `compile_workflow()` enforces strict type normalization: date values are explicitly string-formatted (`sDate`), case indices are converted to integers, coordinates and thresholds are validated as floats.

### 3.4 Challenge 4: Preserving Comments During Programmatic Modification
- **Challenge:** Scripts frequently update workspace paths or parameters before execution using helper functions like `change_json_key_value`. Standard dumping re-writes files without original comments.
- **Solution:**
  - Recommend modifying parameters in-memory via `oPyflowline.pyflowline_change_model_parameter(...)` or runtime overrides in `pyflowline_read_configuration_file(..., sWorkspace_output_in=...)`.
  - Upgraded `change_config_key_value` to support both YAML and JSON files safely.

---

## 4. Layered Architecture

```text
User File (workflow_config.yaml / case.yaml / case.json)
        │
        ▼
Config Loader: load_pyflowline_config()
  ├── Resolves relative paths against file directory
  ├── Resolves OS-specific paths (Linux, Darwin, Windows)
  ├── Strips / parses comments seamlessly via yaml.safe_load()
  ├── Auto-discovers stem / directory files
  └── Normalizes run date (YYYYMMDD)
        │
        ▼
Tier 1 Spec: WorkflowSpec
  ├── workflow: mesh_type, case_index, model, job, date, standalone
  ├── resolutions: degree, meter, resolution_index, accumulation_threshold
  ├── region: bbox, outlet (longitude, latitude)
  ├── flags: create_mesh, save_mesh, simplification, intersect, use_mesh_dem
  ├── paths: workspace_input, workspace_output, mesh_boundary, dem, basins
  └── basins: inlined or resolved list of basin dicts
        │
        ▼
Compiler Boundary: compile_workflow(spec_or_config)
  ├── Maps sectioned keys to backend-native keys
  └── Normalizes basin dictionaries
        │
        ▼
Tier 2 Compiled Representation (Backend-Native Config Dict)
  ├── sMesh_type, iCase_index, dResolution_meter, iFlag_create_mesh, ...
  └── aBasin (list of normalized basin dictionaries)
        │
        ▼
Runtime Case Adapter: flowlinecase & pybasin
  ├── case.setup()
  ├── case.flowline_simplification()
  ├── case.mesh_generation()
  ├── case.reconstruct_topological_relationship()
  └── case.export()
```

---

## 5. Commented Workflow YAML Examples

### 5.1 Case Specification (`pyflowline_susquehanna_hexagon.yaml`)

```yaml
# ==============================================================================
# PyFlowline Case Specification: Hexagonal Mesh (YAML)
# ==============================================================================
workflow:
  date: "20220113"                              # Simulation date (YYYYMMDD)
  mesh_type: hexagon                            # Target mesh: hexagon
  case_index: 1                                 # Simulation run index
  model: pyflowline                             # Model engine
  job: hex                                      # Job name
  region: susquehanna                           # Study region
  standalone: true                              # Standalone case execution
  debug: false                                  # Diagnostic logging

resolutions:
  degree: 0.5                                   # Resolution in degrees
  meter: 50000.0                                # Hexagon cell diameter (50 km)
  accumulation_threshold: 100000.0              # Upstream drainage accumulation threshold
  small_river_threshold: 10000.0                # Small river pruning threshold

region:
  bbox:
    left: -79.44374                             # Western bounding longitude
    right: -74.24774                            # Eastern bounding longitude
    bottom: 39.00                               # Southern bounding latitude
    top: 43.00334                               # Northern bounding latitude
  outlet:
    longitude: -76.0093                         # Main basin outlet longitude
    latitude: 39.4620                           # Main basin outlet latitude

flags:
  create_mesh: true                             # Generate mesh geometry
  save_mesh: true                               # Export mesh to disk
  simplification: true                          # Simplify river geometry
  intersect: true                               # Intersect river reaches with mesh cells
  flowline: true                                # Enable river flowline processing
  use_mesh_dem: false                           # Sample elevation directly from mesh DEM
  rotation: false                               # Grid orientation rotation
  stream_burning_topology: true                 # Stream burning topological constraint
  resample_method: 2                            # Flowline resampling method

paths:
  workspace_input: .                            # Input directory (relative to this file)
  workspace_output: ../output                   # Output directory
  mesh_boundary: ./boundary_wgs.geojson         # Study area boundary mask
  spatial_reference: ./boundary_proj_buff.shp   # Spatial coordinate reference geometry
  dem: ./dem_buff_ext.tif                       # Digital Elevation Model raster
  basins: ./pyflowline_susquehanna_basins.yaml  # Child basin specification file
```

### 5.2 Basin Specification (`pyflowline_susquehanna_basins.yaml`)

```yaml
# ==============================================================================
# PyFlowline Basin Specification (YAML)
# ==============================================================================
- id: 1                                         # Basin unique identifier
  outlet:
    longitude: -76.0093                         # Outlet longitude in decimal degrees
    latitude: 39.4620                           # Outlet latitude in decimal degrees
  accumulation_threshold: 100000.0              # Upstream drainage area threshold
  small_river_threshold: 10000.0                # River reach length threshold for pruning
  flags:
    dam: 0                                      # Include dam locations (0: disabled, 1: enabled)
    debug: 1                                    # Diagnostic output flag
    disconnected: 0                             # Preserve disconnected reaches
  paths:
    dam: ./ICoM_dams.csv                        # Dam points dataset
    flowline_filter: ./flowline.geojson         # Filtered stream network vector
    flowline_raw: ./allflowline.shp             # Raw stream network dataset
    flowline_topo: ./flowline.csv               # River topology table
```

---

## 6. Verification & Acceptance Criteria

1. **Comment Support:** Both case and basin YAML files load cleanly when annotated with `# comments`.
2. **Dual-Format Parity:** Equivalent `.yaml` and `.json` files compile to identical backend configurations.
3. **Directory and Stem Auto-Discovery:** Directory paths and extensionless filenames resolve transparently.
4. **Clean Serialization:** `dump_pyflowline_config()` exports valid YAML or JSON files that reload identically.
5. **No Regressions:** All existing execution scripts, methods, and tests maintain 100% backward compatibility.
