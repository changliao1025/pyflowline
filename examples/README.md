# Examples

This directory contains examples demonstrating how to configure and execute PyFlowline simulations using modern YAML configurations.

## Overview

YAML is the primary configuration format for PyFlowline, offering human-friendly syntax, comment support, and structured sections (`workflow`, `resolutions`, `region`, `flags`, `paths`, and `basins`).

- **`create_model_configuration.py`**: Generates a template YAML configuration file for a given mesh type and resolution. You can edit this file and then use it for actual simulations.
  ```bash
  python create_model_configuration.py --sMesh_type hexagon --dResolution_meter 50000 --sDate 20230101
  ```

- **`susquehanna/run_simulation_<meshtype>.py`**: Demonstrates how to load a YAML domain configuration, configure basin outlets, and run simulations for different mesh families:
  - `run_simulation_hexagon.py`: Hexagonal mesh simulation with `pyflowline_susquehanna_hexagon.yaml`
  - `run_simulation_square.py`: Square mesh simulation with `pyflowline_susquehanna_square.yaml`
  - `run_simulation_latlon.py`: Regular latitude-longitude mesh with `pyflowline_susquehanna_latlon.yaml`
  - `run_simulation_mpas.py`: Variable-resolution MPAS Voronoi mesh with `pyflowline_susquehanna_mpas.yaml` and programmatic parameter updates.

- **`huanghe/run_simulation_mpas.py`**: Clean, parameterless simulation orchestration for the Yellow River (Huanghe) basin using `pyflowline_huanghe_mpas.yaml` and `pyflowline_huanghe_basins.yaml`, driving an MPAS Voronoi mesh workflow via `load_configuration`.
