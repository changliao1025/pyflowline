# Yellow River (Huanghe) Basin Example

This directory contains PyFlowline configuration and simulation scripts for the Yellow River (Huanghe) basin using a variable-resolution MPAS Voronoi mesh.

## Configuration Files

- **`pyflowline_huanghe_mpas.yaml`**: Domain specification defining workflow orchestration, target multi-scale resolutions, bounding box, execution flags, and paths across platforms (Darwin, Linux, Windows).
- **`pyflowline_huanghe_basins.yaml`**: Basin specification with the main Yellow River outlet coordinate (at Dongying, Bohai Sea mouth: `119.1729°E, 37.7687°N`), upstream area thresholds, and cropped flowline shapefile paths.

## Running the Simulation

Execute the MPAS simulation script:
```bash
python run_simulation_mpas.py
```
