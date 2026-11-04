#!/usr/bin/env python3
"""Yellow River (Huanghe) MPAS mesh PyFlowline simulation workflow.

This example demonstrates loading the YAML configuration using load_configuration
and orchestrating PyFlowline model execution (setup, simplification, mesh generation,
and topological reconstruction).
"""

import argparse
import logging
import sys
from pathlib import Path

from pyflowline.configuration import load_configuration, path_manager

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        "-c",
        default=None,
        help="Path to workflow YAML configuration file (default: pyflowline_huanghe_mpas.yaml)",
    )
    parser.add_argument(
        "--simulate",
        action="store_true",
        help="Run simulation pipeline (simplification, mesh generation, topological reconstruction)",
    )
    parser.add_argument(
        "--visualize",
        action="store_true",
        help="Generate flowline visualization map",
    )
    args = parser.parse_args()

    # ── Step 1: Resolve configuration path ────────────────────────────────────
    if args.config:
        sFilename_workflow = Path(args.config).expanduser().resolve()
    else:
        sFilename_workflow = Path(__file__).with_name("pyflowline_huanghe_mpas.yaml")

    logger.info("Step 1: Loading configuration from: %s", sFilename_workflow)
    case = load_configuration(sFilename_workflow)
    if case is None:
        logger.error("Failed to load configuration from: %s", sFilename_workflow)
        sys.exit(1)

    logger.info("   * Region: %s", case.sRegion)
    logger.info("   * Mesh type: %s", case.sMesh_type)
    logger.info("   * Case ID: %s", case.sCase)
    logger.info("   * Output workspace: %s", case.sWorkspace_output)
    if case.sFilename_mpas_mesh_netcdf:
        logger.info("   * MPAS mesh NetCDF: %s", case.sFilename_mpas_mesh_netcdf)
    if case.sFilename_mesh_boundary:
        logger.info("   * Mesh boundary: %s", case.sFilename_mesh_boundary)

    # ── Step 2: Model setup ───────────────────────────────────────────────────
    logger.info("Step 2: Setting up PyFlowline model...")
    case.setup()

    # ── Step 3: Visualization ─────────────────────────────────────────────────
    if args.visualize or getattr(case, "iFlag_visualization", 0) == 1:
        logger.info("Step 3: Generating initial filtered flowline visualization...")
        case.plot(sVariable_in="flowline_filter", sFilename_output_in="filter_flowline.png")

    # ── Step 4: Simulation Execution ──────────────────────────────────────────
    if args.simulate or getattr(case, "iFlag_simulation", 0) == 1:
        logger.info("Step 4: Running flowline simplification...")
        case.flowline_simplification()

        logger.info("Step 5: Generating mesh...")
        aCell = case.mesh_generation()

        logger.info("Step 6: Reconstructing topological relationships...")
        case.reconstruct_topological_relationship(aCell)

        logger.info("Step 7: Exporting simulation results...")
        case.export()

    logger.info("Done.")


if __name__ == "__main__":
    main()
