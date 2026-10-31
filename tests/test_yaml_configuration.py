"""Unit tests for YAML configuration support in pyflowline.

Tests verify:
- Parsing YAML workflow and basin configs with comments
- Auto-discovery of extensions (.yaml, .yml, .json) and directories
- Dual-format parity between YAML and JSON
- Inline basin specifications inside YAML case files
- Safe key-value modification in YAML
- Exporting configurations to YAML and JSON
"""

import tempfile
import unittest
from pathlib import Path

from pyflowline.configuration.change_json_key_value import change_config_key_value
from pyflowline.configuration.config_manager import (
    BasinConfigManager,
    FlowlineConfigManager,
    load_pyflowline_config,
)
from pyflowline.configuration.read_configuration_file import pyflowline_read_configuration_file


class TestYAMLConfiguration(unittest.TestCase):
    """Test suite for YAML and JSON configuration loading and compiling."""

    def test_yaml_case_loading_with_comments(self):
        """Test loading a case YAML file containing inline and block comments."""
        yaml_content = """# Top-level case comment
workflow:
  date: "20230101" # Simulation date
  mesh_type: hexagon # Target mesh family
  case_index: 3
  model: pyflowline
  standalone: true
  debug: false

# Grid resolutions
resolutions:
  degree: 0.25 # Resolution in degrees
  meter: 25000.0 # Resolution in meters
  accumulation_threshold: 50000.0
  small_river_threshold: 5000.0

# Geographic bounds
region:
  bbox:
    left: -78.0 # West
    right: -75.0 # East
    bottom: 40.0 # South
    top: 42.0 # North
  outlet:
    longitude: -76.5
    latitude: 40.5

# Execution flags
flags:
  create_mesh: true
  save_mesh: true
  simplification: true
  intersect: true
  use_mesh_dem: false

paths:
  workspace_output: ./test_output
"""
        with tempfile.TemporaryDirectory() as tmpdir:
            yaml_path = Path(tmpdir) / "case_config.yaml"
            yaml_path.write_text(yaml_content, encoding="utf-8")

            config = load_pyflowline_config(yaml_path)
            self.assertEqual(config["sMesh_type"], "hexagon")
            self.assertEqual(config["iCase_index"], 3)
            self.assertEqual(config["sDate"], "20230101")
            self.assertEqual(config["dResolution_degree"], 0.25)
            self.assertEqual(config["dResolution_meter"], 25000.0)
            self.assertEqual(config["dLongitude_left"], -78.0)
            self.assertEqual(config["dLongitude_right"], -75.0)
            self.assertEqual(config["dLatitude_bot"], 40.0)
            self.assertEqual(config["dLatitude_top"], 42.0)
            self.assertEqual(config["dLongitude_outlet_degree"], -76.5)
            self.assertEqual(config["dLatitude_outlet_degree"], 40.5)
            self.assertEqual(config["iFlag_create_mesh"], 1)
            self.assertEqual(config["iFlag_standalone"], 1)

    def test_yaml_basin_loading_with_comments(self):
        """Test loading a basin YAML file containing comments and friendly structure."""
        yaml_content = """# Basin specification
- id: 1 # Main basin
  outlet:
    longitude: -76.0093
    latitude: 39.4620
  accumulation_threshold: 100000.0
  small_river_threshold: 10000.0
  flags:
    dam: 0
    debug: 1
    disconnected: 0
  paths:
    dam: ./dams.csv
    flowline_filter: ./flowline.geojson
    flowline_raw: ./allflowline.shp
    flowline_topo: ./flowline.csv
"""
        with tempfile.TemporaryDirectory() as tmpdir:
            basin_path = Path(tmpdir) / "basins.yaml"
            basin_path.write_text(yaml_content, encoding="utf-8")

            basins = BasinConfigManager.load_config(basin_path)
            self.assertEqual(len(basins), 1)
            b = basins[0]
            self.assertEqual(b["lBasinID"], 1)
            self.assertAlmostEqual(b["dLongitude_outlet_degree"], -76.0093)
            self.assertAlmostEqual(b["dLatitude_outlet_degree"], 39.4620)
            self.assertEqual(b["dAccumulation_threshold"], 100000.0)
            self.assertEqual(b["dThreshold_small_river"], 10000.0)
            self.assertEqual(b["iFlag_dam"], 0)
            self.assertEqual(b["iFlag_debug"], 1)
            self.assertEqual(b["sFilename_flowline_filter"], "./flowline.geojson")

    def test_inline_basins_in_case_yaml(self):
        """Test defining basins directly inline within the case YAML file."""
        yaml_content = """workflow:
  mesh_type: square
  case_index: 2
  date: "20240101"

resolutions:
  degree: 0.1
  meter: 10000.0

region:
  bbox:
    left: -77.0
    right: -76.0
    bottom: 40.0
    top: 41.0

paths:
  workspace_output: ./out

basins:
  - id: 1
    outlet:
      longitude: -76.5
      latitude: 40.5
    accumulation_threshold: 80000.0
    paths:
      flowline_filter: ./stream.geojson
  - id: 2
    outlet:
      longitude: -76.8
      latitude: 40.8
    accumulation_threshold: 90000.0
    paths:
      flowline_filter: ./stream2.geojson
"""
        with tempfile.TemporaryDirectory() as tmpdir:
            yaml_path = Path(tmpdir) / "case_with_inline_basins.yaml"
            yaml_path.write_text(yaml_content, encoding="utf-8")

            case = pyflowline_read_configuration_file(yaml_path)
            self.assertIsNotNone(case)
            self.assertEqual(case.sMesh_type, "square")
            self.assertEqual(case.nOutlet, 2)
            self.assertEqual(len(case.aBasin), 2)
            self.assertEqual(case.aBasin[0].lBasinID, 1)
            self.assertEqual(case.aBasin[1].lBasinID, 2)
            self.assertAlmostEqual(case.aBasin[0].dLongitude_outlet_degree, -76.5)
            self.assertAlmostEqual(case.aBasin[1].dLongitude_outlet_degree, -76.8)

    def test_extensionless_and_directory_discovery(self):
        """Test auto-discovery of .yaml / .yml / .json and directory configurations."""
        with tempfile.TemporaryDirectory() as tmpdir:
            case_dir = Path(tmpdir) / "my_case"
            case_dir.mkdir()
            config_file = case_dir / "workflow_config.yaml"
            config_file.write_text(
                "workflow:\n  mesh_type: latlon\nresolutions:\n  degree: 0.5\n",
                encoding="utf-8",
            )

            # Discover via directory path
            cfg1 = load_pyflowline_config(case_dir)
            self.assertEqual(cfg1["sMesh_type"], "latlon")

            # Discover via stem without extension
            stem_path = case_dir / "workflow_config"
            cfg2 = load_pyflowline_config(stem_path)
            self.assertEqual(cfg2["sMesh_type"], "latlon")

    def test_yaml_and_json_parity(self):
        """Test that equivalent YAML and JSON files compile to identical case configurations."""
        base_dir = Path(__file__).resolve().parent.parent / "data" / "susquehanna" / "input"
        mesh_types = ["latlon", "hexagon", "square", "mpas"]

        for mesh in mesh_types:
            yaml_path = base_dir / f"pyflowline_susquehanna_{mesh}.yaml"
            json_path = base_dir / f"pyflowline_susquehanna_{mesh}.json"

            self.assertTrue(yaml_path.is_file(), f"Missing {yaml_path}")
            self.assertTrue(json_path.is_file(), f"Missing {json_path}")

            cfg_yaml = load_pyflowline_config(yaml_path)
            cfg_json = load_pyflowline_config(json_path)

            # Verify core modeling keys match
            self.assertEqual(cfg_yaml["sMesh_type"], cfg_json["sMesh_type"])
            self.assertEqual(cfg_yaml["iCase_index"], cfg_json["iCase_index"])
            self.assertEqual(cfg_yaml["dResolution_meter"], cfg_json["dResolution_meter"])
            self.assertEqual(cfg_yaml["iFlag_create_mesh"], cfg_json["iFlag_create_mesh"])
            self.assertEqual(cfg_yaml["iFlag_simplification"], cfg_json["iFlag_simplification"])
            self.assertEqual(cfg_yaml["iFlag_intersect"], cfg_json["iFlag_intersect"])

    def test_change_config_key_value_yaml(self):
        """Test modifying keys in a YAML file using change_config_key_value."""
        yaml_content = """workflow:
  mesh_type: hexagon
  case_index: 1

resolutions:
  degree: 0.5
  meter: 50000.0

paths:
  workspace_output: ./old_output
"""
        with tempfile.TemporaryDirectory() as tmpdir:
            cfg_path = Path(tmpdir) / "test_modify.yaml"
            cfg_path.write_text(yaml_content, encoding="utf-8")

            change_config_key_value(cfg_path, "workspace_output", "./new_output")
            change_config_key_value(cfg_path, "case_index", 5)

            reloaded = load_pyflowline_config(cfg_path)
            self.assertEqual(reloaded["iCase_index"], 5)
            self.assertTrue(reloaded["sWorkspace_output"].endswith("new_output"))

    def test_case_and_basin_yaml_export(self):
        """Test exporting case and basin objects to YAML and reloading."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_case_yaml = Path(tmpdir) / "exported_case.yaml"
            out_basin_yaml = Path(tmpdir) / "exported_basin.yaml"

            # Create template case
            case = FlowlineConfigManager.create_template_config(
                out_case_yaml,
                custom_values={"sMesh_type": "hexagon", "iCase_index": 2},
            )
            self.assertIsNotNone(case)
            self.assertTrue(out_case_yaml.is_file())

            # Create template basins
            basins = BasinConfigManager.create_multiple_basins(out_basin_yaml, 2)
            self.assertEqual(len(basins), 2)
            self.assertTrue(out_basin_yaml.is_file())

            # Reload
            reloaded_case = load_pyflowline_config(out_case_yaml)
            self.assertEqual(reloaded_case["sMesh_type"], "hexagon")
            self.assertEqual(reloaded_case["iCase_index"], 2)

            reloaded_basins = BasinConfigManager.load_config(out_basin_yaml)
            self.assertEqual(len(reloaded_basins), 2)
            self.assertEqual(reloaded_basins[0]["lBasinID"], 1)
            self.assertEqual(reloaded_basins[1]["lBasinID"], 2)

            # Test pyflowline_export_config_to_yaml and pyflowline_export_basin_config_to_yaml on case object
            from pyflowline.classes.pycase import flowlinecase

            oCase = flowlinecase(reloaded_case, iFlag_create_directory_in=0)
            oCase.aBasin = basins

            case_export_file = Path(tmpdir) / "case_method_export.yaml"
            basin_export_file = Path(tmpdir) / "basin_method_export.yaml"

            oCase.pyflowline_export_config_to_yaml(case_export_file)
            self.assertTrue(case_export_file.is_file())

            oCase.pyflowline_export_basin_config_to_yaml(basin_export_file)
            self.assertTrue(basin_export_file.is_file())


if __name__ == "__main__":
    unittest.main()
