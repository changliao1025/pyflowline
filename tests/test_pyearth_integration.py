"""Unit tests for pyearthmesh and pyearthriver integration in pyflowline."""

import json
import os
import tempfile
import unittest

from pyearthmesh.classes.hexagon import pyhexagon
from pyearthmesh.classes.latlon import pylatlon
from pyearthmesh.classes.link import pycelllink
from pyearthmesh.classes.mpas import pympas
from pyearthmesh.classes.square import pysquare
from pyearthmesh.meshes.structured.hexagon.create_hexagon_mesh import create_hexagon_mesh
from pyearthmesh.meshes.structured.latlon.create_latlon_mesh import create_latlon_mesh
from pyearthmesh.meshes.structured.square.create_square_mesh import create_square_mesh
from pyearthriver.classes.edge import pyedge
from pyearthriver.classes.flowline import pyflowline
from pyearthriver.classes.rivergraph import pyrivergraph
from pyearthriver.classes.vertex import pyvertex
from pyearthriver.formats.load_flowline import read_flowlines_from_vector

from pyflowline.classes.pycase import flowlinecase
from pyflowline.formats.convert_attributes import convert_gcs_attributes_to_cell
from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_cell
from pyflowline.formats.convert_flowline_to_geojson import convert_flowline_to_geojson
from pyflowline.formats.read_flowline import read_flowline_vector


class TestPyEarthIntegration(unittest.TestCase):
    """Test suite for integration with pyearthmesh and pyearthriver."""

    def test_imports(self):
        """Test that all required classes and generators are imported."""
        self.assertIsNotNone(pyhexagon)
        self.assertIsNotNone(pysquare)
        self.assertIsNotNone(pylatlon)
        self.assertIsNotNone(pympas)
        self.assertIsNotNone(pycelllink)
        self.assertIsNotNone(create_hexagon_mesh)
        self.assertIsNotNone(create_square_mesh)
        self.assertIsNotNone(create_latlon_mesh)
        self.assertIsNotNone(pyrivergraph)
        self.assertIsNotNone(pyflowline)
        self.assertIsNotNone(pyedge)
        self.assertIsNotNone(pyvertex)
        self.assertIsNotNone(flowlinecase)
        self.assertIsNotNone(convert_gcs_attributes_to_cell)
        self.assertIsNotNone(read_flowline_vector)

    def test_convert_coordinates_to_square_cell(self):
        """Test converting coordinates to square cell using pyearthmesh."""
        coords = [
            [-76.0, 40.0],
            [-75.0, 40.0],
            [-75.0, 41.0],
            [-76.0, 41.0],
            [-76.0, 40.0],
        ]
        cell = convert_gcs_coordinates_to_cell(2, -75.5, 40.5, coords)
        self.assertIsInstance(cell, pysquare)
        self.assertEqual(cell.dLongitude_center_degree, -75.5)
        self.assertEqual(cell.dLatitude_center_degree, 40.5)

    def test_convert_coordinates_to_latlon_cell(self):
        """Test converting coordinates to latlon cell using pyearthmesh."""
        coords = [
            [-76.0, 40.0],
            [-75.0, 40.0],
            [-75.0, 41.0],
            [-76.0, 41.0],
            [-76.0, 40.0],
        ]
        cell = convert_gcs_coordinates_to_cell(3, -75.5, 40.5, coords)
        self.assertIsInstance(cell, pylatlon)
        self.assertEqual(cell.dLongitude_center_degree, -75.5)
        self.assertEqual(cell.dLatitude_center_degree, 40.5)

    def test_read_flowline_vector_geojson(self):
        """Test reading flowlines from GeoJSON using unified vector reader."""
        geojson_data = {
            "type": "FeatureCollection",
            "features": [
                {
                    "type": "Feature",
                    "geometry": {
                        "type": "LineString",
                        "coordinates": [
                            [-76.5, 41.0],
                            [-76.4, 40.9],
                            [-76.3, 40.8],
                        ],
                    },
                    "properties": {
                        "lineid": 101,
                        "stream_order": 2,
                        "stream_segment": 5,
                        "drainage_area": 1250000.0,
                    },
                },
                {
                    "type": "Feature",
                    "geometry": {
                        "type": "LineString",
                        "coordinates": [
                            [-76.3, 40.8],
                            [-76.2, 40.7],
                            [-76.1, 40.6],
                        ],
                    },
                    "properties": {
                        "lineid": 102,
                        "stream_order": 3,
                        "stream_segment": 6,
                        "drainage_area": 2500000.0,
                    },
                },
            ],
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            geojson_path = os.path.join(tmpdir, "test_network.geojson")
            with open(geojson_path, "w", encoding="utf-8") as f:
                json.dump(geojson_data, f)

            flowlines, proj_wkt = read_flowline_vector(geojson_path)
            self.assertEqual(len(flowlines), 2)
            self.assertEqual(flowlines[0].lFlowlineID, 101)
            self.assertEqual(flowlines[0].iStream_order, 2)
            self.assertEqual(flowlines[0].iStream_segment, 5)
            self.assertEqual(flowlines[0].dDrainage_area, 1250000.0)
            self.assertEqual(flowlines[1].lFlowlineID, 102)
            self.assertEqual(flowlines[1].iStream_order, 3)

    def test_read_flowlines_from_dict(self):
        """Test reading flowlines directly from in-memory GeoJSON dictionary."""
        geojson_data = {
            "type": "FeatureCollection",
            "features": [
                {
                    "type": "Feature",
                    "geometry": {
                        "type": "LineString",
                        "coordinates": [[-76.0, 40.0], [-76.1, 40.1]],
                    },
                    "properties": {"lineid": 1},
                }
            ],
        }
        flowlines = read_flowlines_from_vector(geojson_data)
        self.assertEqual(len(flowlines), 1)
        self.assertEqual(flowlines[0].lFlowlineID, 1)

    def test_convert_flowline_to_geojson(self):
        """Test converting a vector flowline dataset to GeoJSON."""
        geojson_data = {
            "type": "FeatureCollection",
            "features": [
                {
                    "type": "Feature",
                    "geometry": {
                        "type": "LineString",
                        "coordinates": [[-76.5, 41.0], [-76.4, 40.9]],
                    },
                    "properties": {
                        "lineid": 42,
                        "stream_order": 2,
                        "drainage_area": 50000.0,
                    },
                }
            ],
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, "input.geojson")
            output_path = os.path.join(tmpdir, "output.geojson")
            with open(input_path, "w", encoding="utf-8") as f:
                json.dump(geojson_data, f)

            flowlines = convert_flowline_to_geojson(input_path, output_path)
            self.assertEqual(len(flowlines), 1)
            self.assertEqual(flowlines[0].lFlowlineID, 42)
            self.assertTrue(os.path.isfile(output_path))

            # Verify exported content can be read back
            reloaded, _ = read_flowline_vector(output_path)
            self.assertEqual(len(reloaded), 1)
            self.assertEqual(reloaded[0].lFlowlineID, 42)
            self.assertEqual(reloaded[0].iStream_order, 2)
            self.assertEqual(reloaded[0].dDrainage_area, 50000.0)

    def test_convert_flowline_to_geojson_file_not_found(self):
        """Test FileNotFoundError when input file does not exist."""
        with self.assertRaises(FileNotFoundError):
            convert_flowline_to_geojson(
                "/path/to/non_existent_flowline.geojson", "/path/to/out.geojson"
            )

    def test_find_nearest_outlet(self):
        """Test snapping approximate coordinates to the nearest outlet vertex."""
        from pyearthriver.algorithms.topology.connectivity import find_nearest_outlet

        v1 = pyvertex({"dLongitude_degree": -77.0, "dLatitude_degree": 40.0})
        v2 = pyvertex({"dLongitude_degree": -76.5, "dLatitude_degree": 40.0})
        e = pyedge(v1, v2)
        fl = pyflowline([e])

        target = pyvertex({"dLongitude_degree": -76.49, "dLatitude_degree": 40.01})
        snapped = find_nearest_outlet([fl], target)
        self.assertEqual(snapped, v2)

        rg = pyrivergraph([fl])
        snapped_rg = rg.find_nearest_outlet(target)
        self.assertEqual(snapped_rg, v2)


if __name__ == "__main__":
    unittest.main()
