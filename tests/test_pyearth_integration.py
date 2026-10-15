"""Unit tests for pyearthmesh and pyearthriver integration in pyflowline."""

import unittest
from pyflowline.classes.pycase import flowlinecase
from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_cell
from pyflowline.formats.convert_attributes import convert_gcs_attributes_to_cell
from pyearthmesh.classes.hexagon import pyhexagon
from pyearthmesh.classes.square import pysquare
from pyearthmesh.classes.latlon import pylatlon
from pyearthmesh.classes.mpas import pympas
from pyearthmesh.classes.link import pycelllink
from pyearthmesh.meshes.structured.hexagon.create_hexagon_mesh import (
    create_hexagon_mesh,
)
from pyearthmesh.meshes.structured.square.create_square_mesh import (
    create_square_mesh,
)
from pyearthmesh.meshes.structured.latlon.create_latlon_mesh import (
    create_latlon_mesh,
)
from pyearthriver.classes.rivergraph import pyrivergraph
from pyearthriver.classes.flowline import pyflowline
from pyearthriver.classes.vertex import pyvertex
from pyearthriver.classes.edge import pyedge


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


if __name__ == "__main__":
    unittest.main()
