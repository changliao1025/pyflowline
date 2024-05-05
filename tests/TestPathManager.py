import unittest
from pathlib import Path
from pyflowline.configuration import path_manager

class TestPathManager(unittest.TestCase):

    def test_root_path_from_pyflowline_project_root(self):
        """Test the dynamic root path identification using both setup.py and pkg_resources."""
        expected_root = Path(__file__).resolve().parents[1]
        returned_root = path_manager.pyflowline_project_root()
        self.assertEqual(returned_root, expected_root, "Root path is incorrectly identified from pyflowline project root.")

    def test_root_path_from_pyflowline_package_root(self):
        """Test the package root path identification using pkg_resources."""
        expected_root = Path(__file__).resolve().parents[1]
        returned_root = path_manager.root_path_from_pyflowline_package_root()
        self.assertEqual(returned_root, expected_root, "Root path is incorrectly identified from pyflowline package root.")

    def test_root_path_from_setup_file(self):
        """Test the root path identification specifically by locating setup.py."""
        expected_root = Path(__file__).resolve().parents[1]
        returned_root = path_manager.root_path_from_setup_file()
        self.assertEqual(returned_root, expected_root, "Root path is incorrectly identified from setup file location.")

    def test_root_path_from_path_manager_location(self):
        """Test the root path based on the relative location of the path_manager module."""
        expected_root = Path(__file__).resolve().parents[1]
        returned_root = path_manager.root_path_from_path_manager_location()
        self.assertEqual(returned_root, expected_root, "Root path is incorrectly identified from path_manager location.")

if __name__ == '__main__':
    unittest.main()
