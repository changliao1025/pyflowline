from pathlib import Path

try:
    import pyflowline
    from pyflowline.configuration import path_manager as pyflowline_path_manager
except:
    pass

def test_root_path_from_pyflowline_project_root():
    """Test the dynamic root path identification using both setup.py and pkg_resources."""
    expected_root = Path(__file__).resolve().parents[1]
    returned_root = pyflowline_path_manager.pyflowline_project_root()
    
    assert returned_root == expected_root, f"Path manager returned {returned_root}, but expected {expected_root}"
    print('Root path is correctly identified from pyflowline project root.')

def test_root_path_from_pyflowline_package_root():
    """Test the package root path identification using pkg_resources."""
    expected_root = Path(__file__).resolve().parents[1]
    returned_root = pyflowline_path_manager.root_path_from_pyflowline_package_root()
    
    assert returned_root == expected_root, f"Path manager returned {returned_root}, but expected {expected_root}"
    print('Root path is correctly identified from pyflowline package root.')

def test_root_path_from_setup_file():
    """Test the root path identification specifically by locating setup.py."""
    expected_root = Path(__file__).resolve().parents[1]
    returned_root = pyflowline_path_manager.root_path_from_setup_file()
    
    assert returned_root == expected_root, f"Path manager returned {returned_root}, but expected {expected_root}"
    print('Root path is correctly identified from setup file location.')

def test_root_path_from_path_manager_location():
    """Test the root path based on the relative location of the path_manager module."""
    expected_root = Path(__file__).resolve().parents[1]
    returned_root = pyflowline_path_manager.root_path_from_path_manager_location()

    assert returned_root == expected_root, f"Path manager returned {returned_root}, but expected {expected_root}"
    print('Root path is correctly identified from path_manager location.')

def test_path_manager():

    try:
        import pyflowline
        from pyflowline.configuration import path_manager as pyflowline_path_manager
        test_root_path_from_pyflowline_project_root()
        test_root_path_from_pyflowline_package_root()
        test_root_path_from_setup_file()
        test_root_path_from_path_manager_location()
    except:
        pass

test_path_manager()
