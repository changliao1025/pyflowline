from pathlib import Path
from pkg_resources import resource_filename

def pyflowline_project_root() -> Path:
    """
    Attempt to find the root path of the project dynamically.
    It first tries to find a 'setup.py' file. If not found, it uses pkg_resources
    to find the installation path.
    """
    try:
        return root_path_from_setup_file()
    except FileNotFoundError:
        return root_path_from_pyflowline_package_root()
    except Exception as e:
        raise RuntimeError(f"An unexpected error occurred while locating the project root: {str(e)}")

def pyflowline_package_root(package_name='pyflowline') -> Path:
    """Get the installation root directory of the package."""
    return Path(resource_filename(package_name, '')).resolve()

def root_path_from_setup_file() -> Path:
    """Return the top-level (root) path of the pyflowline project.
    This function navigates upwards from this file's path until it finds a directory
    with a specific marker (e.g., a setup.py file) indicating the root of the project.
    """
    this_path = Path(__file__).resolve()
    for parent in this_path.parents:
        if (parent / 'setup.py').exists():
            return parent
    raise FileNotFoundError("setup.py not found. Is this the right project structure?")

def root_path_from_pyflowline_package_root() -> Path:
    """Get the installation root directory of the package."""
    return pyflowline_package_root().parent

# Example additional function for getting root from this module's location, adjusted for best practice
def root_path_from_path_manager_location() -> Path:
    """Return the top-level (root) path of the pyflowline project
    This function assumes the module is located at a fixed level:
        root_path/pyflowline/configuration/path_manager
    """
    return Path(__file__).resolve().parents[2]
