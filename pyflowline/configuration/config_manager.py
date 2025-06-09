import os
from pathlib import Path
#use this function to generate an initial json file for hexwatershed
import json
#once it's generated, you can modify it and use it for different simulations
from pyflowline.classes.pycase import flowlinecase
from pyflowline.classes.basin import pybasin

class BasinConfigManager:
    """Configuration manager for basin objects that handles defaults and serialization"""

    @staticmethod
    def get_default_config():
        """Returns a dictionary with all default basin configuration values"""
        return {
            "iFlag_dam": 0,
            "iFlag_disconnected": 0,
            "lBasinID": 1,
            "dLatitude_outlet_degree": -180,
            "dLongitude_outlet_degree": 180,
            "dAccumulation_threshold": -90,
            "dThreshold_small_river": 90,
            "sFilename_dam": "ICoM_dams.csv",
            "sFilename_flowline_filter": "streamord7above.shp",
            "sFilename_flowline_raw": "allflowline.shp",
            "sFilename_flowline_topo": "flowline.csv",
            "sWorkspace_output_basin": ""
        }

    @staticmethod
    def create_template_config(basin_id=1, custom_values=None):
        """Create a configuration dictionary for a basin with default values, optionally customized

        Args:
            basin_id (int): Basin identifier
            custom_values (dict, optional): Dictionary of values to override defaults

        Returns:
            dict: The created configuration
        """
        from pathlib import Path

        # Get defaults
        config = BasinConfigManager.get_default_config()

        # Update basin ID
        config["lBasinID"] = basin_id

        # Handle workspace paths
        sWorkspace_input = Path.cwd()


        config["sFilename_dam"] = str(sWorkspace_input / config["sFilename_dam"])
        config["sFilename_flowline_filter"] = str(sWorkspace_input / config["sFilename_flowline_filter"])
        config["sFilename_flowline_raw"] = str(sWorkspace_input / config["sFilename_flowline_raw"])
        config["sFilename_flowline_topo"] = str(sWorkspace_input / config["sFilename_flowline_topo"])

        sWorkspace_output = Path.cwd()
        basin_folder = "{:09d}".format(basin_id)
        config["sWorkspace_output_basin"] = str(sWorkspace_output / basin_folder)

        # Apply customizations
        if custom_values:
            for key, value in custom_values.items():
                if isinstance(value, Path):
                    value = str(value)
                config[key] = value

        return config

    @staticmethod
    def create_multiple_basins(sFilename_basins_json, nBasin, custom_values=None):
        """Generate configuration for multiple basins and save to JSON file

        Args:
            sFilename_basins_json (str or Path): Output filename for basin configurations
            nBasin (int): Number of basins to create
            sWorkspace_input (str or Path): Input workspace path
            sWorkspace_output (str or Path): Output workspace path
            custom_values (dict, optional): Dictionary with custom values to apply to all basins

        Returns:
            list: List of basin objects
        """
        from pyflowline.classes.basin import pybasin
        import json
        import os
        from pathlib import Path

        # Ensure output directory exists
        os.makedirs(os.path.dirname(os.path.abspath(sFilename_basins_json)), exist_ok=True)

        aBasin_out = []
        for i in range(nBasin):
            # Create basin config with basin-specific ID
            aConfig_basin = BasinConfigManager.create_template_config(
                basin_id=i+1,
                custom_values=custom_values
            )

            # Create basin object
            pBasin = pybasin(aConfig_basin, iFlag_create_directory_in = 0)
            aBasin_out.append(pBasin)

        # Export basin config to a file
        with open(sFilename_basins_json, 'w', encoding='utf-8') as f:
            sJson = json.dumps([json.loads(ob.tojson()) for ob in aBasin_out], indent=4)
            f.write(sJson)

        return aBasin_out

class FlowlineConfigManager:
    """Configuration manager for pyflowline that handles defaults and serialization"""

    @staticmethod
    def get_default_config():
        """Returns a dictionary with all default configuration values"""
        return {
            "iFlag_use_mesh_dem": 0,
            "iFlag_save_mesh": 1,
            "iFlag_simplification": 1,
            "iFlag_create_mesh": 1,
            "iFlag_intersect": 1,
            "iFlag_resample_method": 1,
            "iFlag_global": 0,
            "iFlag_multiple_outlet": 0,
            "iFlag_elevation_profile": 1,
            "iFlag_rotation": 0,
            "iFlag_stream_burning_topology": 1,
            "iFlag_save_elevation": 1,
            "nOutlet": 1,
            "dResolution_degree": 0.5,
            "dResolution_meter": 50000,
            "dLongitude_left": -180,
            "dLongitude_right": 180,
            "dLatitude_bot": -90,
            "dLatitude_top": 90,
            "sRegion": "susquehanna",
            "sModel": "pyflowline",
            "iCase_index": 1,
            "sMesh_type": "hexagon",
            "sJob": "pyflowline",
            "sDate": "20220202",
            "flowline_info": "flowline_info.json",
            "sFilename_mesh_info": "mesh_info.json",
            "sFilename_elevation": "elevation.json"
        }

    @staticmethod
    def create_template_config(output_filename, custom_values=None):
        """Create a configuration file with default values, optionally customized

        Args:
            output_filename (str or Path): Path to save the configuration file
            custom_values (dict, optional): Dictionary of values to override defaults

        Returns:
            dict: The created configuration
        """
        # Get defaults
        config = FlowlineConfigManager.get_default_config()

        # Handle workspace paths
        config["sWorkspace_input"] = str(Path.cwd())
        config["sWorkspace_output"] = str(Path.cwd())
        config["sFilename_model_configuration"] = str(output_filename)

        # Apply customizations
        if custom_values:
            for key, value in custom_values.items():
                if isinstance(value, Path):
                    value = str(value)
                config[key] = value

        # Ensure directory exists
        os.makedirs(os.path.dirname(os.path.abspath(output_filename)), exist_ok=True)

        # Write to file
        with open(output_filename, 'w') as f:
            json.dump(config, f, indent=4)

        return config

class JigsawConfigManager:
    """Configuration manager for JIGSAW that handles defaults and serialization"""

    @staticmethod
    def get_default_config():
        """Returns a dictionary with all default JIGSAW configuration values"""
        return {
            # Grid resolution parameters
            "ncolumn_space": 360,        # Number of columns in spacing grid (longitude)
            "nrow_space": 180,           # Number of rows in spacing grid (latitude)
            "dSpac_value": 50.0,         # Default spacing value

            # Feature flags for geometry generation
            "iFlag_geom": False,         # Enable geometry generation
            "iFlag_spac": False,         # Enable spacing function generation
            "iFlag_init": False,         # Enable initialization mesh
            "iFlag_opts": False,         # Enable custom options

            # Environment type flags
            "iFlag_ocean": False,        # Apply ocean-specific spacing
            "iFlag_land": False,         # Apply land-specific spacing

            # Point feature geometry and spacing flags
            "iFlag_geom_dam": False,     # Include dam geometries
            "iFlag_spac_dam": False,     # Apply dam-specific spacing
            "iFlag_geom_city": False,    # Include city geometries
            "iFlag_spac_city": False,    # Apply city-specific spacing

            # Line feature geometry and spacing flags
            "iFlag_geom_river_network": False,   # Include river network geometries
            "iFlag_spac_river_network": False,   # Apply river network-specific spacing
            "iFlag_geom_coastline": False,       # Include coastline geometries
            "iFlag_spac_coastline": False,       # Apply coastline-specific spacing

            # Polygon feature geometry and spacing flags
            "iFlag_geom_watershed_boundary": False,  # Include watershed boundary geometries
            "iFlag_spac_watershed_boundary": False,  # Apply watershed boundary-specific spacing
            "iFlag_geom_lake_boundary": False,       # Include lake boundary geometries
            "iFlag_spac_lake_boundary": False,       # Apply lake boundary-specific spacing

            # Resolution parameters for different features (in degrees)
            "dResolution_land": 45.0,                # Land feature spacing
            "dResolution_dam": 3.0,                  # Dam feature spacing
            "dResolution_city": 3.0,                 # City feature spacing
            "dResolution_river_network": 3.0,        # River network feature spacing
            "dResolution_coastline": 3.0,            # Coastline feature spacing
            "dResolution_watershed_boundary": 3.0,    # Watershed boundary feature spacing
            "dResolution_lake_boundary": 3.0,        # Lake boundary feature spacing

            # Mesh type identifiers
            "geom_mshID": "ellipsoid-mesh",          # Geometry mesh type
            "spac_mshID": "ellipsoid-grid",          # Spacing grid type

            # Earth/sphere parameters
            "FULL_SPHERE_RADIUS": 6371.0,            # Earth radius in km

            # Gradient limiting
            "dhdx_lim": 0.25,                        # Gradient limit for mesh sizing

            # File paths for features
            "sFilename_dam_vector": None,            # Dam vector file
            "sFilename_dam_raster": None,            # Dam raster file
            "sFilename_city_vector": None,           # City vector file
            "sFilename_city_raster": None,           # City raster file
            "sFilename_river_network_vector": None,  # River network vector file
            "sFilename_river_network_raster": None,  # River network raster file
            "sFilename_coastline_vector": None,      # Coastline vector file
            "sFilename_coastline_raster": None,      # Coastline raster file
            "sFilename_watershed_boundary_vector": None,  # Watershed boundary vector file
            "sFilename_watershed_boundary_raster": None,  # Watershed boundary raster file
            "sFilename_lake_boundary_vector": None,  # Lake boundary vector file
            "sFilename_lake_boundary_raster": None,  # Lake boundary raster file

            # Mesh sizing parameters
            "hfun_hmax": float("inf"),     # Max. refinement function value
            "hfun_hmin": 0.0,              # Min. refinement function value
            "hfun_scal": "absolute",       # Scaling type: "relative" or "absolute"
            "mesh_dims": 2,                # Mesh dimension (2 for surface)
            "bisection": -1,               # Bisection method (-1 for heuristic)

            # Optimization parameters
            "optm_qlim": 0.95,             # Quality limit for optimization
            "optm_iter": 32,               # Number of optimization iterations
            "optm_qtol": 1.0E-05,          # Quality tolerance

            # Core mesh sizing and quality parameters from JigsawConfigManager
            "mesh_rad2": 1.5,            # Max. radius-edge ratio
            "mesh_rad3": 2.0,            # Max. radius-circumsphere ratio for tetras
            "mesh_eps1": 0.333,          # Min. mesh quality threshold
            "mesh_eps2": 0.333,          # Min. mesh quality threshold for tetra
            "mesh_top": 1,               # Mesh topology (1 for manifold surface)
            "mesh_iter": 3,              # Mesh iteration limit

            # Verbosity and iterations
            "verbosity": 0,              # Verbosity level (0-3)

            # File paths (will be populated based on workspace)
            "geom_file": None,           # Input geometry file
            "hfun_file": None,           # Input mesh-size file
            "mesh_file": None,           # Output mesh file

            # Algorithm selection
            "mesh_kern": "delfront",     # Meshing kernel: "delfront" or "delaunay"
            "optm_kern": "odt+dqdx",     # Optimisation kernel

            # Region boundary
            "geom_feat": True,           # Detect sharp features in geometry

            # Output options
            "mesh_type": "euclidean-mesh",  # Mesh type (euclidean-mesh or ellipsoid-mesh)
            "output_formats": ["vtk", "gmsh"]  # Output formats to generate
        }

    @staticmethod
    def create_template_config(output_filename,  custom_values=None):
        """Create a JIGSAW configuration file with default values, optionally customized

        Args:
            output_filename (str or Path): Path to save the configuration file
            custom_values (dict, optional): Dictionary of values to override defaults

        Returns:
            dict: The created configuration
        """
        # Get defaults
        config = JigsawConfigManager.get_default_config()


        # Apply customizations
        if custom_values:
            for key, value in custom_values.items():
                if isinstance(value, Path):
                    value = str(value)
                config[key] = value

        # Ensure directory exists
        os.makedirs(os.path.dirname(os.path.abspath(output_filename)), exist_ok=True)

        # Write to file
        with open(output_filename, 'w') as f:
            json.dump(config, f, indent=4)

        return config

    @staticmethod
    def load_config(filename):
        """Load a JIGSAW configuration from a JSON file

        Args:
            filename (str or Path): Path to the configuration file

        Returns:
            dict: The loaded configuration
        """
        with open(filename, 'r') as f:
            return json.load(f)

# Update the original create_template_basin_configuration_file function to use the new manager
def create_template_basin_configuration_file(
        sFilename_basins_json,
        nBasin):
    """Generate basin configuration using the BasinConfigManager

    Args:
        sFilename_basins_json (str or Path): the filename
        nBasin (int): the total number of basin
        sWorkspace_input_in (str or Path): the input data path
        sWorkspace_output (str or Path): the output path

    Returns:
        list: a list of basin objects
    """
    return BasinConfigManager.create_multiple_basins(
        sFilename_basins_json,
        nBasin
    )
# Then you can simplify your create_template_configuration_file function
def create_template_configuration_file(sFilename_json, **kwargs):
    """Generate pyflowline config template file using parameter keywords"""

    # Create the configuration
    config = FlowlineConfigManager.create_template_config(sFilename_json, kwargs)

    # Initialize the model with the config
    oModel = flowlinecase(config, iFlag_create_directory_in = 0)

    # Generate basin config
    sDirname = os.path.dirname(sFilename_json)
    sFilename = Path(sFilename_json).stem + '_basins.json'
    sFilename_basins_json = os.path.join(sDirname, sFilename)

    aBasin = create_template_basin_configuration_file(
        sFilename_basins_json,
        config.get("nOutlet", 1))

    oModel.aBasin = aBasin
    oModel.sFilename_basins = sFilename_basins_json
    oModel.pyflowline_export_config_to_json(sFilename_json, iFlag_export_basin_in=0)

    return oModel

def create_template_jigsaw_configuration_file(sFilename_json, **kwargs):
    """Generate JIGSAW config template file using parameter keywords"""

    # Create the configuration
    config = JigsawConfigManager.create_template_config(sFilename_json, custom_values=kwargs)

    # Ensure the output directory exists
    os.makedirs(os.path.dirname(sFilename_json), exist_ok=True)

    # Write the configuration to the specified JSON file
    with open(sFilename_json, 'w') as f:
        json.dump(config, f, indent=4)

    return config


