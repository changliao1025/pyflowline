from pyflowline.mesh.cubicsphere.create_cubicsphere_mesh import create_cubicsphere_mesh
from pyflowline.mesh.triangular.create_triangular_mesh import create_triangular_mesh
from pyflowline.mesh.tin.create_tin_mesh import create_tin_mesh
from pyflowline.mesh.dggrid.create_dggrid_mesh import dggrid_find_resolution_by_index
from pyflowline.mesh.dggrid.create_dggrid_mesh import create_dggrid_mesh
from pyflowline.mesh.mpas.create_mpas_mesh import create_mpas_mesh
from pyflowline.mesh.square.create_square_mesh import create_square_mesh
from pyflowline.mesh.latlon.create_latlon_mesh import create_latlon_mesh
from pyflowline.mesh.hexagon.create_hexagon_mesh import create_hexagon_mesh
import os
import shutil
import stat
from pathlib import Path
import json
from json import JSONEncoder
import datetime
import importlib.util
from shutil import copy2
import numpy as np
from osgeo import ogr, osr, gdal

from pyearth.system.define_global_variables import *
from pyearth.gis.gdal.gdal_check_file_type import gdal_check_file_type
from pyearth.gis.spatialref.convert_between_degree_and_meter import degree_to_meter
from pyearth.gis.spatialref.convert_between_degree_and_meter import meter_to_degree
from pyearth.gis.gdal.read.vector.gdal_read_geojson_boundary import gdal_read_geojson_boundary
from pyearth.toolbox.data.geoparquet.convert_geojson_to_geoparquet import convert_geojson_to_geoparquet
from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import gdal_read_geotiff_file
from pyearth.gis.gdal.read.raster.gdal_get_raster_extent import gdal_get_raster_extent
from pyearth.gis.spatialref.get_utm_spatial_reference import get_utm_spatial_reference_wkt
from pyearth.gis.spatialref.reproject_coordinates import reproject_coordinates
from pyearth.gis.gdal.read.vector.gdal_get_vector_spatial_reference import gdal_get_vector_spatial_reference_wkt
from pyearth.gis.gdal.read.raster.gdal_get_raster_spatial_reference import gdal_get_raster_spatial_reference_wkt

from pyflowline.classes.timer import pytimer
from pyflowline.classes.mpas import pympas
from pyflowline.classes.hexagon import pyhexagon
from pyflowline.classes.latlon import pylatlon
from pyflowline.classes.square import pysquare
from pyflowline.classes.dggrid import pydggrid
from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.basin import pybasin
from pyflowline.classes.flowline import pyflowline
from pyflowline.classes.edge import pyedge
from pyflowline.formats.read_mesh import read_mesh_json, read_mesh_json_w_topology
from pyflowline.formats.convert_boundary_to_geojson import convert_boundary_to_geojson

iFlag_kml = importlib.util.find_spec("simplekml")

iFlag_cython = importlib.util.find_spec("cython")
if iFlag_cython is not None:
    from tinyr import RTree
    iFlag_use_rtree = 1
else:
    iFlag_use_rtree = 0
    pass


gdal.UseExceptions()
pDate = datetime.datetime.today()
sDate_default = "{:04d}".format(
    pDate.year) + "{:02d}".format(pDate.month) + "{:02d}".format(pDate.day)


class CaseClassEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.float32):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, list):
            pass
        if isinstance(obj, pyvertex):
            return json.loads(obj.tojson())
        if isinstance(obj, pyedge):
            return obj.lEdgeID
        if isinstance(obj, pyflowline):
            return obj.lFlowlineID
        if isinstance(obj, pyhexagon):
            return obj.lCellID
        if isinstance(obj, pysquare):
            return obj.lCellID
        if isinstance(obj, pylatlon):
            return obj.lCellID
        if isinstance(obj, pympas):
            return obj.lCellID
        if isinstance(obj, pydggrid):
            return obj.lCellID
        if isinstance(obj, pybasin):
            return obj.lBasinID

        return JSONEncoder.default(self, obj)


class flowlinecase(object):
    """
    The flowline case class

    Args:
        object (obj): None

    Returns:
        flowlinecase: A flowlinecase object
    """
    iCase_index = 0

    sMesh_type = 1
    iFlag_standalone = 1
    iFlag_flowline = 1
    iFlag_watershed_boundary = None
    iFlag_global = 0
    iFlag_antarctic = 0
    iFlag_arctic = 0
    iFlag_multiple_outlet = 0
    iFlag_mesh_boundary = 0
    iFlag_land_ocean_mask = 0
    iFlag_run_jigsaw = 0
    iFlag_user_provided_binary = 0
    # iFlag_use_shapefile_extent=1
    iFlag_use_mesh_dem = 0
    iFlag_save_mesh = 0
    iFlag_simplification = 1  # user can turn on/off
    iFlag_create_mesh = 1
    iFlag_intersect = 1
    iFlag_rotation = 0
    iFlag_break_by_distance = 0  # if the distance between two vertice are far,
    iFlag_analysis = 0
    iFlag_evaluation = 0  # run evaluation or not
    iFlag_fill_hole = 0  # fill holes in the mesh
    iResolution_index = 10
    iFlag_run_dggrid = 0
    nOutlet = 1  # by default , there shoule ne only one ouelet
    dResolution_degree = 0.0
    dResolution_meter = 0.0
    dThreshold_small_river = 0.0
    dThreshold_break_by_distance = 5000.0
    dLongitude_left = -180
    dLongitude_right = 180
    dLatitude_bot = -90
    dLatitude_top = 90

    dLongitude_mean = 0.0  # the average location of the domain
    dLatitude_mean = 0.0

    dX_lowerleft = 0.0  # these coordinate only used for projection based meshes, and should be set during the parameter check function
    dY_lowerleft = 0.0

    dX_upperleft = 0.0
    dY_upperleft = 0.0

    dX_lowerright = 0.0
    dY_lowerright = 0.0

    dElevation_mean = -9999.0
    sFilename_model_configuration = ''

    sFilename_watershed_boundary = ''
    sFilename_watershed_boundary_geojson = ''
    sFilename_mesh_boundary = ''
    sFilename_mesh_boundary_geojson = ''
    sFilename_coastal_boundary = ''
    sFilename_coastal_boundary_geojson = ''
    sWorkspace_input = ''
    sWorkspace_output = ''
    sWorkspace_jigsaw = ''
    sWorkspace_dggrid = ''
    # sWorkspace_output_case=''
    sRegion = ''
    sModel = ''
    sMesh_type = 'hexagon'
    sDggrid_type = 'ISEA3H'

    sCase = ''
    sDate = ''

    sFilename_dggrid = ''

    sFilename_spatial_reference = ''
    sFilename_dem = ''
    # target mesh desired projection, can be any projection for local simulation
    pProjection_reference = ''
    # pyflowline internal projection, alwasy wgs84, all data will be re-projected to this projection
    pProjection_model = ''

    sFilename_mesh = ''
    sFilename_mesh_info = ''
    sFilename_mpas_mesh_netcdf = ''
    sFilename_jigsaw_mesh_netcdf = ''
    sFilename_mesh_kml = ''
    sFilename_mesh_parquet = ''

    sFilename_basins = ''
    sFilename_jigsaw_configuration = ''
    sFilename_dggrid_configuration = ''
    aConfig_jigsaw = None

    aBasin = list()
    aFlowline_simplified = list()
    aFlowline_conceptual = list()
    aCellID_outlet = list()
    aCell = list()

    pRTree_mesh = None  # only when rtree is used

    iFlag_visual = importlib.util.find_spec("cartopy")
    if iFlag_visual is not None:
        from ._visual import plot
        from ._visual import _plot_study_area
        from ._visual import _plot_mesh
        from ._visual import _plot_mesh_with_flowline
        from ._visual import _compare_with_raster_dem_method
    else:
        pass

    from ._hpc import _pyflowline_create_hpc_job

    def __init__(self, aConfig_in,
                 iFlag_standalone_in=None,
                 iFlag_create_directory_in = 1,
                 sModel_in=None,
                 sDate_in=None,
                 sWorkspace_output_in=None):
        """
        Initialize a flowlinecase object

        Args:
            aConfig_in (dict): A dictionary of parameters
            iFlag_standalone_in (int, optional): Flag for whether run the case standalone. Defaults to None.
            sModel_in (str, optional): The model name. Defaults to None.
            sDate_in (str, optional): The case date. Defaults to None.
            sWorkspace_output_in (str, optional): The output workspace. Defaults to None.
        """
        # flags
        if iFlag_standalone_in is not None:
            self.iFlag_standalone = iFlag_standalone_in
        else:
            if 'iFlag_standalone' in aConfig_in:
                self.iFlag_standalone = int(aConfig_in['iFlag_standalone'])
            else:
                self.iFlag_standalone = 1

        if 'iFlag_flowline' in aConfig_in:
            self.iFlag_flowline = int(aConfig_in['iFlag_flowline'])
        else:
            # without iFlag_flowline the model is a mesh generator
            self.iFlag_flowline = 1

        if 'iFlag_global' in aConfig_in:
            self.iFlag_global = int(aConfig_in['iFlag_global'])

        if 'iFlag_analysis' in aConfig_in:
            self.iFlag_analysis = int(aConfig_in['iFlag_analysis'])

        if 'iFlag_evaluation' in aConfig_in:
            self.iFlag_evaluation = int(aConfig_in['iFlag_evaluation'])

        if 'iFlag_antarctic' in aConfig_in:
            self.iFlag_antarctic = int(aConfig_in['iFlag_antarctic'])

        if 'iFlag_arctic' in aConfig_in:
            self.iFlag_arctic = int(aConfig_in['iFlag_arctic'])

        if 'iFlag_mesh_boundary' in aConfig_in:
            self.iFlag_mesh_boundary = int(aConfig_in['iFlag_mesh_boundary'])
        else:
            self.iFlag_mesh_boundary = 0

        if 'iFlag_multiple_outlet' in aConfig_in:
            self.iFlag_multiple_outlet = int(
                aConfig_in['iFlag_multiple_outlet'])

        if 'iFlag_simplification' in aConfig_in:
            self.iFlag_simplification = int(aConfig_in['iFlag_simplification'])
            # if simplification is desired, then we must activate flowline flag
            # it can be turn off only when a previous simplification was done and you dont want to redo it
        else:
            self.iFlag_simplification = 1
        if self.iFlag_simplification == 1:
            self.iFlag_flowline = 1
        else:
            self.iFlag_flowline = 0

        if 'iFlag_create_mesh' in aConfig_in:
            self.iFlag_create_mesh = int(aConfig_in['iFlag_create_mesh'])
        else:
            self.iFlag_create_mesh = 1

        if 'iFlag_save_mesh' in aConfig_in:
            self.iFlag_save_mesh = int(aConfig_in['iFlag_save_mesh'])

        if 'iFlag_run_jigsaw' in aConfig_in:
            self.iFlag_run_jigsaw = int(aConfig_in['iFlag_run_jigsaw'])

        if 'iFlag_use_mesh_dem' in aConfig_in:
            self.iFlag_use_mesh_dem = int(aConfig_in['iFlag_use_mesh_dem'])

        # if 'iFlag_use_shapefile_extent' in aConfig_in:
        #    self.iFlag_use_shapefile_extent = int(aConfig_in['iFlag_use_shapefile_extent'])

        if 'iFlag_user_provided_binary' in aConfig_in:
            self.iFlag_user_provided_binary = int(
                aConfig_in['iFlag_user_provided_binary'])
        else:
            self.iFlag_user_provided_binary = 0

        if 'iFlag_rotation' in aConfig_in:
            self.iFlag_rotation = int(aConfig_in['iFlag_rotation'])

        if 'iFlag_intersect' in aConfig_in:
            self.iFlag_intersect = int(aConfig_in['iFlag_intersect'])
        else:
            self.iFlag_intersect = 1

        if self.iFlag_flowline == 1:
            pass
        else:
            self.iFlag_intersect = 0

        if 'iFlag_break_by_distance' in aConfig_in:
            self.iFlag_break_by_distance = int(
                aConfig_in['iFlag_break_by_distance'])
        else:
            self.iFlag_break_by_distance = 0

        if 'iFlag_run_dggrid' in aConfig_in:
            self.iFlag_run_dggrid = int(aConfig_in['iFlag_run_dggrid'])
        else:
            self.iFlag_run_dggrid = 0

        if 'iFlag_fill_hole' in aConfig_in:
            self.iFlag_fill_hole = int(aConfig_in['iFlag_fill_hole'])
        else:
            self.iFlag_fill_hole = 0

        if 'iResolution_index' in aConfig_in:
            self.iResolution_index = int(aConfig_in['iResolution_index'])
        else:
            self.iResolution_index = 10

        if 'sDggrid_type' in aConfig_in:
            self.sDggrid_type = aConfig_in['sDggrid_type']
        else:
            self.sDggrid_type = 'ISEA3H'

        if 'nOutlet' in aConfig_in:
            self.nOutlet = int(aConfig_in['nOutlet'])

        if 'iCase_index' in aConfig_in:
            iCase_index = int(aConfig_in['iCase_index'])
        else:
            iCase_index = 1

        sCase_index = "{:03d}".format(iCase_index)
        self.iCase_index = iCase_index

        if 'dResolution_degree' in aConfig_in:
            self.dResolution_degree = float(aConfig_in['dResolution_degree'])

        if 'dResolution_meter' in aConfig_in:
            self.dResolution_meter = float(aConfig_in['dResolution_meter'])

        if 'dThreshold_break_by_distance' in aConfig_in:
            self.dThreshold_break_by_distance = float(
                aConfig_in['dThreshold_break_by_distance'])

        if 'dLongitude_left' in aConfig_in:
            self.dLongitude_left = float(aConfig_in['dLongitude_left'])

        if 'dLongitude_right' in aConfig_in:
            self.dLongitude_right = float(aConfig_in['dLongitude_right'])

        if 'dLatitude_bot' in aConfig_in:
            self.dLatitude_bot = float(aConfig_in['dLatitude_bot'])

        if 'dLatitude_top' in aConfig_in:
            self.dLatitude_top = float(aConfig_in['dLatitude_top'])

        if 'sFilename_model_configuration' in aConfig_in:
            self.sFilename_model_configuration = aConfig_in['sFilename_model_configuration']

        if 'sFilename_mesh_boundary' in aConfig_in:
            self.sFilename_mesh_boundary = aConfig_in['sFilename_mesh_boundary']

        if 'sFilename_coastal_boundary' in aConfig_in:
            self.sFilename_coastal_boundary = aConfig_in['sFilename_coastal_boundary']
            self.iFlag_land_ocean_mask = 1
        else:
            self.iFlag_land_ocean_mask = 0

        if 'sFilename_spatial_reference' in aConfig_in:
            self.sFilename_spatial_reference = aConfig_in['sFilename_spatial_reference']

        if 'sFilename_dem' in aConfig_in:
            self.sFilename_dem = aConfig_in['sFilename_dem']

        if 'sFilename_mpas_mesh_netcdf' in aConfig_in:
            self.sFilename_mpas_mesh_netcdf = aConfig_in['sFilename_mpas_mesh_netcdf']

        if 'sFilename_jigsaw_mesh_netcdf' in aConfig_in:
            self.sFilename_jigsaw_mesh_netcdf = aConfig_in['sFilename_jigsaw_mesh_netcdf']

        if 'sWorkspace_bin' in aConfig_in:
            self.sWorkspace_bin = aConfig_in['sWorkspace_bin']

        if 'sWorkspace_input' in aConfig_in:
            self.sWorkspace_input = aConfig_in['sWorkspace_input']

        if sWorkspace_output_in is not None:
            self.sWorkspace_output = sWorkspace_output_in
        else:
            if 'sWorkspace_output' in aConfig_in:
                self.sWorkspace_output = aConfig_in['sWorkspace_output']

        if 'sJob' in aConfig_in:
            self.sJob = aConfig_in['sJob']
        if 'sRegion' in aConfig_in:
            self.sRegion = aConfig_in['sRegion']

        if sModel_in is not None:
            self.sModel = sModel_in
        else:
            if 'sModel' in aConfig_in:
                self.sModel = aConfig_in['sModel']

        sDate = aConfig_in['sDate']
        if sDate is not None:
            self.sDate = sDate
        else:
            self.sDate = sDate_default

        sCase = self.sModel + self.sDate + sCase_index
        self.sCase = sCase

        if 'sMesh_type' in aConfig_in:
            self.sMesh_type = aConfig_in['sMesh_type']
        else:
            self.sMesh_type = 'hexagon'

        sMesh_type = self.sMesh_type
        if sMesh_type == 'hexagon':  # hexagon
            self.iMesh_type = 1
        else:
            if sMesh_type == 'square':  # square
                self.iMesh_type = 2
            else:
                if sMesh_type == 'latlon':  # latlon
                    self.iMesh_type = 3
                else:
                    if sMesh_type == 'mpas':  # mpas
                        self.iMesh_type = 4
                    else:
                        if sMesh_type == 'dggrid':
                            self.iMesh_type = 5
                            self.iFlag_run_dggrid = 1
                        else:
                            if sMesh_type == 'tin':
                                self.iMesh_type = 6
                            else:
                                print('Unsupported mesh type?')

        if self.sMesh_type == 'dggrid':
            self.dResolution_meter = dggrid_find_resolution_by_index(
                self.sDggrid_type, self.iResolution_index)

        if self.iFlag_run_dggrid == 1:
            if 'sFilename_dggrid' in aConfig_in:
                self.sFilename_dggrid = aConfig_in['sFilename_dggrid']

        # the model can be run as part of hexwatershed or standalone
        if self.iFlag_standalone == 1:
            # in standalone case, will add case information and update output path
            sPath = str(Path(self.sWorkspace_output) / sCase)
            self.sWorkspace_output = sPath
        else:
            # use specified output path, also do not add output or input tag
            sPath = self.sWorkspace_output

        if iFlag_create_directory_in == 1:
            Path(sPath).mkdir(parents=True, exist_ok=True)

        self.aBasin = list()
        if self.iFlag_flowline == 1:
            if 'sFilename_basins' in aConfig_in:
                self.sFilename_basins = aConfig_in['sFilename_basins']
                if os.path.isfile(self.sFilename_basins):
                    with open(self.sFilename_basins) as json_file:
                        dummy_data = json.load(json_file)
                        for i in range(self.nOutlet):
                            sBasin = "{:08d}".format(i+1)
                            dummy_basin = dummy_data[i]
                            dummy_basin['sWorkspace_output_basin'] = str(
                                Path(self.sWorkspace_output) / sBasin)
                            pBasin = pybasin(dummy_basin)
                            self.aBasin.append(pBasin)
                else:
                    print('This basin configuration file does not exist: ',
                          self.sFilename_basins)
                    print('Please update this parameter before setting up the model!')
            else:
                pass

        if self.iFlag_run_jigsaw == 1:
            if 'sFilename_jigsaw_configuration' in aConfig_in:
                self.sFilename_jigsaw_configuration = aConfig_in['sFilename_jigsaw_configuration']
                if os.path.isfile(self.sFilename_jigsaw_configuration):
                    with open(self.sFilename_jigsaw_configuration) as json_file:
                        self.aConfig_jigsaw = json.load(json_file)  # dict
            else:
                print('Please provide the jigsaw binary file!')
                return

        # converted geojson files
        self.sFilename_watershed_boundary_geojson = os.path.join(
            str(Path(self.sWorkspace_output)), 'watershed_boundary.geojson')
        self.sFilename_mesh_boundary_geojson = os.path.join(
            str(Path(self.sWorkspace_output)), 'mesh_boundary.geojson')
        self.sFilename_coastal_boundary_geojson = os.path.join(
            str(Path(self.sWorkspace_output)), 'coastal_boundary.geojson')
        # model generated files
        self.sFilename_mesh = os.path.join(
            str(Path(self.sWorkspace_output)), sMesh_type + ".geojson")
        self.sFilename_mesh_info = os.path.join(
            str(Path(self.sWorkspace_output)), sMesh_type + "_mesh_info.json")
        self.sFilename_mesh_kml = os.path.join(
            # for google service
            str(Path(self.sWorkspace_output)), sMesh_type + ".kml")
        self.sFilename_mesh_parquet = os.path.join(
            str(Path(self.sWorkspace_output)), sMesh_type + ".parquet")
        return

    def pyflowline_set_boundary_based_on_dem(self, sFilename_dem):
        if os.path.isfile(sFilename_dem):
            dummy = gdal_read_geotiff_file(
                self.sFilename_dem, iFlag_metadata_only=1)
            dPixelWidth = dummy['pixelWidth']
            pPixelHeight = dummy['pixelHeight']
            dOriginX = dummy['originX']
            dOriginY = dummy['originY']
            nrow = dummy['nrow']
            ncolumn = dummy['ncolumn']
            pProjection_dem = dummy['projection']
            dX_lowerleft = dOriginX
            dY_lowerleft = dOriginY + nrow * pPixelHeight
            dLongitude_left0,  dLatitude_bot0 = reproject_coordinates(
                dX_lowerleft, dY_lowerleft, pProjection_dem, self.pProjection_model)
            dX_upperright = dOriginX + ncolumn * dPixelWidth
            dY_upperright = dOriginY
            dLongitude_right0, dLatitude_top0 = reproject_coordinates(
                dX_upperright, dY_upperright, pProjection_dem, self.pProjection_model)
            dX_lowerright = dOriginX + ncolumn * dPixelWidth
            dY_lowerright = dOriginY + nrow * pPixelHeight
            dLongitude_right1,  dLatitude_bot1 = reproject_coordinates(
                dX_lowerright, dY_lowerright, pProjection_dem, self.pProjection_model)
            dX_upperleft = dOriginX
            dY_upperleft = dOriginY
            dLongitude_left1, dLatitude_top1 = reproject_coordinates(
                dX_upperleft, dY_upperleft, pProjection_dem,  self.pProjection_model)
            self.dLatitude_top = np.max([dLatitude_top0, dLatitude_top1])
            self.dLatitude_bot = np.min([dLatitude_bot0, dLatitude_bot1])
            self.dLongitude_left = np.min([dLongitude_left0, dLongitude_left1])
            self.dLongitude_right = np.max(
                [dLongitude_right0, dLongitude_right1])
            self.dLongitude_mean = 0.5 * \
                (self.dLongitude_left + self.dLongitude_right)
            self.dLatitude_mean = 0.5 * \
                (self.dLatitude_bot + self.dLatitude_top)
            self.dX_lowerleft = dX_lowerleft
            self.dY_lowerleft = dY_lowerleft
            self.dX_upperleft = dX_upperleft
            self.dY_upperleft = dY_upperleft
            self.dX_lowerright = dX_lowerright
            self.dY_lowerright = dY_lowerright
            if self.dResolution_meter is not None:  # convert user provide meter to degree
                self.dResolution_degree = meter_to_degree(
                    self.dResolution_meter, self.dLatitude_mean)
        else:
            print('The DEM file does not exist!')

        return

    def pyflowline_setup(self):
        """
        Set up the flowlinecase
        """
        self.pyflowline_check_parameters()
        system = platform.system()

        if self.iFlag_flowline == 1:
            self.pyflowline_convert_flowline_to_geojson()
            pass
        else:
            # at least we will creat a folder

            pass

        if self.iFlag_watershed_boundary == 1:
            self.pyflowline_convert_watershed_boundary_to_geojson()
            pass

        if self.iFlag_mesh_boundary == 1:
            self.pyflowline_convert_mesh_boundary_to_geojson()
            pass

        if self.iFlag_run_jigsaw == 1:
            # create dggrid output folder
            if self.iFlag_standalone == 1:
                sWorkspace_output = os.path.join(self.sWorkspace_output , 'jigsaw')
            else:
                sWorkspace_output = os.path.join(self.sWorkspace_output, '..', 'jigsaw')
                sWorkspace_output = os.path.abspath(sWorkspace_output)

            # check if the folder exists
            # if os.path.exists(sWorkspace_output):
            #    #delete the folder
            #    shutil.rmtree(sWorkspace_output)

            Path(sWorkspace_output).mkdir(parents=True, exist_ok=True)
            self.sWorkspace_jigsaw = sWorkspace_output
            # then copy the binary file to the folder
            # copy execulate
            if self.iFlag_user_provided_binary == 1:
                sFilename_new = os.path.join(sWorkspace_output , 'jigsaw')
                copy2(self.sFilename_dggrid, sFilename_new)
                os.chmod(sFilename_new, stat.S_IREAD |
                         stat.S_IWRITE | stat.S_IXUSR)
                # make this binary as the default one and overwrite other binaries
                # Add the destination directory ahead of the system PATH
                os.environ["PATH"] = sWorkspace_output + \
                    os.pathsep + os.environ["PATH"]
                # Verify that the binary is in the PATH and is the default
                pass
            else:
                if system == 'Windows':
                    print('Windows system is not supported for jigsaw mesh generation!')
                    exit()
                else:
                    sFilename_executable = 'jigsaw'
                # search for system wide binary in the system path
                iFlag_found_binary = 0
                for folder in os.environ['PATH'].split(os.pathsep):
                    sFilename_jigsaw_bin = os.path.join(
                        folder, sFilename_executable)
                    if os.path.isfile(sFilename_jigsaw_bin):
                        print('Found binary at:', sFilename_jigsaw_bin)
                        iFlag_found_binary = 1
                        break
                else:
                    print('Binary not found in system path.')
                if iFlag_found_binary == 1:
                    sFilename_new = os.path.join(sWorkspace_output , 'jigsaw')
                    copy2(sFilename_jigsaw_bin, sFilename_new)
                    os.chmod(sFilename_new, stat.S_IREAD |
                             stat.S_IWRITE | stat.S_IXUSR)
                else:
                    print('Binary not found.')
                    return
            pass
        if self.iFlag_run_dggrid == 1:
            # create dggrid output folder
            if self.iFlag_standalone == 1:
                sWorkspace_output = os.path.join(self.sWorkspace_output , 'dggrid')
            else:
                sWorkspace_output = os.path.join(self.sWorkspace_output, '..', 'dggrid')
                sWorkspace_output = os.path.abspath(sWorkspace_output)

            if os.path.exists(sWorkspace_output):
                # delete the folder
                shutil.rmtree(sWorkspace_output)

            Path(sWorkspace_output).mkdir(parents=True, exist_ok=True)
            self.sWorkspace_dggrid = sWorkspace_output
            # then copy the binary file to the folder
            # copy execulate
            if self.iFlag_user_provided_binary == 1:
                sFilename_new = os.path.join(sWorkspace_output, 'dggrid')
                copy2(self.sFilename_dggrid, sFilename_new)
                os.chmod(sFilename_new, stat.S_IREAD |
                         stat.S_IWRITE | stat.S_IXUSR)
                pass
            else:
                if system == 'Windows':
                    sFilename_executable = 'dggrid.exe'
                else:
                    sFilename_executable = 'dggrid'
                # search for system wide binary in the system path
                iFlag_found_binary = 0
                for folder in os.environ['PATH'].split(os.pathsep):
                    sFilename_dggrid_bin = os.path.join(
                        folder, sFilename_executable)
                    if os.path.isfile(sFilename_dggrid_bin):
                        print('Found binary at:', sFilename_dggrid_bin)
                        iFlag_found_binary = 1
                        break
                else:
                    print('Binary not found in system path.')
                if iFlag_found_binary == 1:
                    sFilename_new = os.path.join(sWorkspace_output, 'dggrid')
                    copy2(sFilename_dggrid_bin, sFilename_new)
                    os.chmod(sFilename_new, stat.S_IREAD |
                             stat.S_IWRITE | stat.S_IXUSR)
                else:
                    print('Binary not found in system path.')
                    return

        return

    def pyflowline_check_parameters(self):
        # this function check the user input
        pSpatial_reference_model = osr.SpatialReference()
        pSpatial_reference_model.ImportFromEPSG(4326)
        self.pProjection_model = pSpatial_reference_model.ExportToWkt()
        sMesh_type = self.sMesh_type
        if self.iFlag_global == 1:  # a global mesh does not require boundary
            # other requirements maybe needed
            pass
        else:
            if self.iFlag_antarctic == 1 and self.iFlag_arctic == 1:
                print(
                    'Both antarctic and arctic flags are on, please check the parameters!')

            if self.iFlag_mesh_boundary == 1:
                if not os.path.isfile(self.sFilename_mesh_boundary):
                    print(
                        "The mesh boundary file does not exist, you should update this parameter before running the model!")
                else:
                    if not os.path.isfile(self.sFilename_coastal_boundary):
                        print(
                            "The coastal boundary file does not exist, the mesh boundary will be used instead!")
                        self.sFilename_coastal_boundary = self.sFilename_mesh_boundary
                    pass

            if sMesh_type == 'hexagon':  # hexagon #need spatial referece
                # check boundary
                if self.iFlag_mesh_boundary == 1:
                    # still need spatial reference, no
                    if os.path.isfile(self.sFilename_spatial_reference):
                        # check vector or raster
                        if (gdal_check_file_type(self.sFilename_spatial_reference) == 'raster'):
                            self.pProjection_reference = gdal_get_raster_spatial_reference_wkt(
                                self.sFilename_spatial_reference)
                        else:
                            if (gdal_check_file_type(self.sFilename_spatial_reference) == 'vector'):
                                self.pProjection_reference = gdal_get_vector_spatial_reference_wkt(
                                    self.sFilename_spatial_reference)
                    else:
                        if os.path.isfile(self.sFilename_dem):
                            self.pProjection_reference = gdal_get_raster_spatial_reference_wkt(
                                self.sFilename_dem)
                        else:
                            # use utm in this case, will depeends on longitude
                            pass
                else:
                    # check DEM or spatial reference file
                    if os.path.isfile(self.sFilename_dem):
                        # get extent
                        self.pyflowline_set_boundary_based_on_dem(
                            self.sFilename_dem)

                    else:
                        # use user provided extent
                        # get the spatial reference file
                        if os.path.isfile(self.sFilename_spatial_reference):
                            # check vector or raster
                            if (gdal_check_file_type(self.sFilename_spatial_reference) == 'raster'):
                                self.pProjection_reference = gdal_get_raster_spatial_reference_wkt(
                                    self.sFilename_spatial_reference)
                            else:
                                if (gdal_check_file_type(self.sFilename_spatial_reference) == 'vector'):
                                    self.pProjection_reference = gdal_get_vector_spatial_reference_wkt(
                                        self.sFilename_spatial_reference)
                        else:
                            print(
                                'No spatial reference information can be defined because there is no DEM or spatial reference file!')

            else:
                if sMesh_type == 'square':  # square #need spatial referece
                    if self.iFlag_mesh_boundary == 1:
                        # still need spatial reference, no
                        if os.path.isfile(self.sFilename_spatial_reference):
                            # check vector or raster
                            if (gdal_check_file_type(self.sFilename_spatial_reference) == 'raster'):
                                self.pProjection_reference = gdal_get_raster_spatial_reference_wkt(
                                    self.sFilename_spatial_reference)
                            else:
                                if (gdal_check_file_type(self.sFilename_spatial_reference) == 'vector'):
                                    self.pProjection_reference = gdal_get_vector_spatial_reference_wkt(
                                        self.sFilename_spatial_reference)
                        else:
                            if os.path.isfile(self.sFilename_dem):
                                self.pProjection_reference = gdal_get_raster_spatial_reference_wkt(
                                    self.sFilename_dem)
                            else:
                                # use utm in this case, will depeends on longitude
                                pass
                    else:
                        # check DEM or spatial reference file
                        if os.path.isfile(self.sFilename_dem):
                            # get extent
                            self.pyflowline_set_boundary_based_on_dem(
                                self.sFilename_dem)

                        else:
                            # use user provided extent
                            # get the spatial reference file
                            if os.path.isfile(self.sFilename_spatial_reference):
                                # check vector or raster
                                if (gdal_check_file_type(self.sFilename_spatial_reference) == 'raster'):
                                    self.pProjection_reference = gdal_get_raster_spatial_reference_wkt(
                                        self.sFilename_spatial_reference)
                                else:
                                    if (gdal_check_file_type(self.sFilename_spatial_reference) == 'vector'):
                                        self.pProjection_reference = gdal_get_vector_spatial_reference_wkt(
                                            self.sFilename_spatial_reference)
                            else:
                                print(
                                    'No spatial reference information can be defined because there is no DEM or spatial reference file!')

                else:
                    if sMesh_type == 'latlon':  # latlon #do not need spatial referece
                        if self.iFlag_mesh_boundary == 1:
                            pass
                        else:
                            if os.path.isfile(self.sFilename_dem):
                                self.pyflowline_set_boundary_based_on_dem(
                                    self.sFilename_dem)
                            else:
                                # use user provided extent
                                pass
                    else:
                        if sMesh_type == 'mpas':  # mpas #do not need spatial referece
                            if self.iFlag_run_jigsaw == 1:
                                # if the mesh is generated by jigsaw, then we do not need to check the mesh file
                                # but do we need the boundary file? yes
                                pass
                            else:
                                if not os.path.isfile(self.sFilename_mpas_mesh_netcdf):
                                    print("The MPAS mesh file does not exist!")

                            if self.iFlag_mesh_boundary == 1:
                                pass
                            else:
                                if os.path.isfile(self.sFilename_dem):
                                    self.pyflowline_set_boundary_based_on_dem(
                                        self.sFilename_dem)

                        else:
                            if sMesh_type == 'dggrid':  # do not need spatial referece
                                self.dResolution_meter = dggrid_find_resolution_by_index(
                                    self.sDggrid_type, self.iResolution_index)
                                if self.iFlag_user_provided_binary == 1:
                                    if not os.path.isfile(self.sFilename_dggrid):
                                        print(
                                            "The dggrid binary file does not exist, you need to update this parameter before running the model!")
                                    pass
                                if self.iFlag_mesh_boundary == 1:
                                    pass
                                else:
                                    if os.path.isfile(self.sFilename_dem):
                                        self.pyflowline_set_boundary_based_on_dem(
                                            self.sFilename_dem)
                            else:
                                if sMesh_type == 'tin':  # need spatial referece
                                    if self.iFlag_mesh_boundary == 1:
                                        pass
                                    else:
                                        if os.path.isfile(self.sFilename_dem):
                                            self.pyflowline_set_boundary_based_on_dem(
                                                self.sFilename_dem)
                                        pass
                                else:
                                    print('Unsupported mesh type?')

        return

    def pyflowline_change_model_parameter(self, sVariable_in, dValue, iFlag_basin_in=None):
        if iFlag_basin_in is None:
            if hasattr(self, sVariable_in):
                # get default data type
                sType_default = type(getattr(self, sVariable_in))
                # get the data type of the input value
                sType_input = type(dValue)
                if sType_default == sType_input:
                    setattr(self, sVariable_in, dValue)
                    pass
                else:
                    print('Incorrect data type for the input value: ' + sVariable_in)
                return True
            else:
                print(
                    "This model parameter is unknown, please check the full parameter list in the documentation: " + sVariable_in)
                return False

        else:
            # this mean the variable is in the basin object
            for pBasin in self.aBasin:
                if hasattr(pBasin, sVariable_in):
                    # get default data type
                    sType_default = type(getattr(pBasin, sVariable_in))
                    sType_input = type(dValue)
                    if sType_default == sType_input:
                        setattr(pBasin, sVariable_in, dValue)
                    else:
                        print(
                            'Incorrect data type for the input value: ' + sVariable_in)
                        return False
                else:
                    print(
                        "This model parameter is unknown, please check the full parameter list in the documentation: " + sVariable_in)
                    return False

    def pyflowline_convert_flowline_to_geojson(self):
        if self.iFlag_flowline == 1:
            for pBasin in self.aBasin:
                pBasin.basin_convert_flowline_to_geojson()
                pass
        return

    def pyflowline_convert_watershed_boundary_to_geojson(self):
        # convert watershed boundary to geojson
        # this boundary is currently only used for mesh generation purpose
        # in the future, there may be a list of boundary file for mesh generation
        # we can either provide a merged boundary file or a list of boundary files

        iFlag_future = 0
        if iFlag_future == 1:
            pass
        else:
            # single watershed boundary
            sFilename_raw = self.sFilename_watershed_boundary
            sFilename_out = self.sFilename_watershed_boundary_geojson
            convert_boundary_to_geojson(sFilename_raw, sFilename_out)
            pass

        # if self.iFlag_watershed_boundary == 1:
        #    for pBasin in self.aBasin:
        #        pBasin.basin_convert_boundary_to_geojson()
        #        pass

        return

    def pyflowline_convert_mesh_boundary_to_geojson(self):
        sFilename_raw = self.sFilename_mesh_boundary
        sFilename_out = self.sFilename_mesh_boundary_geojson
        # check whether the file exists
        if os.path.isfile(sFilename_raw):
            convert_boundary_to_geojson(sFilename_raw, sFilename_out)
        else:
            print('The mesh boundary file does not exist!')
            return

        # get mean location
        pBoundary_wkt, aExtent = gdal_read_geojson_boundary(
            self.sFilename_mesh_boundary_geojson)
        self.dLongitude_mean = 0.5 * (aExtent[0] + aExtent[2])
        self.dLatitude_mean = 0.5 * (aExtent[1] + aExtent[3])

        #now do the coastal boundary
        if self.iFlag_land_ocean_mask == 1:
            sFilename_raw = self.sFilename_coastal_boundary
            sFilename_out = self.sFilename_coastal_boundary_geojson
            # check whether the file exists
            if os.path.isfile(sFilename_raw):
                convert_boundary_to_geojson(sFilename_raw, sFilename_out)
            else:
                print('The coastal boundary file does not exist!')
        return

    def pyflowline_flowline_simplification(self):
        aFlowline_out = list()  # store all the flowline
        # export the outlet into a single file
        aOutlet = list()
        if self.iFlag_simplification == 1:
            for i in range(self.nOutlet):
                pBasin = self.aBasin[i]
                aFlowline_basin = pBasin.basin_flowline_simplification()
                aFlowline_out = aFlowline_out + aFlowline_basin
                aOutlet.append(pBasin.pVertex_outlet)

            self.aFlowline_simplified = aFlowline_out

        return aFlowline_out

    def pyflowline_mesh_generation(self, iFlag_antarctic_in=None, iFlag_arctic_in=None):
        """
        The mesh generation operation

        Returns:
            list [pycell]: A list of cell object
        """
        print('Start mesh generation.')

        if iFlag_antarctic_in is None:
            iFlag_antarctic = self.iFlag_antarctic
        else:
            iFlag_antarctic = iFlag_antarctic_in

        if iFlag_arctic_in is None:
            iFlag_arctic = self.iFlag_arctic
        else:
            iFlag_arctic = iFlag_arctic_in

        if self.iFlag_create_mesh == 1:
            iFlag_global = self.iFlag_global
            iMesh_type = self.iMesh_type
            iFlag_save_mesh = self.iFlag_save_mesh
            iFlag_rotation = self.iFlag_rotation
            iFlag_mesh_boundary = self.iFlag_mesh_boundary
            dLatitude_mean = self.dLatitude_mean
            dResolution_degree = self.dResolution_degree
            dResolution_meter = self.dResolution_meter
            sFilename_mesh = self.sFilename_mesh
            if iMesh_type == 1:  # hexagon
                # hexagon edge
                dResolution_meter = degree_to_meter(
                    dResolution_degree, dLatitude_mean)
                dArea = np.power(dResolution_meter, 2.0)
                dLength_edge = np.sqrt(2.0 * dArea / (3.0 * np.sqrt(3.0)))
                if iFlag_rotation == 0:
                    dX_spacing = dLength_edge * np.sqrt(3.0)
                    dY_spacing = dLength_edge * 1.5
                    ncolumn = int(
                        (self.dX_lowerright - self.dX_lowerleft) / dX_spacing)
                    nrow = int(
                        (self.dY_upperleft - self.dY_lowerleft) / dY_spacing)
                else:
                    dX_spacing = dLength_edge * 1.5
                    dY_spacing = dLength_edge * np.sqrt(3.0)
                    ncolumn = int(
                        (self.dX_lowerright - self.dX_lowerleft) / dX_spacing)+1
                    nrow = int(
                        (self.dY_upperleft - self.dY_lowerleft) / dY_spacing)

                if iFlag_mesh_boundary == 1:
                    # create a polygon based on real boundary
                    pBoundary_wkt, aExtent = gdal_read_geojson_boundary(
                        self.sFilename_mesh_boundary_geojson)
                    aHexagon = create_hexagon_mesh(iFlag_rotation, self.dX_lowerleft, self.dY_lowerleft, dResolution_meter, ncolumn, nrow,
                                                   sFilename_mesh, self.pProjection_reference, pBoundary_wkt)
                    pass
                else:
                    pRing = ogr.Geometry(ogr.wkbLinearRing)
                    pRing.AddPoint(dLongitude_left, dLatitude_top)
                    pRing.AddPoint(dLongitude_right, dLatitude_top)
                    pRing.AddPoint(dLongitude_right, dLatitude_bot)
                    pRing.AddPoint(dLongitude_left, dLatitude_bot)
                    pRing.AddPoint(dLongitude_left, dLatitude_top)
                    pBoundary = ogr.Geometry(ogr.wkbPolygon)
                    pBoundary.AddGeometry(pRing)
                    pBoundary_wkt = pBoundary.ExportToWkt()  # wkt format
                    aHexagon = create_hexagon_mesh(iFlag_rotation, self.dX_lowerleft, self.dY_lowerleft, dResolution_meter,
                                                   ncolumn, nrow,
                                                   sFilename_mesh, self.pProjection_reference, pBoundary_wkt)
                    pass
                self.aCell = aHexagon
            else:
                if iMesh_type == 2:  # sqaure
                    ncolumn = int(
                        (self.dX_lowerright - self.dX_lowerleft) / dResolution_meter)
                    nrow = int(
                        (self.dY_upperleft - self.dY_lowerleft) / dResolution_meter)
                    if iFlag_mesh_boundary == 1:
                        pBoundary_wkt, aExtent = gdal_read_geojson_boundary(
                            self.sFilename_mesh_boundary_geojson)
                        aSquare = create_square_mesh(self.dX_lowerleft, self.dY_lowerleft, dResolution_meter,
                                                     ncolumn, nrow,
                                                     sFilename_mesh, self.pProjection_reference, pBoundary_wkt)
                        pass
                    else:
                        pRing = ogr.Geometry(ogr.wkbLinearRing)
                        pRing.AddPoint(dLongitude_left, dLatitude_top)
                        pRing.AddPoint(dLongitude_right, dLatitude_top)
                        pRing.AddPoint(dLongitude_right, dLatitude_bot)
                        pRing.AddPoint(dLongitude_left, dLatitude_bot)
                        pRing.AddPoint(dLongitude_left, dLatitude_top)
                        pBoundary = ogr.Geometry(ogr.wkbPolygon)
                        pBoundary.AddGeometry(pRing)
                        pBoundary_wkt = pBoundary.ExportToWkt()
                        aSquare = create_square_mesh(self.dX_lowerleft, self.dY_lowerleft, dResolution_meter,
                                                     ncolumn, nrow,
                                                     sFilename_mesh, self.pProjection_reference, pBoundary_wkt)
                        pass

                    self.aCell = aSquare
                else:
                    if iMesh_type == 3:  # latlon
                        dResolution_meter = degree_to_meter(
                            dResolution_degree, self.dLatitude_mean)
                        dArea = np.power(dResolution_meter, 2.0)
                        ncolumn = int(
                            (dLongitude_right - dLongitude_left) / dResolution_degree)
                        nrow = int((dLatitude_top - dLatitude_bot) /
                                   dResolution_degree)
                        if iFlag_mesh_boundary == 1:
                            # create a polygon based on real boundary
                            # already produced
                            pBoundary_wkt, aExtent = gdal_read_geojson_boundary(
                                self.sFilename_mesh_boundary_geojson)
                            aLatlon = create_latlon_mesh(dLongitude_left, dLatitude_bot, dResolution_degree,
                                                         ncolumn, nrow,
                                                         sFilename_mesh, pBoundary_wkt)
                            pass
                        else:
                            pRing = ogr.Geometry(ogr.wkbLinearRing)
                            pRing.AddPoint(dLongitude_left, dLatitude_top)
                            pRing.AddPoint(dLongitude_right, dLatitude_top)
                            pRing.AddPoint(dLongitude_right, dLatitude_bot)
                            pRing.AddPoint(dLongitude_left, dLatitude_bot)
                            pRing.AddPoint(dLongitude_left, dLatitude_top)
                            pBoundary = ogr.Geometry(ogr.wkbPolygon)
                            pBoundary.AddGeometry(pRing)
                            pBoundary_wkt = pBoundary.ExportToWkt()
                            aLatlon = create_latlon_mesh(dLongitude_left, dLatitude_bot, dResolution_degree,
                                                         ncolumn, nrow,
                                                         sFilename_mesh, pBoundary_wkt)

                        self.aCell = aLatlon
                    else:
                        if iMesh_type == 4:  # mpas
                            sFilename_mpas_mesh_netcdf = self.sFilename_mpas_mesh_netcdf
                            sFilename_jigsaw_mesh_netcdf = self.sFilename_jigsaw_mesh_netcdf
                            iFlag_use_mesh_dem = self.iFlag_use_mesh_dem
                            dLatitude_top = self.dLatitude_top
                            dLatitude_bot = self.dLatitude_bot
                            dLongitude_left = self.dLongitude_left
                            dLongitude_right = self.dLongitude_right
                            if iFlag_antarctic == 1 or iFlag_arctic == 1:
                                aMpas = create_mpas_mesh(sFilename_mesh,
                                                         iFlag_global_in=iFlag_global,
                                                         iFlag_use_mesh_dem_in=iFlag_use_mesh_dem,
                                                         iFlag_save_mesh_in=iFlag_save_mesh,
                                                         sFilename_mpas_mesh_netcdf_in=sFilename_mpas_mesh_netcdf,

                                                         iFlag_run_jigsaw_in=self.iFlag_run_jigsaw,
                                                         iFlag_antarctic_in=iFlag_antarctic,
                                                         iFlag_arctic_in=iFlag_arctic,
                                                         aConfig_jigsaw_in=self.aConfig_jigsaw,
                                                         sWorkspace_jigsaw_in=self.sWorkspace_jigsaw)
                                pass
                            else:
                                if iFlag_mesh_boundary == 1:
                                    # create a polygon based on
                                    # read boundary
                                    pBoundary_wkt, aExtent = gdal_read_geojson_boundary(
                                        self.sFilename_mesh_boundary_geojson)
                                    print(sFilename_mpas_mesh_netcdf)
                                    # if user wants to use the coastal boundary, then we can replace the default land ocean mask
                                    if self.iFlag_land_ocean_mask == 1:
                                        sFilename_land_ocean_mask_in = self.sFilename_coastal_boundary_geojson
                                    else:
                                        sFilename_land_ocean_mask_in = None
                                        pass
                                    aMpas = create_mpas_mesh(sFilename_mesh,
                                                             iFlag_global_in=iFlag_global,
                                                             iFlag_use_mesh_dem_in=iFlag_use_mesh_dem,
                                                             iFlag_save_mesh_in=iFlag_save_mesh,
                                                             iFlag_run_jigsaw_in=self.iFlag_run_jigsaw,
                                                             iFlag_antarctic_in=iFlag_antarctic_in,
                                                             iFlag_arctic_in=iFlag_arctic_in,
                                                             sWorkspace_jigsaw_in=self.sWorkspace_jigsaw,
                                                             sFilename_mpas_mesh_netcdf_in=sFilename_mpas_mesh_netcdf,
                                                             sFilename_jigsaw_mesh_netcdf_in=sFilename_jigsaw_mesh_netcdf,
                                                             sFilename_land_ocean_mask_in=sFilename_land_ocean_mask_in,
                                                             aConfig_jigsaw_in=self.aConfig_jigsaw,
                                                             pBoundary_in=pBoundary_wkt)

                                else:
                                    pRing = ogr.Geometry(ogr.wkbLinearRing)
                                    pRing.AddPoint(
                                        dLongitude_left, dLatitude_top)
                                    pRing.AddPoint(
                                        dLongitude_right, dLatitude_top)
                                    pRing.AddPoint(
                                        dLongitude_right, dLatitude_bot)
                                    pRing.AddPoint(
                                        dLongitude_left, dLatitude_bot)
                                    pRing.AddPoint(
                                        dLongitude_left, dLatitude_top)
                                    pBoundary = ogr.Geometry(ogr.wkbPolygon)
                                    pBoundary.AddGeometry(pRing)
                                    pBoundary_wkt = pBoundary.ExportToWkt()
                                    # new method using polygon object
                                    aMpas = create_mpas_mesh(sFilename_mesh,
                                                             iFlag_global_in=iFlag_global,
                                                             iFlag_use_mesh_dem_in=iFlag_use_mesh_dem,
                                                             iFlag_save_mesh_in=iFlag_save_mesh,
                                                             iFlag_run_jigsaw_in=self.iFlag_run_jigsaw,
                                                             iFlag_antarctic_in=iFlag_antarctic_in,
                                                             iFlag_arctic_in=iFlag_arctic_in,
                                                             sWorkspace_jigsaw_in=self.sWorkspace_jigsaw,
                                                             sFilename_mpas_mesh_netcdf_in=sFilename_mpas_mesh_netcdf,
                                                             sFilename_jigsaw_mesh_netcdf_in=sFilename_jigsaw_mesh_netcdf,
                                                             aConfig_jigsaw_in=self.aConfig_jigsaw,
                                                             pBoundary_in=pBoundary_wkt)
                                    pass

                            self.aCell = aMpas
                        else:
                            if iMesh_type == 5:  # dggrid
                                dLatitude_top = self.dLatitude_top
                                dLatitude_bot = self.dLatitude_bot
                                dLongitude_left = self.dLongitude_left
                                dLongitude_right = self.dLongitude_right
                                if self.iFlag_standalone == 1:
                                    sWorkspace_output = os.path.join(self.sWorkspace_output, 'dggrid')
                                    pass
                                else:
                                    sWorkspace_output = os.path.join(self.sWorkspace_output, '..', 'dggrid')
                                    sWorkspace_output = os.path.abspath(
                                        sWorkspace_output)
                                    pass

                                if not os.path.exists(sWorkspace_output):
                                    os.makedirs(sWorkspace_output)

                                if iFlag_mesh_boundary == 1:
                                    # create a polygon based on
                                    aDggrid = create_dggrid_mesh(iFlag_global,
                                                                 iFlag_save_mesh,
                                                                 sFilename_mesh,
                                                                 sWorkspace_output,
                                                                 iResolution_index_in=self.iResolution_index,
                                                                 iFlag_antarctic_in=iFlag_antarctic_in,
                                                                 iFlag_arctic_in=iFlag_arctic_in,
                                                                 sDggrid_type_in=self.sDggrid_type,
                                                                 sFilename_boundary_in=self.sFilename_mesh_boundary_geojson)

                                    pass
                                else:
                                    aDggrid = create_dggrid_mesh(iFlag_global,
                                                                 iFlag_save_mesh,
                                                                 sFilename_mesh,
                                                                 sWorkspace_output,
                                                                 iResolution_index_in=self.iResolution_index,
                                                                 sDggrid_type_in=self.sDggrid_type)

                                    pass

                                self.aCell = aDggrid
                            else:
                                if iMesh_type == 6:  # structured tin this one need to be updated because central location issue?
                                    # tin edge
                                    dArea = np.power(dResolution_meter, 2.0)
                                    dLength_edge = np.sqrt(
                                        4.0 * dArea / np.sqrt(3.0))
                                    dX_shift = 0.5 * dLength_edge
                                    dY_shift = 0.5 * \
                                        dLength_edge * np.sqrt(3.0)
                                    dX_spacing = dX_shift * 2
                                    dY_spacing = dY_shift
                                    ncolumn = int(
                                        (self.dX_lowerright - self.dX_lowerleft) / dX_shift)
                                    nrow = int(
                                        (self.dY_upperleft - self.dY_lowerleft) / dY_spacing)
                                    aTriangular = create_triangular_mesh(self.dX_lowerleft, self.dY_lowerleft, dResolution_meter, ncolumn, nrow,
                                                                         sFilename_mesh,
                                                                         self.pProjection_reference)
                                    self.aCell = aTriangular
                                else:
                                    if iMesh_type == 7:  # tin this one need to be updated to use the jigsaw mesh method
                                        sFilename_jigsaw_mesh_netcdf = self.sFilename_jigsaw_mesh_netcdf
                                        iFlag_use_mesh_dem = self.iFlag_use_mesh_dem
                                        dLatitude_top = self.dLatitude_top
                                        dLatitude_bot = self.dLatitude_bot
                                        dLongitude_left = self.dLongitude_left
                                        dLongitude_right = self.dLongitude_right
                                        if iFlag_antarctic == 1 or iFlag_arctic == 1:
                                            aMpas = create_tin_mesh(iFlag_global, iFlag_use_mesh_dem, iFlag_save_mesh,
                                                                    sFilename_jigsaw_mesh_netcdf, sFilename_mesh,
                                                                    iFlag_run_jigsaw_in=self.iFlag_run_jigsaw,
                                                                    iFlag_antarctic_in=iFlag_antarctic,
                                                                    iFlag_arctic_in=iFlag_arctic,
                                                                    sWorkspace_jigsaw_in=self.sWorkspace_jigsaw)
                                            pass
                                        else:
                                            if iFlag_mesh_boundary == 1:
                                                # create a polygon based on
                                                # read boundary
                                                pBoundary_wkt, aExtent = gdal_read_geojson_boundary(
                                                    self.sFilename_mesh_boundary_geojson)
                                                print(
                                                    sFilename_jigsaw_mesh_netcdf)
                                                aMpas = create_tin_mesh(iFlag_global, iFlag_use_mesh_dem, iFlag_save_mesh,
                                                                        sFilename_jigsaw_mesh_netcdf,  sFilename_mesh,
                                                                        iFlag_run_jigsaw_in=self.iFlag_run_jigsaw,
                                                                        iFlag_antarctic_in=iFlag_antarctic_in,
                                                                        iFlag_arctic_in=iFlag_arctic_in,
                                                                        sWorkspace_jigsaw_in=self.sWorkspace_jigsaw,
                                                                        pBoundary_in=pBoundary_wkt)
                                                pass
                                            else:
                                                pRing = ogr.Geometry(
                                                    ogr.wkbLinearRing)
                                                pRing.AddPoint(
                                                    dLongitude_left, dLatitude_top)
                                                pRing.AddPoint(
                                                    dLongitude_right, dLatitude_top)
                                                pRing.AddPoint(
                                                    dLongitude_right, dLatitude_bot)
                                                pRing.AddPoint(
                                                    dLongitude_left, dLatitude_bot)
                                                pRing.AddPoint(
                                                    dLongitude_left, dLatitude_top)
                                                pBoundary = ogr.Geometry(
                                                    ogr.wkbPolygon)
                                                pBoundary.AddGeometry(pRing)
                                                pBoundary_wkt = pBoundary.ExportToWkt()
                                                # new method using polygon object
                                                aMpas = create_tin_mesh(iFlag_global, iFlag_use_mesh_dem, iFlag_save_mesh,
                                                                        sFilename_jigsaw_mesh_netcdf, sFilename_mesh,
                                                                        iFlag_run_jigsaw_in=self.iFlag_run_jigsaw,
                                                                        iFlag_antarctic_in=iFlag_antarctic_in,
                                                                        iFlag_arctic_in=iFlag_arctic_in,
                                                                        sWorkspace_jigsaw_in=self.sWorkspace_jigsaw,
                                                                        pBoundary_in=pBoundary_wkt)
                                                pass

                                        self.aCell = aMpas
                                    else:
                                        if iMesh_type == 8:  # cubic sphere mesh
                                            dLatitude_top = self.dLatitude_top
                                            dLatitude_bot = self.dLatitude_bot
                                            dLongitude_left = self.dLongitude_left
                                            dLongitude_right = self.dLongitude_right
                                            if self.iFlag_standalone == 1:
                                                sWorkspace_output = self.sWorkspace_output
                                                pass
                                            else:
                                                sWorkspace_output = os.path.join(self.sWorkspace_output, '..', 'tempestremap')
                                                sWorkspace_output = os.path.abspath(
                                                    sWorkspace_output)
                                                pass

                                            if not os.path.exists(sWorkspace_output):
                                                os.makedirs(sWorkspace_output)

                                            if iFlag_mesh_boundary == 1:
                                                # create a polygon based on
                                                aCubicsphere = create_cubicsphere_mesh(iFlag_global,
                                                                                       iFlag_save_mesh,
                                                                                       dResolution_meter,
                                                                                       sFilename_mesh,
                                                                                       sWorkspace_output,
                                                                                       iFlag_antarctic_in=iFlag_antarctic_in,
                                                                                       iFlag_arctic_in=iFlag_arctic_in,
                                                                                       sFilename_boundary_in=self.sFilename_mesh_boundary_geojson)

                                                pass
                                            else:
                                                aCubicsphere = create_cubicsphere_mesh(iFlag_global,
                                                                                       iFlag_save_mesh,
                                                                                       dResolution_meter,
                                                                                       sFilename_mesh,
                                                                                       sWorkspace_output
                                                                                       )

                                                pass
                                            pass
                                        else:

                                            print('Unsupported mesh type?')
                                return

            # no matter what type of mash, we will convert it to geoparquet for easy visualization
            convert_geojson_to_geoparquet(
                sFilename_mesh, self.sFilename_mesh_parquet)
        else:
            # this means the mesh is already generated, we need a function that can read an existing pyflowline compatible mesh
            pass

        # build rtree
        if iFlag_use_rtree == 1:
            interleaved = True
            self.pRTree_mesh = RTree(
                interleaved=interleaved, max_cap=5, min_cap=2)
            for lCellIndex in range(len(self.aCell)):
                pBound = self.aCell[lCellIndex].pBound
                self.pRTree_mesh.insert(lCellIndex, pBound)  #

        print('Finish mesh generation.')
        return self.aCell

    def pyflowline_reconstruct_topological_relationship(self):
        """
        The topological relationship reconstruction operation

        Args:
            aCell_raw (list [pycell]): A list of intersected cell objects

        Returns:
            tuple [list [pycell], list [pyflowline], list [long]]: A list of cells, flowlines, and outlet cell IDs.
        """
        print('Start topology reconstruction.')
        iFlag_intersect = self.iFlag_intersect
        if iFlag_intersect == 1:
            iMesh_type = self.iMesh_type
            sFilename_mesh = self.sFilename_mesh
            aFlowline_conceptual = list()  # store all the flowline
            aCellID_outlet = list()
            aBasin = list()
            aCell_intersect = list()
            ncell = len(self.aCell)
            # there is a one-to-one match between cell id and cell center because each cell has only one center

            for pBasin in self.aBasin:
                aCell_intersect_basin = pBasin.basin_reconstruct_topological_relationship(
                    iMesh_type, sFilename_mesh)
                aFlowline_conceptual = aFlowline_conceptual + pBasin.aFlowline_basin_conceptual
                aBasin.append(pBasin)
                aCellID_outlet.append(pBasin.lCellID_outlet)
                aCell_intersect = aCell_intersect + aCell_intersect_basin
                if iFlag_use_rtree == 1:
                    # use rtree to update topology and length
                    for pFlowline in pBasin.aFlowline_basin_conceptual:
                        iStream_segment = pFlowline.iStream_segment
                        iStream_order = pFlowline.iStream_order
                        for pEdge in pFlowline.aEdge:
                            pVertex_start = pEdge.pVertex_start
                            pVertex_end = pEdge.pVertex_end
                            pBound = pEdge.pBound
                            aIntersect = list(self.pRTree_mesh.search(pBound))
                            for k in aIntersect:
                                pCell = self.aCell[k]
                                if pVertex_start.calculate_distance(pCell.pVertex_center) < 1.0E-6:
                                    self.aCell[k].iStream_segment_burned = iStream_segment
                                    self.aCell[k].iStream_order_burned = iStream_order
                                    lCellIndex_upstream = k

                                if pVertex_end.calculate_distance(pCell.pVertex_center) < 1.0E-6:
                                    self.aCell[k].iStream_segment_burned = iStream_segment
                                    self.aCell[k].iStream_order_burned = iStream_order
                                    lCellIndex_downstream = k

                            self.aCell[lCellIndex_upstream].lCellID_downstream_burned = self.aCell[lCellIndex_downstream].lCellID

                    for pCell2 in aCell_intersect_basin:
                        pBound = pCell2.pBound
                        aIntersect = list(self.pRTree_mesh.search(pBound))
                        for k in aIntersect:
                            pCell = self.aCell[k]
                            if pCell2.lCellID == pCell.lCellID:
                                pCell.dLength_flowline = pCell2.dLength_flowline
                        pass
                else:
                    # set topology here
                    for pFlowline in pBasin.aFlowline_basin_conceptual:
                        iStream_segment = pFlowline.iStream_segment
                        iStream_order = pFlowline.iStream_order
                        for pEdge in pFlowline.aEdge:
                            try:
                                pVertex_start = pEdge.pVertex_start
                                pVertex_end = pEdge.pVertex_end
                                iFlag_found_start = 0
                                for k in range(ncell):
                                    pVertex_center = self.aCell[k].pVertex_center
                                    if pVertex_center == pVertex_start:
                                        iFlag_found_start = 1
                                        self.aCell[k].iStream_segment_burned = iStream_segment
                                        self.aCell[k].iStream_order_burned = iStream_order
                                        iFlag_found_end = 0
                                        for l in range(ncell):
                                            pVertex_center2 = self.aCell[l].pVertex_center
                                            lCellID = self.aCell[l].lCellID
                                            if pVertex_center2 == pVertex_end:
                                                iFlag_found_end = 1
                                                self.aCell[k].lCellID_downstream_burned = lCellID
                                                if lCellID == pBasin.lCellID_outlet:
                                                    self.aCell[l].iStream_segment_burned = iStream_segment
                                                    self.aCell[l].iStream_order_burned = iStream_order
                                                break
                                        if iFlag_found_end == 0:
                                            print('End not found')
                                            pass
                                if iFlag_found_start == 0:
                                    print('Start not found')
                                    pass
                            except:
                                print("error in step")
                                print(pFlowline.tojson())
                                print(pEdge.tojson())
                                print(pVertex_start.tojson())
                                print(pVertex_end.tojson())
                                pass

                    # update length using rtree
                    for pCell in self.aCell:
                        for pCell2 in aCell_intersect_basin:
                            if pCell2.lCellID == pCell.lCellID:
                                pCell.dLength_flowline = pCell2.dLength_flowline

            self.aFlowline_conceptual = aFlowline_conceptual
            self.aCellID_outlet = aCellID_outlet
            print('Finish topology reconstruction.')
            return self.aCell, aFlowline_conceptual, aCellID_outlet
        else:

            return None

    def pyflowline_merge_cell_info(self, aCell_raw):
        """
        Merge cell information after reconstruction

        Args:
            aCell_raw (list [pycell]): The original cell information that contains neighbor definition
            This information is defined in the mesh generation function, so mesh generation must be run.

        Returns:
            list [pycell]: The updated list of cell objects.
        """

        for pCell in self.aCell:
            for pCell2 in aCell_raw:
                if pCell.lCellID == pCell2.lCellID:
                    pCell.aNeighbor = pCell2.aNeighbor
                    pCell.nNeighbor = pCell2.nNeighbor
                    pCell.aNeighbor_land = pCell2.aNeighbor_land
                    pCell.nNeighbor_land = pCell2.nNeighbor_land
                    pCell.aNeighbor_distance = pCell2.aNeighbor_distance

                    break

        return self.aCell

    def pyflowline_run(self):
        """
        Run the flowlinecase simulation

        Returns:
            list: A list of cell objects
        """
        aCell_out = None
        if self.iFlag_flowline == 1:
            self.pyflowline_flowline_simplification()
        else:
            pass

        if self.iFlag_create_mesh:
            self.aCell = self.pyflowline_mesh_generation(
                iFlag_antarctic_in=self.iFlag_antarctic, iFlag_arctic_in=self.iFlag_arctic)
            aCell_out = self.aCell
            pass
        else:
            # may be read mesh
            iMesh_type = self.iMesh_type
            # there must be some auxiliary file associated with the mesh file
            # self.aCell = read_mesh_json_w_topology(iMesh_type, self.sFilename_mesh)
            # aCell_out = self.aCell
            pass

        if self.iFlag_intersect == 1:
            self.aCell, a, b = self.pyflowline_reconstruct_topological_relationship()
            aCell_out = self.aCell
            pass
        else:
            pass

        return aCell_out

    def pyflowline_analyze(self):
        """
        Analyze the domain results for every watershed
        """
        print('PyFlowline started analyze results')
        ptimer = pytimer()
        ptimer.start()
        if self.iFlag_flowline == 1:
            if self.iFlag_analysis == 1:
                for pBasin in self.aBasin:
                    pBasin.basin_analyze()
                    pass
        ptimer.stop()
        return

    def pyflowline_evaluate(self):
        """
        Evaluate the model performance
        """
        print('PyFlowline started evaluate results')
        ptimer = pytimer()
        ptimer.start()
        if self.iFlag_evaluation == 1:
            for pBasin in self.aBasin:
                pBasin.basin_evaluate(self.iMesh_type, self.sMesh_type)
        ptimer.stop()
        return

    def pyflowline_export(self):
        """
        Export the model outputs
        """
        print('PyFlowline started export results')
        ptimer = pytimer()
        ptimer.start()
        self.pyflowline_export_mesh_info_to_json()
        # convert the mesh into the kml format so it can be visualized in google earth and google map
        # shoule move this to the export function
        if iFlag_kml is not None:
            # convert_geojson_to_kml(self.sFilename_mesh, self.sFilename_mesh_kml)
            pass

        if self.iFlag_flowline == 1:
            for pBasin in self.aBasin:
                pBasin.basin_export()

        self.tojson()
        ptimer.stop()
        return

    def pyflowline_export_mesh_info_to_json(self):
        """
        Export the mesh information to a json file
        """

        sFilename_json = self.sFilename_mesh_info

        # with open(sFilename_json, 'w', encoding='utf-8') as f:
        #    sJson = json.dumps([json.loads(ob.tojson()) for ob in self.aCell], indent = 4)
        #    f.write(sJson)
        #    f.close()

        # a more efficient way to handle large datasets
        with open(sFilename_json, 'w', encoding='utf-8') as f:
            f.write('[\n')
            for i, ob in enumerate(self.aCell):
                sJson = json.dumps(json.loads(ob.tojson()), indent=4)
                if i < len(self.aCell) - 1:
                    f.write(sJson + ',\n')
                else:
                    f.write(sJson + '\n')
            f.write(']\n')

        return

    def tojson(self):
        """
        Convert the flowline case object to a json string

        Returns:
            json str: A json string
        """
        aSkip = ['aBasin',
                 'aFlowline_simplified', 'aFlowline_conceptual', 'aCellID_outlet',
                 'aCell', 'pRTree_mesh']

        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)
            pass

        sJson = json.dumps(obj,
                           sort_keys=True,
                           indent=4,
                           ensure_ascii=True,
                           cls=CaseClassEncoder)
        return sJson

    def pyflowline_print(self):
        """
        Print the flowline case object
        """
        print(self.tojson())
        return

    def pyflowline_export_config_to_json(self,
                                         sFilename_output_in=None,
                                         iFlag_export_basin_in = 1):
        """
        Export the configuration to a json file

        Args:
            sFilename_output_in (str, optional): The json filename. Defaults to None.
        """

        if self.iFlag_standalone == 1:
            if sFilename_output_in is not None:
                sFilename_output = sFilename_output_in
                sPath_current = Path(sFilename_output).parent
            else:
                # use current output path
                sFilename_output = os.path.join(
                    self.sWorkspace_output, 'configuration.json')
                sPath_current = Path(self.sWorkspace_output)

            # all basins
            if iFlag_export_basin_in == 1:
                sName = 'configuration_basin.json'
                sFilename_configuration = os.path.join(
                    sPath_current, sName)
                with open(sFilename_configuration, 'w', encoding='utf-8') as f:
                    sJson = json.dumps([json.loads(ob.tojson()) for ob in self.aBasin],
                                       sort_keys=True,
                                       indent=4)
                    f.write(sJson)
                    f.close()
                    # update
                self.sFilename_basins = sFilename_configuration
        else:
            if sFilename_output_in is not None:
                sFilename_output = sFilename_output_in
                sPath_current = Path(sFilename_output).parent
            else:
                # use parent path
                sPath_current = Path(self.sWorkspace_output)
                sFilename_output = os.path.join(
                    sPath_current, 'configuration.json')

            # all basins
            if iFlag_export_basin_in == 1:
                sName = 'configuration_basin.json'
                sFilename_configuration = os.path.join(
                    sPath_current, sName)
                with open(sFilename_configuration, 'w', encoding='utf-8') as f:
                    sJson = json.dumps([json.loads(ob.tojson()) for ob in self.aBasin],
                                       sort_keys=True,
                                       indent=4)
                    f.write(sJson)
                    f.close()
                    # update for pyhexwatershed
                self.sFilename_basins = sFilename_configuration

            self.sWorkspace_output = sPath_current.parent.absolute()

        aSkip = ['aBasin',
                 'aFlowline_simplified', 'aFlowline_conceptual', 'aCellID_outlet',
                 'aCell', 'pRTree_mesh']

        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)

        with open(sFilename_output, 'w', encoding='utf-8') as f:
            json.dump(obj, f, sort_keys=True,
                      ensure_ascii=False,
                      indent=4,
                      cls=CaseClassEncoder)
        return

    def pyflowline_export_basin_config_to_json(self, sFilename_output_in=None):
        """
        Export the member basin configuration to a json file

        Args:
            sFilename_output_in (str, optional): The json filename. Defaults to None.
        """
        if self.iFlag_standalone == 1:
            if sFilename_output_in is not None:
                sFilename_output = sFilename_output_in
            else:
                # use current output path
                sName = 'configuration_basin.json'
                sFilename_output = os.path.join(self.sWorkspace_output, sName)

            # all basins
            with open(sFilename_output, 'w', encoding='utf-8') as f:
                sJson = json.dumps([json.loads(ob.tojson()) for ob in self.aBasin],
                                   sort_keys=True,
                                   indent=4)
                f.write(sJson)
                f.close()
                # update
            self.sFilename_basins = sFilename_output

        else:
            if sFilename_output_in is not None:
                sFilename_output = sFilename_output_in
            else:
                # use current output path
                sPath = Path(self.sWorkspace_output)
                sName = 'configuration_basin.json'
                sFilename_output = os.path.join(sPath.parent.absolute(), sName)

            # all basins
            with open(sFilename_output, 'w', encoding='utf-8') as f:
                sJson = json.dumps([json.loads(ob.tojson()) for ob in self.aBasin],
                                   sort_keys=True,
                                   indent=4)
                f.write(sJson)
                f.close()
                # update for pyhexwatershed
            self.sFilename_basins = sFilename_output

        return
