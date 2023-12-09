import os
import stat
from pathlib import Path
import json
from json import JSONEncoder
import datetime
import importlib
from shutil import copy2
import numpy as np
from osgeo import ogr, osr, gdal

from pyearth.system.define_global_variables import *
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


from pyearth.gis.spatialref.reproject_coodinates import reproject_coordinates
from pyearth.gis.spatialref.conversion_between_degree_and_meter  import degree_to_meter
from pyearth.gis.spatialref.conversion_between_degree_and_meter  import meter_to_degree
from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import gdal_read_geotiff_file
from pyearth.gis.gdal.read.vector.gdal_read_geojson_boundary import gdal_read_geojson_boundary
from pyearth.gis.spatialref.get_utm_spatial_reference import get_utm_spatial_reference
from pyearth.toolbox.data.geoparquet.convert_geojson_to_geoparquet import convert_geojson_to_geoparquet

iFlag_kml = importlib.util.find_spec("simplekml") 


iFlag_cython = importlib.util.find_spec("cython") 
if iFlag_cython is not None:
    from pyflowline.external.tinyr.tinyr.tinyr import RTree
    iFlag_use_rtree = 1
else:
    iFlag_use_rtree =0
    pass

from pyflowline.mesh.hexagon.create_hexagon_mesh import create_hexagon_mesh
from pyflowline.mesh.latlon.create_latlon_mesh import create_latlon_mesh
from pyflowline.mesh.square.create_square_mesh import create_square_mesh
from pyflowline.mesh.mpas.create_mpas_mesh import create_mpas_mesh
from pyflowline.mesh.dggrid.create_dggrid_mesh import create_dggrid_mesh
from pyflowline.mesh.dggrid.create_dggrid_mesh import dggrid_find_resolution_by_index
from pyflowline.mesh.tin.create_tin_mesh import create_tin_mesh
gdal.UseExceptions()   
pDate = datetime.datetime.today()
sDate_default = "{:04d}".format(pDate.year) + "{:02d}".format(pDate.month) + "{:02d}".format(pDate.day)

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
    iCase_index= 0

    sMesh_type = 1
    iFlag_standalone=1
    iFlag_flowline = 1
    iFlag_global = 0
    iFlag_antarctic = 0
    iFlag_multiple_outlet = 0
    iFlag_mesh_boundary = 0
    #iFlag_use_shapefile_extent=1
    iFlag_use_mesh_dem=0
    iFlag_save_mesh = 0
    iFlag_simplification = 1 #user can turn on/off
    iFlag_create_mesh=1
    iFlag_intersect = 1
    iFlag_rotation=0
    iFlag_break_by_distance = 0  #if the distance between two vertice are far,
    iResolution_index = 10
    iFlag_dggrid = 0
    nOutlet = 1 #by default , there shoule ne only one ouelet
    dResolution_degree=0.0
    dResolution_meter=0.0
    dThreshold_small_river=0.0
    dThreshold_break_by_distance = 5000.0
    dLongitude_left = -180
    dLongitude_right = 180
    dLatitude_bot = -90
    dLatitude_top = 90
    dElevation_mean=-9999.0
    sFilename_model_configuration=''
    sFilename_mesh_boundary = ''
    sWorkspace_input=''
    sWorkspace_output=''
    #sWorkspace_output_case=''
    sRegion=''
    sModel=''
    sMesh_type ='hexagon'
    sDggrid_type = 'ISEA3H'

    sCase=''
    sDate=''

    sFilename_dggrid=''

    sFilename_spatial_reference=''
    sFilename_dem=''

    sFilename_mesh=''
    sFilename_mesh_info=''
    sFilename_mesh_netcdf=''
    sFilename_mesh_kml=''

    aBasin = list()
    aFlowline_simplified=list()
    aFlowline_conceptual=list()
    aCellID_outlet = list()
    aCell=list()

    pRTree_mesh = None #only when rtree is used

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
                 iFlag_standalone_in= None,
                 sModel_in = None,
                 sDate_in = None,
                 sWorkspace_output_in = None):
        """
        Initialize a flowlinecase object

        Args:
            aConfig_in (dict): A dictionary of parameters
            iFlag_standalone_in (int, optional): Flag for whether run the case standalone. Defaults to None.
            sModel_in (str, optional): The model name. Defaults to None.
            sDate_in (str, optional): The case date. Defaults to None.
            sWorkspace_output_in (str, optional): The output workspace. Defaults to None.
        """
        #flags
        if iFlag_standalone_in is not None:
            self.iFlag_standalone = iFlag_standalone_in
        else:
            if 'iFlag_standalone' in aConfig_in:
                self.iFlag_standalone = int(aConfig_in['iFlag_standalone'])
            else:
                self.iFlag_standalone=1


        if 'iFlag_flowline' in aConfig_in:
            self.iFlag_flowline             = int(aConfig_in[ 'iFlag_flowline'])
        else:
            #without iFlag_flowline the model is a mesh generator
            self.iFlag_flowline=1

        if 'iFlag_global' in aConfig_in:
            self.iFlag_global             = int(aConfig_in[ 'iFlag_global'])

        if 'iFlag_antarctic' in aConfig_in:
            self.iFlag_antarctic             = int(aConfig_in[ 'iFlag_antarctic'])

        if 'iFlag_mesh_boundary' in aConfig_in:
            self.iFlag_mesh_boundary             = int(aConfig_in[ 'iFlag_mesh_boundary'])
        else:
            self.iFlag_mesh_boundary=0

        if 'iFlag_multiple_outlet' in aConfig_in:
            self.iFlag_multiple_outlet             = int(aConfig_in[ 'iFlag_multiple_outlet'])

        if 'iFlag_simplification' in aConfig_in:
            self.iFlag_simplification = int(aConfig_in['iFlag_simplification'])
            #if simplification is desired, then we must activate flowline flag
            #it can be turn off only when a previous simplification was done and you dont want to redo it
        else:
            self.iFlag_simplification  = 1
        if self.iFlag_simplification ==1:
            self.iFlag_flowline = 1

        if 'iFlag_create_mesh' in aConfig_in:
            self.iFlag_create_mesh = int(aConfig_in['iFlag_create_mesh'])
        else:
            self.iFlag_create_mesh=1

        if 'iFlag_save_mesh' in aConfig_in:
            self.iFlag_save_mesh             = int(aConfig_in[ 'iFlag_save_mesh'])

        if 'iFlag_use_mesh_dem' in aConfig_in:
            self.iFlag_use_mesh_dem = int(aConfig_in['iFlag_use_mesh_dem'])

        #if 'iFlag_use_shapefile_extent' in aConfig_in:
        #    self.iFlag_use_shapefile_extent = int(aConfig_in['iFlag_use_shapefile_extent'])

        if 'iFlag_rotation' in aConfig_in:
            self.iFlag_rotation = int(aConfig_in['iFlag_rotation'])


        if 'iFlag_intersect' in aConfig_in:
            self.iFlag_intersect = int(aConfig_in['iFlag_intersect'])
        else:
            self.iFlag_intersect=1
        
        if self.iFlag_flowline == 1:
            pass
        else:
            self.iFlag_intersect=0

        if 'iFlag_break_by_distance' in aConfig_in:
            self.iFlag_break_by_distance = int(aConfig_in['iFlag_break_by_distance'])
        else:
            self.iFlag_break_by_distance=0

        if 'iFlag_dggrid' in aConfig_in:
            self.iFlag_dggrid = int(aConfig_in['iFlag_dggrid'])
        else:
            self.iFlag_dggrid=0

        if 'iResolution_index' in aConfig_in:
            self.iResolution_index = int(aConfig_in['iResolution_index'])
        else:
            self.iResolution_index=10
        
        if 'sDggrid_type' in aConfig_in:
            self.sDggrid_type = aConfig_in['sDggrid_type']
        else:
            self.sDggrid_type='ISEA3H'

        if 'nOutlet' in aConfig_in:
            self.nOutlet = int(aConfig_in['nOutlet'])

        if 'iCase_index' in aConfig_in:
            iCase_index = int(aConfig_in['iCase_index'])
        else:
            iCase_index = 1
        
        sCase_index = "{:03d}".format( iCase_index )
        self.iCase_index =   iCase_index

        if 'dResolution_degree' in aConfig_in:
            self.dResolution_degree = float(aConfig_in['dResolution_degree'])

        if 'dResolution_meter' in aConfig_in:
            self.dResolution_meter = float(aConfig_in['dResolution_meter'])
        else:
            print('Please specify resolution.')

        

        if 'dThreshold_break_by_distance' in aConfig_in:
            self.dThreshold_break_by_distance = float(aConfig_in['dThreshold_break_by_distance'])

        if 'dLongitude_left' in aConfig_in:
            self.dLongitude_left = float(aConfig_in['dLongitude_left'])

        if 'dLongitude_right' in aConfig_in:
            self.dLongitude_right = float(aConfig_in['dLongitude_right'])

        if 'dLatitude_bot' in aConfig_in:
            self.dLatitude_bot = float(aConfig_in['dLatitude_bot'])

        if 'dLatitude_top' in aConfig_in:
            self.dLatitude_top = float(aConfig_in['dLatitude_top'])

        if 'sFilename_model_configuration' in aConfig_in:
            self.sFilename_model_configuration    = aConfig_in[ 'sFilename_model_configuration']

        if 'sFilename_mesh_boundary' in aConfig_in:
            self.sFilename_mesh_boundary    = aConfig_in[ 'sFilename_mesh_boundary']

            if self.iFlag_mesh_boundary==1:
                if not os.path.isfile(self.sFilename_mesh_boundary ):
                    print("The mesh boundary file does not exist, you should update this parameter before running the model!")
                    exit()
                pass

        if 'sFilename_spatial_reference' in aConfig_in:
            self.sFilename_spatial_reference = aConfig_in['sFilename_spatial_reference']


        if 'sFilename_dem' in aConfig_in:
            self.sFilename_dem = aConfig_in['sFilename_dem']
        
        if self.iFlag_dggrid == 1:
            if 'sFilename_dggrid' in aConfig_in:
                self.sFilename_dggrid = aConfig_in['sFilename_dggrid']
                if not os.path.isfile(self.sFilename_dggrid ):
                    print("The dggrid binary file does not exist, you need to update this parameter before running the model!")
                    exit()
                pass

            else:
                print('Please specify the dggrid binary file.')
                exit()

        if 'sFilename_mesh_netcdf' in aConfig_in:
            self.sFilename_mesh_netcdf = aConfig_in['sFilename_mesh_netcdf']

        if 'sWorkspace_bin' in aConfig_in:
            self.sWorkspace_bin= aConfig_in[ 'sWorkspace_bin']

        if 'sWorkspace_input' in aConfig_in:
            self.sWorkspace_input = aConfig_in[ 'sWorkspace_input']


        if sWorkspace_output_in is not None:
            self.sWorkspace_output = sWorkspace_output_in
        else:
            if 'sWorkspace_output' in aConfig_in:
                self.sWorkspace_output = aConfig_in[ 'sWorkspace_output']

        if 'sJob' in aConfig_in:
            self.sJob =  aConfig_in['sJob']
        if 'sRegion' in aConfig_in:
            self.sRegion               = aConfig_in[ 'sRegion']

        if sModel_in is not None:
            self.sModel = sModel_in
        else:
            if 'sModel' in aConfig_in:
                self.sModel                = aConfig_in[ 'sModel']


        sDate   = aConfig_in[ 'sDate']
        if sDate is not None:
            self.sDate= sDate
        else:
            self.sDate = sDate_default


        sCase = self.sModel  + self.sDate + sCase_index
        self.sCase = sCase

        if 'sMesh_type' in aConfig_in:
            self.sMesh_type =  aConfig_in['sMesh_type']
        else:
            self.sMesh_type = 'hexagon'



        sMesh_type = self.sMesh_type
        if sMesh_type =='hexagon': #hexagon
            self.iMesh_type = 1
        else:
            if sMesh_type =='square': #square
                self.iMesh_type = 2
            else:
                if sMesh_type =='latlon': #latlon
                    self.iMesh_type = 3
                else:
                    if sMesh_type =='mpas': #mpas
                        self.iMesh_type = 4
                    else:
                        if sMesh_type =='dggrid':
                            self.iMesh_type = 5
                        else:
                            if sMesh_type =='tin': #
                                self.iMesh_type = 6
                            else:
                                print('Unsupported mesh type?')

        if self.iMesh_type == 5:
            #update resolution
            self.dResolution_meter = dggrid_find_resolution_by_index(self.sDggrid_type, self.iResolution_index)
            
        #the model can be run as part of hexwatershed or standalone
        if self.iFlag_standalone == 1:
            #in standalone case, will add case information and update output path
            sPath = str(Path(self.sWorkspace_output)  /  sCase)
            self.sWorkspace_output = sPath
        else:
            #use specified output path, also do not add output or input tag
            sPath = self.sWorkspace_output

        Path(sPath).mkdir(parents=True, exist_ok=True)

        if self.iMesh_type == 4:
            if not os.path.isfile(self.sFilename_mesh_netcdf ):
                print("The MPAS mesh file does not exist!")
                exit()
        else:
            if not os.path.isfile(self.sFilename_dem ): #why DEM is required?
                print("The DEM file does not exist!")
                #exit()

        self.aBasin = list()

        if self.iFlag_flowline==1:
            if 'sFilename_basins' in aConfig_in:
                self.sFilename_basins = aConfig_in['sFilename_basins']
                if os.path.isfile(self.sFilename_basins):
                    pass
                else:
                    print('This basin configuration file does not exist: ', self.sFilename_basins )
                    print('Please update this parameter before running the model!')
                    exit()
                    
                with open(self.sFilename_basins) as json_file:
                    dummy_data = json.load(json_file)
                    for i in range(self.nOutlet):
                        sBasin =  "{:08d}".format(i+1)
                        dummy_basin = dummy_data[i]
                        dummy_basin['sWorkspace_output_basin'] = str(Path(self.sWorkspace_output) / sBasin )
                        pBasin = pybasin(dummy_basin)
                        self.aBasin.append(pBasin)
            else:
                pass
        else:
            
            pass


        #model generated files
        
        self.sFilename_mesh = os.path.join(str(Path(self.sWorkspace_output)  ) , sMesh_type + ".geojson" )
        self.sFilename_mesh_info= os.path.join(str(Path(self.sWorkspace_output)  ) , sMesh_type + "_mesh_info.json"  )
        self.sFilename_mesh_kml = os.path.join(str(Path(self.sWorkspace_output)  ) , sMesh_type + ".kml" ) #for google service

        return


    def pyflowline_convert_flowline_to_geojson(self):
        if self.iFlag_flowline == 1:
            for pBasin in self.aBasin:
                pBasin.basin_convert_flowline_to_geojson()
                pass

        return

    def pyflowline_flowline_simplification(self):
        aFlowline_out = list()   #store all the flowline
        #export the outlet into a single file
        aOutlet = list()
        if self.iFlag_simplification == 1:
            for i in range(self.nOutlet):
                pBasin = self.aBasin[i]
                aFlowline_basin = pBasin.basin_flowline_simplification()
                aFlowline_out = aFlowline_out + aFlowline_basin
                aOutlet.append(pBasin.pVertex_outlet)

            self.aFlowline_simplified = aFlowline_out


        return aFlowline_out

    def pyflowline_mesh_generation(self, iFlag_antarctic_in=None):
        """
        The mesh generation operation

        Returns:
            list [pycell]: A list of cell object
        """
        print('Start mesh generation.')
        if iFlag_antarctic_in is None:
            iFlag_antarctic = 0
        else:
            iFlag_antarctic = iFlag_antarctic_in

     
        if self.iFlag_create_mesh ==1:
            iFlag_global =  self.iFlag_global
            iMesh_type = self.iMesh_type
            iFlag_save_mesh = self.iFlag_save_mesh
            iFlag_rotation = self.iFlag_rotation
            iFlag_mesh_boundary = self.iFlag_mesh_boundary
            dResolution_degree = self.dResolution_degree
            dResolution_meter = self.dResolution_meter
            sFilename_dem = self.sFilename_dem
            sFilename_spatial_reference = self.sFilename_spatial_reference
            sFilename_mesh = self.sFilename_mesh 

            if iFlag_global ==1: #a global mesh does not require boundary
                pass
            else:

                if iFlag_mesh_boundary ==1:
                    #create a polygon based on real boundary
                    pBoundary_wkt, aExtent = gdal_read_geojson_boundary(self.sFilename_mesh_boundary)           

                    if iMesh_type != 4: #not mpas
                        spatial_reference_target = osr.SpatialReference()
                        spatial_reference_target.ImportFromEPSG(4326)

                        #check whether DEM exists
                        if os.path.isfile(sFilename_dem):
                            #dPixelWidth, pPixelHeight, dOriginX, dOriginY, nrow, ncolumn, pSpatialRef_dem, pProjection, pGeotransform\
                            #    = retrieve_geotiff_metadata(sFilename_dem)
                            dummy = gdal_read_geotiff_file(sFilename_dem, iFlag_metadata_only= 1)
                            dPixelWidth = dummy['pixelWidth']                        
                            pPixelHeight = dummy['pixelHeight']
                            dOriginX = dummy['originX']
                            dOriginY = dummy['originY']
                            nrow = dummy['nrow']
                            ncolumn = dummy['ncolumn']
                            pSpatialRef_dem = dummy['spatialReference']
                            #pProjection= dummy['projection']
                            #pGeotransform = dummy['geotransform']

                            #lower left
                            dX_lowerleft  = dOriginX
                            dY_lowerleft = dOriginY + (nrow+1) * pPixelHeight
                            dLongitude_left0,  dLatitude_bot0= reproject_coordinates(dX_lowerleft, dY_lowerleft,pSpatialRef_dem,    spatial_reference_target)

                            #upper right
                            dX_upperright = dOriginX + (ncolumn +1) * dPixelWidth
                            dY_upperright = dOriginY
                            dLongitude_right0, dLatitude_top0= reproject_coordinates(dX_upperright, dY_upperright,pSpatialRef_dem,  spatial_reference_target)

                            #lower right
                            dX_lowerright = dOriginX + (ncolumn +1) * dPixelWidth
                            dY_lowerright = dOriginY - (nrow+1) * dPixelWidth

                            dLongitude_right1,  dLatitude_bot1= reproject_coordinates(dX_lowerright, dY_lowerright,pSpatialRef_dem,    spatial_reference_target)

                            #uppler left
                            dX_upperleft   = dOriginX
                            dY_upperleft   =  dOriginY
                            dLongitude_left1, dLatitude_top1= reproject_coordinates(dX_upperleft, dY_upperleft, pSpatialRef_dem,  spatial_reference_target)

                            dLatitude_top = np.max([dLatitude_top0, dLatitude_top1])
                            dLatitude_bot = np.min([dLatitude_bot0, dLatitude_bot1])

                            dLongitude_left = np.min([dLongitude_left0, dLongitude_left1])
                            dLongitude_right = np.max([dLongitude_right0, dLongitude_right1])

                            dLatitude_mean = 0.5 * (dLatitude_top + dLatitude_bot)
                            #pass

                            if dResolution_meter is None:
                                #not used
                                pass
                            else:
                                dResolution_degree = meter_to_degree(dResolution_meter, dLatitude_mean)

                            dX_lowerleft = dOriginX
                            dY_upperleft = dOriginY
                        else: #use the spatial extent and utm zone to obatin the mesh boundary
                            dLongitude_left, dLongitude_right,dLatitude_bot,  dLatitude_top = aExtent
                            dLatitude_mean = 0.5 * (dLatitude_top + dLatitude_bot)
                            dLongitude_mean = 0.5 * (dLongitude_left + dLongitude_right)

                            if iMesh_type == 3: #latlon

                                pass
                            else:
                                if dResolution_meter is None:
                                    #not used
                                    pass
                                else:
                                    dResolution_degree = meter_to_degree(dResolution_meter, dLatitude_mean)
                                pass

                            pUTM = get_utm_spatial_reference(dLongitude_mean)

                            #calculate the lower left and upper right in utm projection
                            dX_lowerleft, dY_lowerleft = reproject_coordinates(dLongitude_left, dLatitude_bot, spatial_reference_target,  pUTM)

                            dX_lowerright, dY_lowerright = reproject_coordinates(dLongitude_right, dLatitude_bot, spatial_reference_target,  pUTM)

                            dX_upperright, dY_upperright = reproject_coordinates(dLongitude_right, dLatitude_top, spatial_reference_target,  pUTM)

                            dX_upperleft, dY_upperright = reproject_coordinates(dLongitude_left, dLatitude_top, spatial_reference_target,  pUTM)

                    else:

                        pass
                else:
                    #if the boundary is not specified, a DEM file is required
                    if iMesh_type != 4: #mpas
                        spatial_reference_target = osr.SpatialReference()
                        spatial_reference_target.ImportFromEPSG(4326)

                        #check whether DEM exists
                        if os.path.isfile(sFilename_dem):
                            dummy = gdal_read_geotiff_file(sFilename_dem, iFlag_metadata_only= 1)
                            dPixelWidth = dummy['pixelWidth']                        
                            pPixelHeight = dummy['pixelHeight']
                            dOriginX = dummy['originX']
                            dOriginY = dummy['originY']
                            nrow = dummy['nrow']
                            ncolumn = dummy['ncolumn']
                            pSpatialRef_dem = dummy['spatialReference']
                            pProjection= dummy['projection']
                            pGeotransform = dummy['geotransform']

                            #lower left
                            dX_lowerleft  = dOriginX
                            dY_lowerleft = dOriginY - (nrow+1) * dPixelWidth
                            dLongitude_left0,  dLatitude_bot0= reproject_coordinates(dX_lowerleft, dY_lowerleft,pSpatialRef_dem,    spatial_reference_target)

                            #upper right
                            dX_upperright = dOriginX + (ncolumn +1) * dPixelWidth
                            dY_upperright = dOriginY
                            dLongitude_right0, dLatitude_top0= reproject_coordinates(dX_upperright, dY_upperright,pSpatialRef_dem,  spatial_reference_target)

                            #lower right
                            dX_lowerright = dOriginX + (ncolumn +1) * dPixelWidth
                            dY_lowerright = dOriginY - (nrow+1) * dPixelWidth

                            dLongitude_right1,  dLatitude_bot1= reproject_coordinates(dX_lowerright, dY_lowerright,pSpatialRef_dem,    spatial_reference_target)

                            #uppler left
                            dX_upperleft   = dOriginX
                            dY_upperleft   =  dOriginY
                            dLongitude_left1, dLatitude_top1= reproject_coordinates(dX_upperleft, dY_upperleft, pSpatialRef_dem,  spatial_reference_target)

                            dLatitude_top = np.max([dLatitude_top0, dLatitude_top1])
                            dLatitude_bot = np.min([dLatitude_bot0, dLatitude_bot1])

                            dLongitude_left = np.min([dLongitude_left0, dLongitude_left1])
                            dLongitude_right = np.max([dLongitude_right0, dLongitude_right1])

                            dLatitude_mean = 0.5 * (dLatitude_top + dLatitude_bot)
                            #pass

                            if dResolution_meter < 0:
                                #not used
                                pass
                            else:
                                dResolution_degree = meter_to_degree(dResolution_meter, dLatitude_mean)

                            dX_lowerleft = dOriginX
                            dY_upperleft = dOriginY
                        else:
                            print("Error: DEM file does not exist: " + sFilename_dem)
                            return
                    pass

            if iMesh_type ==1: #hexagon
                #hexagon edge
                dResolution_meter = degree_to_meter(dLatitude_mean, dResolution_degree )
                dArea = np.power(dResolution_meter,2.0)
                dLength_edge = np.sqrt(  2.0 * dArea / (3.0* np.sqrt(3.0))  )
                if iFlag_rotation ==0:
                    dX_spacing = dLength_edge * np.sqrt(3.0)
                    dY_spacing = dLength_edge * 1.5
                    ncolumn= int( (dX_lowerright - dX_lowerleft) / dX_spacing )
                    nrow= int( (dY_upperleft - dY_lowerleft) / dY_spacing )
                else:
                    dX_spacing = dLength_edge * 1.5
                    dY_spacing = dLength_edge * np.sqrt(3.0)
                    ncolumn= int( (dX_lowerright - dX_lowerleft) / dX_spacing )+1
                    nrow= int( (dY_upperleft - dY_lowerleft) / dY_spacing )

                if iFlag_mesh_boundary ==1:
                    #create a polygon based on real boundary
                    pBoundary_wkt, aExtent = gdal_read_geojson_boundary(self.sFilename_mesh_boundary)

                    aHexagon = create_hexagon_mesh(iFlag_rotation, dX_lowerleft, dY_lowerleft, dResolution_meter, ncolumn, nrow, 
                                                   sFilename_mesh, sFilename_spatial_reference, pBoundary_wkt)
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
                    pBoundary_wkt = pBoundary.ExportToWkt() #wkt format

                    aHexagon = create_hexagon_mesh(iFlag_rotation, dX_lowerleft, dY_lowerleft, dResolution_meter, ncolumn, nrow,
                                                   sFilename_mesh, sFilename_spatial_reference, pBoundary_wkt)
                    
                    
                    pass

                

                self.aCell = aHexagon
            else:
                if iMesh_type ==2: #sqaure
                    ncolumn= int( (dX_lowerright - dX_lowerleft) / dResolution_meter )
                    nrow= int( (dY_upperleft - dY_lowerleft) / dResolution_meter )
                    if iFlag_mesh_boundary ==1:
                        #create a polygon based on real boundary
                        pBoundary_wkt, aExtent= gdal_read_geojson_boundary(self.sFilename_mesh_boundary)

                        aSquare = create_square_mesh(dX_lowerleft, dY_lowerleft, dResolution_meter, ncolumn, nrow, 
                                                     sFilename_mesh, sFilename_spatial_reference, pBoundary_wkt)
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
                        aSquare = create_square_mesh(dX_lowerleft, dY_lowerleft, dResolution_meter, ncolumn, nrow, 
                                                     sFilename_mesh, sFilename_spatial_reference, pBoundary_wkt)
                        pass
                    
                    self.aCell = aSquare
                else:
                    if iMesh_type ==3: #latlon
                        dResolution_meter = degree_to_meter(dLatitude_mean, dResolution_degree)
                        dArea = np.power(dResolution_meter,2.0)
                        ncolumn= int( (dLongitude_right - dLongitude_left) / dResolution_degree )
                        nrow= int( (dLatitude_top - dLatitude_bot) / dResolution_degree )

                        if iFlag_mesh_boundary ==1:
                            #create a polygon based on real boundary
                            #already produced
                            pBoundary_wkt , aExtent= gdal_read_geojson_boundary(self.sFilename_mesh_boundary)
                            aLatlon = create_latlon_mesh(dLongitude_left, dLatitude_bot, dResolution_degree, ncolumn, nrow, 
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
                            pBoundary_wkt =  pBoundary.ExportToWkt() 
                            aLatlon = create_latlon_mesh(dLongitude_left, dLatitude_bot, dResolution_degree, ncolumn, nrow, \
                                                         sFilename_mesh, pBoundary_wkt)

                            pass
                    
                        self.aCell = aLatlon

                        
                    else:
                        if iMesh_type == 4: #mpas
                            iFlag_use_mesh_dem = self.iFlag_use_mesh_dem
                            sFilename_mesh_netcdf = self.sFilename_mesh_netcdf
                            dLatitude_top    = self.dLatitude_top
                            dLatitude_bot    = self.dLatitude_bot
                            dLongitude_left  = self.dLongitude_left
                            dLongitude_right = self.dLongitude_right

                            if iFlag_antarctic ==1:
                                aMpas = create_mpas_mesh(iFlag_global, iFlag_use_mesh_dem, iFlag_save_mesh,
                                                         sFilename_mesh_netcdf, sFilename_mesh, iFlag_antarctic_in=iFlag_antarctic_in )
                                pass
                            else:

                                if iFlag_mesh_boundary ==1:
                                    #create a polygon based on
                                    #read boundary
                                    pBoundary_wkt, aExtent = gdal_read_geojson_boundary(self.sFilename_mesh_boundary)

                                    aMpas = create_mpas_mesh(iFlag_global, iFlag_use_mesh_dem, iFlag_save_mesh,
                                                             sFilename_mesh_netcdf,  sFilename_mesh, iFlag_antarctic_in=iFlag_antarctic_in, pBoundary_in = pBoundary_wkt)
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
                                    pBoundary_wkt =  pBoundary.ExportToWkt() 

                                    #new method using polygon object
                                    aMpas = create_mpas_mesh(iFlag_global, iFlag_use_mesh_dem, iFlag_save_mesh, 
                                                             sFilename_mesh_netcdf, sFilename_mesh, iFlag_antarctic_in= iFlag_antarctic_in, pBoundary_in = pBoundary_wkt  )
                                    pass

                            self.aCell = aMpas
                        else:
                            if iMesh_type == 5: #dggrid
                                dLatitude_top    = self.dLatitude_top
                                dLatitude_bot    = self.dLatitude_bot
                                dLongitude_left  = self.dLongitude_left
                                dLongitude_right = self.dLongitude_right                                
                                

                                if self.iFlag_standalone == 1:
                                    sWorkspace_output = self.sWorkspace_output + slash + 'dggrid'
                                    pass
                                else:
                                    sWorkspace_output = self.sWorkspace_output + slash + '..'+ slash + 'dggrid'
                                    
                                    sWorkspace_output = os.path.abspath(sWorkspace_output)
                                    pass

                                if not os.path.exists(sWorkspace_output):
                                    os.makedirs(sWorkspace_output)
                                
                                if iFlag_mesh_boundary ==1:
                                        #create a polygon based on                                   

                                    aDggrid = create_dggrid_mesh(iFlag_global,
                                                                     iFlag_save_mesh,
                                                                     sFilename_mesh,
                                                                     sWorkspace_output,
                                                                     iResolution_index_in= self.iResolution_index,
                                                                     iFlag_antarctic_in=iFlag_antarctic_in,
                                                                     sDggrid_type_in = self.sDggrid_type,
                                                                     sFilename_boundary_in = self.sFilename_mesh_boundary)
                                    
                                    
                                                                     
                                    pass
                                else:                                       
                                    aDggrid = create_dggrid_mesh(iFlag_global,
                                                                     iFlag_save_mesh,                                                                   
                                                                     sFilename_mesh,
                                                                     sWorkspace_output,
                                                                      iResolution_index_in= self.iResolution_index,
                                                                          sDggrid_type_in = self.sDggrid_type)
                                    
                                    pass


                                self.aCell = aDggrid

                            else:
                                if iMesh_type == 6: #tin this one need to be updated because central location issue

                                    #tin edge
                                    dArea = np.power(dResolution_meter,2.0)
                                    dLength_edge = np.sqrt( 4.0 * dArea / np.sqrt(3.0) )
                                    dX_shift = 0.5 * dLength_edge
                                    dY_shift = 0.5 * dLength_edge * np.sqrt(3.0)
                                    dX_spacing = dX_shift * 2
                                    dY_spacing = dY_shift
                                    ncolumn= int( (dX_lowerright - dX_lowerleft) / dX_shift )
                                    nrow= int( (dY_upperleft - dY_lowerleft) / dY_spacing )
                                    aTin = create_tin_mesh(dX_lowerleft, dY_lowerleft, dResolution_meter, ncolumn, nrow, sFilename_mesh, sFilename_spatial_reference)
                                    self.aCell = aTin
                                else:

                                    print('Unsupported mesh type?')
                                return

            #no matter what type of mash, we will convert it to geoparquet for easy visualization
            convert_geojson_to_geoparquet(sFilename_mesh, sFilename_mesh.replace('.geojson','.parquet'))    
        else:            
            pass

        #build rtree
        if iFlag_use_rtree == 1:
            interleaved = True
            self.pRTree_mesh = RTree(interleaved=interleaved, max_cap=5, min_cap=2)
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
            sFilename_mesh=self.sFilename_mesh         
            aFlowline_conceptual = list()   #store all the flowline
            aCellID_outlet = list()
            aBasin = list()
            aCell_intersect=list()
            ncell=len(self.aCell)
            #there is a one-to-one match between cell id and cell center because each cell has only one center         
            
            
            for pBasin in self.aBasin:
                aCell_intersect_basin = pBasin.basin_reconstruct_topological_relationship(iMesh_type,sFilename_mesh)
                aFlowline_conceptual = aFlowline_conceptual + pBasin.aFlowline_basin_conceptual
                aBasin.append(pBasin)
                aCellID_outlet.append(pBasin.lCellID_outlet)
                aCell_intersect = aCell_intersect + aCell_intersect_basin

                if iFlag_use_rtree == 1:
                    #use rtree to update topology and length
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
                    #set topology here               
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
                                                if lCellID ==  pBasin.lCellID_outlet:
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

                    #update length using rtree            
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
                    pCell.nNeighbor= pCell2.nNeighbor
                    pCell.aNeighbor_land= pCell2.aNeighbor_land
                    pCell.nNeighbor_land= pCell2.nNeighbor_land
                    pCell.aNeighbor_distance= pCell2.aNeighbor_distance

                    break

        return self.aCell

    def pyflowline_analyze(self):
        """
        Analyze the domain results for every watershed
        """
        if self.iFlag_flowline == 1:
            for pBasin in self.aBasin:
                pBasin.analyze()
        return

    def pyflowline_setup(self):
        """
        Set up the flowlinecase
        """
        if self.iFlag_flowline == 1:
            self.pyflowline_convert_flowline_to_geojson()
            pass
        
        if self.iFlag_dggrid == 1:
            #create dggrid output folder
            if self.iFlag_standalone == 1:
                sWorkspace_output = self.sWorkspace_output + slash + 'dggrid'
            else:
                sWorkspace_output = self.sWorkspace_output + slash + '..'+ slash + 'dggrid'
                sWorkspace_output = os.path.abspath(sWorkspace_output)

         
            Path(sWorkspace_output).mkdir(parents=True, exist_ok=True)
            #then copy the binary file to the folder
            #copy execulate
            sFilename_new = sWorkspace_output + slash + 'dggrid'
            copy2(self.sFilename_dggrid, sFilename_new)
            os.chmod(sFilename_new, stat.S_IREAD | stat.S_IWRITE | stat.S_IXUSR)
            pass
        return

    def pyflowline_change_model_parameter(self, sVariable_in, dValue, iFlag_basin_in = None):
        if iFlag_basin_in is None:
            if hasattr(self, sVariable_in):
                #get default data type                
                sType_default = type(getattr(self, sVariable_in))
                #get the data type of the input value
                sType_input = type(dValue)
                if sType_default == sType_input:      
                    setattr(self, sVariable_in, dValue)                    
                    pass   
                else:
                    print('Incorrect data type for the input value: ' + sVariable_in)
                return True
            else:
                print("This model parameter is unknown, please check the full parameter list in the documentation: " + sVariable_in)
                return False
            
            

        else:
            #this mean the variable is in the basin object
            for pBasin in self.aBasin:
                if hasattr(pBasin, sVariable_in):
                    #get default data type
                    sType_default = type(getattr(pBasin, sVariable_in))
                    sType_input = type(dValue)
                    if sType_default == sType_input:      
                        setattr(pBasin, sVariable_in, dValue)
                    else:
                        print('Incorrect data type for the input value: ' + sVariable_in)                       
                        return False
                else:
                    print("This model parameter is unknown, please check the full parameter list in the documentation: " + sVariable_in)
                    return False

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
            self.aCell = self.pyflowline_mesh_generation(iFlag_antarctic_in= self.iFlag_antarctic)
            aCell_out = self.aCell
            pass
        else:
            #may be read mesh           
            iMesh_type = self.iMesh_type
            #there must be some auxiliary file associated with the mesh file
            self.aCell = read_mesh_json_w_topology(iMesh_type, self.sFilename_mesh)             
            aCell_out = self.aCell
            pass

        if self.iFlag_intersect == 1:
            self.aCell, a, b = self.pyflowline_reconstruct_topological_relationship()
            aCell_out = self.aCell
            pass
        else:
            pass

        return aCell_out

    def pyflowline_evaluate(self):
        """
        Evaluate the model performance
        """
        for pBasin in self.aBasin:
            pBasin.evaluate(self.iMesh_type, self.sMesh_type)
        return

    def pyflowline_export(self):
        """
        Export the model outputs
        """
        print('Started export elevation')
        ptimer = pytimer()
        ptimer.start()
        self.pyflowline_export_mesh_info_to_json()
        #convert the mesh into the kml format so it can be visualized in google earth and google map
        #shoule move this to the export function
        if iFlag_kml is not None:
            #convert_geojson_to_kml(self.sFilename_mesh, self.sFilename_mesh_kml)
            pass

        if self.iFlag_flowline ==1:
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
        
        #if iFlag_flowline == 1: #if there is conceptual flowline
        #    pass
        #else:
        #    #only mesh, no flowline
        #    pass


        with open(sFilename_json, 'w', encoding='utf-8') as f:
            sJson = json.dumps([json.loads(ob.tojson()) for ob in self.aCell], indent = 4)
            f.write(sJson)
            f.close()

        return

    def tojson(self):
        """
        Convert the flowline case object to a json string

        Returns:
            json str: A json string
        """
        aSkip = ['aBasin', \
                 'aFlowline_simplified','aFlowline_conceptual','aCellID_outlet',
                 'aCell', 'pRTree_mesh']

        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)
            pass
        
        sJson = json.dumps(obj,\
                               sort_keys=True, \
                               indent = 4, \
                               ensure_ascii=True, \
                               cls=CaseClassEncoder)
        return sJson

    def pyflowline_print(self):
        """
        Print the flowline case object
        """
        print(self.tojson())
        return

    def pyflowline_export_config_to_json(self, sFilename_output_in = None):
        """
        Export the configuration to a json file

        Args:
            sFilename_output_in (str, optional): The json filename. Defaults to None.
        """

        if self.iFlag_standalone == 1:
            if sFilename_output_in is not None:
                sFilename_output = sFilename_output_in
            else:
                #use current output path
                sFilename_output = os.path.join(self.sWorkspace_output, 'configuration.json' )

            #all basins
            sName = 'configuration_basin.json'
            sFilename_configuration  =  os.path.join( self.sWorkspace_output  , sName)
            with open(sFilename_configuration, 'w', encoding='utf-8') as f:
                sJson = json.dumps([json.loads(ob.tojson()) for ob in self.aBasin],\
                                   sort_keys=True, \
                                   indent = 4)
                f.write(sJson)
                f.close()
                #update
            self.sFilename_basins = sFilename_configuration

        else:
            if sFilename_output_in is not None:
                sFilename_output = sFilename_output_in
            else:
                #use parent path
                sPath = Path(self.sWorkspace_output)
                sFilename_output = os.path.join(sPath.parent.absolute(), 'configuration.json' )
                #all basins
            sPath = Path(self.sWorkspace_output)
            sName = 'configuration_basin.json'
            sFilename_configuration  =  os.path.join( sPath.parent.absolute() , sName)
            with open(sFilename_configuration, 'w', encoding='utf-8') as f:
                sJson = json.dumps([json.loads(ob.tojson()) for ob in self.aBasin],\
                                   sort_keys=True, \
                                   indent = 4)
                f.write(sJson)
                f.close()
                #update for pyhexwatershed
            self.sFilename_basins = sFilename_configuration
            self.sWorkspace_output =   sPath.parent.absolute()


        aSkip = ['aBasin', \
                 'aFlowline_simplified','aFlowline_conceptual','aCellID_outlet',
                 'aCell']

        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)

        with open(sFilename_output, 'w', encoding='utf-8') as f:
            json.dump(obj, f,sort_keys=True, \
                      ensure_ascii=False, \
                      indent=4, \
                      cls=CaseClassEncoder)
        return

    def pyflowline_export_basin_config_to_json(self, sFilename_output_in= None):
        """
        Export the member basin configuration to a json file

        Args:
            sFilename_output_in (str, optional): The json filename. Defaults to None.
        """
        if self.iFlag_standalone == 1:
            if sFilename_output_in is not None:
                sFilename_output = sFilename_output_in
            else:
                #use current output path
                sName = 'configuration_basin.json'
                sFilename_output  =  os.path.join( self.sWorkspace_output, sName)

            #all basins
            with open(sFilename_output, 'w', encoding='utf-8') as f:
                sJson = json.dumps([json.loads(ob.tojson()) for ob in self.aBasin],\
                                   sort_keys=True, \
                                   indent = 4)
                f.write(sJson)
                f.close()
                #update
            self.sFilename_basins = sFilename_output

        else:
            if sFilename_output_in is not None:
                sFilename_output = sFilename_output_in
            else:
                #use current output path
                sPath = Path(self.sWorkspace_output)
                sName = 'configuration_basin.json'
                sFilename_output  =  os.path.join( sPath.parent.absolute() , sName)

            #all basins
            with open(sFilename_output, 'w', encoding='utf-8') as f:
                sJson = json.dumps([json.loads(ob.tojson()) for ob in self.aBasin],\
                                   sort_keys=True, \
                                   indent = 4)
                f.write(sJson)
                f.close()
                #update for pyhexwatershed
            self.sFilename_basins = sFilename_output

        return
