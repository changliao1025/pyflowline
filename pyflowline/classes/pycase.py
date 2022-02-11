import os
from pathlib import Path
from abc import ABCMeta
import json
from json import JSONEncoder
import datetime

import numpy as np
import utm
import requests
from osgeo import ogr, osr, gdal

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
import matplotlib.cm as cm

from shapely.wkt import loads
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.io.shapereader as shpreader
from cartopy.feature import ShapelyFeature
import cartopy.crs as ccrs

from pyflowline.classes.mpas import pympas
from pyflowline.classes.hexagon import pyhexagon
from pyflowline.classes.latlon import pylatlon
from pyflowline.classes.square import pysquare
from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.basin import pybasin
from pyflowline.classes.flowline import pyflowline
from pyflowline.classes.edge import pyedge

from pyflowline.formats.convert_shapefile_to_json import convert_shapefile_to_json_swat
from pyflowline.formats.read_mesh import read_mesh_json

from pyflowline.algorithms.auxiliary.text_reader_string import text_reader_string
from pyflowline.algorithms.auxiliary.gdal_functions import reproject_coordinates, reproject_coordinates_batch
from pyflowline.algorithms.auxiliary.gdal_functions  import degree_to_meter
from pyflowline.algorithms.auxiliary.gdal_functions  import meter_to_degree
from pyflowline.algorithms.auxiliary.gdal_functions import gdal_read_geotiff_file, retrieve_geotiff_metadata
from pyflowline.algorithms.auxiliary.gdal_functions import  Google_MetersPerPixel

from pyflowline.mesh.hexagon.create_hexagon_mesh import create_hexagon_mesh
from pyflowline.mesh.latlon.create_latlon_mesh import create_latlon_mesh
from pyflowline.mesh.square.create_square_mesh import create_square_mesh
from pyflowline.mesh.mpas.create_mpas_mesh import create_mpas_mesh
from pyflowline.mesh.tin.create_tin_mesh import create_tin_mesh

#desired_proj = ccrs.Orthographic(central_longitude=-75, central_latitude=42, globe=None)
desired_proj = ccrs.PlateCarree()

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
            return json.loads(obj.tojson()) #lVertexID
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
        if isinstance(obj, pybasin):
            return obj.lBasinID    
            
        return JSONEncoder.default(self, obj)

class flowlinecase(object):
    __metaclass__ = ABCMeta
    iCase_index= 0
    sMesh_type = 1
    iFlag_standalone=1
    iFlag_global = 0
    iFlag_multiple_outlet = 0
    iFlag_use_mesh_dem=0
    iFlag_save_mesh = 0
    iFlag_simplification = 1 #user can turn on/off
    iFlag_create_mesh=1
    iFlag_intersect = 1
    
    iFlag_rotation=0

    nOutlet = 1 #by default , there shoule ne only one ouelet
    dResolution_degree=0.0
    dResolution_meter=0.0
    dThreshold_small_river=0.0
    dLongitude_left = -180
    dLongitude_right = 180
    dLatitude_bot = -90
    dLatitude_top = 90
    sFilename_model_configuration=''
    sWorkspace_data=''     
    sWorkspace_project=''    
    sWorkspace_output=''    
    sRegion=''
    sModel=''
    sMesh_type ='hexagon'

    sCase=''
    sDate=''    

    sFilename_spatial_reference=''
    sFilename_dem=''   

    sFilename_mesh=''
    sFilename_mesh_info=''
    sFilename_mesh_netcdf=''
    

    aBasin = list()
    aFlowline_simplified=list()
    aFlowline_conceptual=list()
    aCellID_outlet = list()
    aCell=list()

    
    def __init__(self, aParameter,\
            iFlag_standalone_in= None,\
            sModel_in = None,\
            sWorkspace_output_in = None):
        #flags
        if iFlag_standalone_in is not None:
            self.iFlag_standalone = iFlag_standalone_in
        else:
            if 'iFlag_standalone' in aParameter:
                self.iFlag_standalone = int(aParameter['iFlag_standalone'])
            else:
                self.iFlag_standalone=1

        if 'iFlag_flowline' in aParameter:
            self.iFlag_flowline             = int(aParameter[ 'iFlag_flowline'])
        
        

        if 'iFlag_global' in aParameter:
            self.iFlag_global             = int(aParameter[ 'iFlag_global'])
        
        if 'iFlag_multiple_outlet' in aParameter:
            self.iFlag_multiple_outlet             = int(aParameter[ 'iFlag_multiple_outlet'])  
        
        if 'iFlag_simplification' in aParameter:
            self.iFlag_simplification = int(aParameter['iFlag_simplification'])


        if 'iFlag_create_mesh' in aParameter:
            self.iFlag_create_mesh = int(aParameter['iFlag_create_mesh'])    
        
        if 'iFlag_save_mesh' in aParameter:
            self.iFlag_save_mesh             = int(aParameter[ 'iFlag_save_mesh'])
        
        if 'iFlag_use_mesh_dem' in aParameter:
            self.iFlag_use_mesh_dem = int(aParameter['iFlag_use_mesh_dem'])

        if 'iFlag_rotation' in aParameter:
            self.iFlag_rotation = int(aParameter['iFlag_rotation'])
        

        if 'iFlag_intersect' in aParameter:
            self.iFlag_intersect = int(aParameter['iFlag_intersect'])

        if 'nOutlet' in aParameter:
            self.nOutlet = int(aParameter['nOutlet'])
             
        if 'iCase_index' in aParameter:
            iCase_index = int(aParameter['iCase_index'])
        else:
            iCase_index = 1
        sCase_index = "{:03d}".format( iCase_index )
        self.iCase_index =   iCase_index

        if 'dResolution_degree' in aParameter:
            self.dResolution_degree = float(aParameter['dResolution_degree']) 

        if 'dResolution_meter' in aParameter:
            self.dResolution_meter = float(aParameter['dResolution_meter']) 
        else:
            print('Please specify resolution.')

        if 'dLongitude_left' in aParameter:
            self.dLongitude_left = float(aParameter['dLongitude_left']) 

        if 'dLongitude_right' in aParameter:
            self.dLongitude_right = float(aParameter['dLongitude_right']) 

        if 'dLatitude_bot' in aParameter:
            self.dLatitude_bot = float(aParameter['dLatitude_bot']) 

        if 'dLatitude_top' in aParameter:
            self.dLatitude_top = float(aParameter['dLatitude_top']) 

        if 'sFilename_model_configuration' in aParameter:
            self.sFilename_model_configuration    = aParameter[ 'sFilename_model_configuration']

        if 'sFilename_spatial_reference' in aParameter:
            self.sFilename_spatial_reference = aParameter['sFilename_spatial_reference']

        if 'sFilename_dem' in aParameter:
            self.sFilename_dem = aParameter['sFilename_dem']

        if 'sFilename_mesh_netcdf' in aParameter:
            self.sFilename_mesh_netcdf = aParameter['sFilename_mesh_netcdf']
        
        if 'sWorkspace_bin' in aParameter:
            self.sWorkspace_bin= aParameter[ 'sWorkspace_bin']
            
        if 'sWorkspace_data' in aParameter:
            self.sWorkspace_data = aParameter[ 'sWorkspace_data']
        
        if 'sWorkspace_project' in aParameter:
            self.sWorkspace_project= aParameter[ 'sWorkspace_project']
        
        if sWorkspace_output_in is not None:
            self.sWorkspace_output = sWorkspace_output_in
        else:
            if 'sWorkspace_output' in aParameter:
                self.sWorkspace_output = aParameter[ 'sWorkspace_output']
        
        if 'sJob' in aParameter:
            self.sJob =  aParameter['sJob'] 
        if 'sRegion' in aParameter:
            self.sRegion               = aParameter[ 'sRegion']
        
        if sModel_in is not None:
            self.sModel = sModel_in
        else:
            if 'sModel' in aParameter:
                self.sModel                = aParameter[ 'sModel']

                      
        sDate   = aParameter[ 'sDate']
        if sDate is not None:
            self.sDate= sDate
        else:
            self.sDate = sDate_default
        
        
        sCase = self.sModel  + self.sDate + sCase_index
        self.sCase = sCase

        if 'sMesh_type' in aParameter:
            self.sMesh_type =  aParameter['sMesh_type']
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
                        if sMesh_type =='tin': #tin
                            self.iMesh_type = 5
                        else:
                            print('Unsupported mesh type?')
        
        #the model can be run as part of hexwatershed or standalone
        if self.iFlag_standalone ==1:
            sPath = str(Path(self.sWorkspace_output)  /  sCase)
            self.sWorkspace_output = sPath
        else:
            sPath = self.sWorkspace_output
        
        Path(sPath).mkdir(parents=True, exist_ok=True)

        self.aBasin = list()
        if self.iFlag_flowline == 1:
            if 'sFilename_basins' in aParameter:
                self.sFilename_basins = aParameter['sFilename_basins']
                with open(self.sFilename_basins) as json_file:
                    dummy_data = json.load(json_file)     
                    for i in range(self.nOutlet):
                        sBasin =  "{:03d}".format(i+1)   
                        dummy_basin = dummy_data[i]
                        dummy_basin['sWorkspace_output_basin'] = str(Path(self.sWorkspace_output) / sBasin )
                        
                        pBasin = pybasin(dummy_basin)
    
                        self.aBasin.append(pBasin)
            else:
                pass
     

        #model generated files   
        self.sFilename_mesh = os.path.join(str(Path(self.sWorkspace_output)  ) , sMesh_type + ".json" )               
        self.sFilename_mesh_info= os.path.join(str(Path(self.sWorkspace_output)  ) , sMesh_type + "_mesh_info.json"  )    
        self.sWorkspace_data_project = str(Path(self.sWorkspace_data ) / self.sWorkspace_project)
                
        return
        

    def convert_flowline_to_json(self):
        for pBasin in self.aBasin:            
            pBasin.convert_flowline_to_json()
            pass
 
    def flowline_simplification(self):
        aFlowline_out = list()   #store all the flowline
        if self.iFlag_simplification == 1: 
            for pBasin in self.aBasin:
                aFlowline_basin = pBasin.flowline_simplification()                
                aFlowline_out = aFlowline_out + aFlowline_basin

            self.aFlowline_simplified = aFlowline_out

        #export the outlet into a since file
        aOutlet = list()
        if self.iFlag_simplification == 1: 
            for pBasin in self.aBasin:
                aOutlet.append(pBasin.pVertex_outlet)

        return aFlowline_out
    
    def mesh_generation(self):        
        iFlag_global =  self.iFlag_global
        iMesh_type = self.iMesh_type
        iFlag_save_mesh = self.iFlag_save_mesh
        iFlag_rotation = self.iFlag_rotation
        dResolution_degree = self.dResolution_degree
        dResolution_meter = self.dResolution_meter
        sFilename_dem = self.sFilename_dem
        sFilename_spatial_reference = self.sFilename_spatial_reference
        sFilename_mesh = self.sFilename_mesh
        if iMesh_type !=4: #hexagon
            dPixelWidth, dOriginX, dOriginY, nrow, ncolumn, pSpatialRef_dem, pProjection, pGeotransform\
                 = retrieve_geotiff_metadata(sFilename_dem)
            spatial_reference_source = pSpatialRef_dem
            spatial_reference_target = osr.SpatialReference()  
            spatial_reference_target.ImportFromEPSG(4326)

            dY_bot = dOriginY - (nrow+1) * dPixelWidth
            dLongitude_left,  dLatitude_bot= reproject_coordinates(dOriginX, dY_bot,pSpatialRef_dem,spatial_reference_target)
            dX_right = dOriginX + (ncolumn +1) * dPixelWidth

            dLongitude_right, dLatitude_top= reproject_coordinates(dX_right, dOriginY,pSpatialRef_dem,spatial_reference_target)
            dLatitude_mean = 0.5 * (dLatitude_top + dLatitude_bot)


            if dResolution_meter < 0:
                #not used
                pass
            else:
                dResolution_degree = meter_to_degree(dResolution_meter, dLatitude_mean)


            dX_left = dOriginX
            dY_top = dOriginY
        else:
            pass

        if iMesh_type ==1: #hexagon

            #hexagon edge
            dResolution_meter = degree_to_meter(dLatitude_mean, dResolution_degree )
            dArea = np.power(dResolution_meter,2.0)
            dLength_edge = np.sqrt(  2.0 * dArea / (3.0* np.sqrt(3.0))  )
            if iFlag_rotation ==0:            
                dX_spacing = dLength_edge * np.sqrt(3.0)
                dY_spacing = dLength_edge * 1.5
                ncolumn= int( (dX_right - dX_left) / dX_spacing )
                nrow= int( (dY_top - dY_bot) / dY_spacing ) 
            else:            
                dX_spacing = dLength_edge * 1.5
                dY_spacing = dLength_edge * np.sqrt(3.0)    
                ncolumn= int( (dX_right - dX_left) / dX_spacing )+1
                nrow= int( (dY_top - dY_bot) / dY_spacing )

            aHexagon = create_hexagon_mesh(iFlag_rotation, dX_left, dY_bot, dResolution_meter, ncolumn, nrow, \
                sFilename_mesh, sFilename_spatial_reference)
            return aHexagon
        else:
            if iMesh_type ==2: #sqaure
                ncolumn= int( (dX_right - dX_left) / dResolution_meter )
                nrow= int( (dY_top - dY_bot) / dResolution_meter )

                aSquare = create_square_mesh(dX_left, dY_bot, dResolution_meter, ncolumn, nrow, \
                    sFilename_mesh, sFilename_spatial_reference)
                return aSquare
            else:
                if iMesh_type ==3: #latlon
                    dResolution_meter = degree_to_meter(dLatitude_mean, dResolution_degree)
                    dArea = np.power(dResolution_meter,2.0)
                    dLatitude_top    = self.dLatitude_top   
                    dLatitude_bot    = self.dLatitude_bot   
                    dLongitude_left  = self.dLongitude_left 
                    dLongitude_right = self.dLongitude_right
                    ncolumn= int( (dLongitude_right - dLongitude_left) / dResolution_degree )
                    nrow= int( (dLatitude_top - dLatitude_bot) / dResolution_degree )
                    aLatlon = create_latlon_mesh(dLongitude_left, dLatitude_bot, dResolution_degree, ncolumn, nrow, \
                        sFilename_mesh)
                    return aLatlon
                else:
                    if iMesh_type == 4: #mpas
                        iFlag_use_mesh_dem = self.iFlag_use_mesh_dem
                        sFilename_mesh_netcdf = self.sFilename_mesh_netcdf
                        dLatitude_top    = self.dLatitude_top   
                        dLatitude_bot    = self.dLatitude_bot   
                        dLongitude_left  = self.dLongitude_left 
                        dLongitude_right = self.dLongitude_right
                        aMpas = create_mpas_mesh(iFlag_global, iFlag_use_mesh_dem, iFlag_save_mesh, \
                              dLongitude_left, dLongitude_right,  dLatitude_top, dLatitude_bot, \
                                    sFilename_mesh_netcdf,      sFilename_mesh)
                        return aMpas
                    else:
                        if iMesh_type ==5: #tin this one need to be updated because central location issue
                            #tin edge
                            dArea = np.power(dResolution_meter,2.0)
                            dLength_edge = np.sqrt(  4.0 * dArea /  np.sqrt(3.0) )  
                            dX_shift = 0.5 * dLength_edge
                            dY_shift = 0.5 * dLength_edge * np.sqrt(3.0) 
                            dX_spacing = dX_shift * 2
                            dY_spacing = dY_shift
                            ncolumn= int( (dX_right - dX_left) / dX_shift )
                            nrow= int( (dY_top - dY_bot) / dY_spacing ) 
                            aTin = create_tin_mesh(dX_left, dY_bot, dResolution_meter, ncolumn, nrow,sFilename_mesh, sFilename_spatial_reference)
                            return aTin
                        else:
                            print('Unsupported mesh type?')
                            return
        return
    
    def reconstruct_topological_relationship(self):
        iMesh_type = self.iMesh_type
        iFlag_intersect = self.iFlag_intersect
        sWorkspace_output = self.sWorkspace_output
        nOutlet = self.nOutlet
        sFilename_mesh=self.sFilename_mesh
        self.aCell, pSpatial_reference_mesh = read_mesh_json(iMesh_type, sFilename_mesh)
        
        aFlowline_conceptual = list()   #store all the flowline
        aCellID_outlet = list()
        aBasin = list()
        aCell_intersect=list()
        if iFlag_intersect == 1:
            for pBasin in self.aBasin:
                aCell_intersect_basin = pBasin.reconstruct_topological_relationship(iMesh_type,sFilename_mesh)
                aFlowline_conceptual = aFlowline_conceptual + pBasin.aFlowline_basin
                aBasin.append(pBasin)
                aCellID_outlet.append(pBasin.lCellID_outlet)
                aCell_intersect = aCell_intersect + aCell_intersect_basin

                #update length?
            for pCell in self.aCell:
                for pCell2 in aCell_intersect:
                    if pCell2.lCellID == pCell.lCellID:
                        pCell.dLength_flowline = pCell2.dLength_flowline

            
            #save basin json info for hexwatershed model
            sPath = os.path.dirname(self.sFilename_basins)
            sName = str(Path(self.sFilename_basins).stem ) + '_new.json'
            sFilename_configuration  =  os.path.join( sPath  , sName)

            with open(sFilename_configuration, 'w', encoding='utf-8') as f:
                sJson = json.dumps([json.loads(ob.tojson()) for ob in aBasin],\
                    sort_keys=True, \
                    indent = 4)        

                f.write(sJson)    
                f.close()
            
            self.aFlowline_conceptual = aFlowline_conceptual
            self.aCellID_outlet = aCellID_outlet

            return   aFlowline_conceptual, aCellID_outlet

        return

    def analyze(self):
        for pBasin in self.aBasin:
            pBasin.analyze()
        return
        
    def evaluate(self):
        for pBasin in self.aBasin:
            pBasin.evaluate(self.iMesh_type, self.sMesh_type)
        return

    def setup(self):
        self.convert_flowline_to_json()
        return

    def run(self):
        #self.flowline_simplification()
        #self.mesh_generation()
        self.reconstruct_topological_relationship()
        self.export()
        return

    def plot(self, sVariable_in=None, aExtent_in = None):
        if sVariable_in == 'mesh':
            self.plot_mesh()
        else:            
            if sVariable_in == 'overlap':
                self.plot_mesh_with_flowline( aExtent_in= aExtent_in)
            else:
                for pBasin in self.aBasin:            
                    pBasin.plot(sVariable_in= sVariable_in)
                    pass
        
        return
    
    def plot_study_area(self):
        sWorkspace_output_case = self.sWorkspace_output

        fig = plt.figure(dpi=300)
        fig.set_figwidth( 4 )
        fig.set_figheight( 4 )

        #plot dem first        
        sFilename_dem= self.sFilename_dem
         
        dummy = gdal_read_geotiff_file(sFilename_dem)
        dResolution_x = dummy[1]
        pProjection = dummy[8]
        pSpatial_reference_source = dummy[9]
        dem_proj =ccrs.AlbersEqualArea(central_longitude=pSpatial_reference_source.GetProjParm('longitude_of_center'), \
            central_latitude=pSpatial_reference_source.GetProjParm('latitude_of_center'),\
                standard_parallels=(pSpatial_reference_source.GetProjParm('standard_parallel_1'),pSpatial_reference_source.GetProjParm('standard_parallel_2')))
        ax_dem = fig.add_axes([0.3, 0.15, 0.6, 0.55] , projection=dem_proj )
        ax_dem.set_xmargin(0.05)
        ax_dem.set_ymargin(0.10)          
        
        missing_value = dummy[6]

        
        dOriginX=dummy[2]
        dOriginY=dummy[3]
        nrow=dummy[4]
        ncolumn=dummy[5]
        dLon_min = dOriginX
        dLon_max = dOriginX + ncolumn * dResolution_x
        dLat_max = dOriginY
        dLat_min = dOriginY - nrow * dResolution_x
       
        pSpatial_reference_target = osr.SpatialReference()  
        pSpatial_reference_target.ImportFromEPSG(4326)
        aLon= list()
        aLat=list()
        aLon.append(dLon_min)
        aLon.append(dLon_max)
        aLat.append(dLat_min)
        aLat.append(dLat_max)
        aLon, aLat = reproject_coordinates_batch(aLon, aLat, pSpatial_reference_source, pSpatial_reference_target)

        dLongitude_center = np.mean(aLon)
        dLatitude_center = np.mean(aLat)
        aImage_extent =  [dLon_min- dResolution_x ,dLon_max + dResolution_x, dLat_min -dResolution_x,  dLat_max+dResolution_x]
        aImage_in = dummy[0]
        aImage_in[np.where(aImage_in == missing_value)] = np.nan

        # 
        demplot = ax_dem.imshow(aImage_in, origin='upper', extent=aImage_extent,cmap=cm.terrain,transform=dem_proj) 
        # change all spines
        for axis in ['top','bottom','left','right']:
            ax_dem.spines[axis].set_linewidth(0.5)

        # increase tick width
        ax_dem.tick_params(width=0.5)
        u = utm.from_latlon(dLatitude_center, dLongitude_center)
        gl = ax_dem.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.3, linestyle='--')
        gl.left_labels=False
        gl.top_labels=False        
       
         # declare text size
        XTEXT_SIZE = 8
        YTEXT_SIZE = 8
        # to facilitate text rotation at bottom edge, ...
        # text justification: 'ha':'right' is used to avoid clashing with map's boundary
        # default of 'ha' is center, often causes trouble when text rotation is not zero
        gl.xlabel_style = {'size': XTEXT_SIZE, 'color': 'k', 'rotation':0, 'ha':'right'}
        gl.ylabel_style = {'size': YTEXT_SIZE, 'color': 'k', 'rotation':90,'weight': 'normal'}

        sFilename_nhd = '/qfs/people/liao313/data/hexwatershed/susquehanna/vector/hydrology/nhd_proj.shp'
        reader = shpreader.Reader(sFilename_nhd)
        shape_feature = ShapelyFeature(reader.geometries(),
                                dem_proj, facecolor='none')
        ax_dem.add_feature(shape_feature)
        pSpatial_reference_source = osr.SpatialReference()  
        pSpatial_reference_source.ImportFromEPSG(4326)
        pSpatial_reference_target = osr.SpatialReference()  
        pSpatial_reference_target.ImportFromWkt(dem_proj.to_wkt())
        for pBasin in self.aBasin:
            sFilename_dam = pBasin.sFilename_dam

            aData_dam = text_reader_string(sFilename_dam, iSkipline_in =1,cDelimiter_in=',' )
            ndam = len(aData_dam)
            for idam in range(ndam):
                dLon = float(aData_dam[idam, 1])
                dLat = float(aData_dam[idam, 0])
                x,y = reproject_coordinates(dLon, dLat,pSpatial_reference_source, pSpatial_reference_target )
                ax_dem.plot(x,y,marker = 'x', color='red', markersize=4)

                pass

        ax_cb= fig.add_axes([0.2, 0.2, 0.02, 0.5])
        
        cb = plt.colorbar(demplot, cax = ax_cb, extend = 'both')
        cb.ax.get_yaxis().set_ticks_position('left')
        cb.ax.get_yaxis().labelpad = 10
        cb.ax.set_ylabel('Unit: meter', rotation=270)
        cb.ax.tick_params(labelsize=6) 

        #google earth
        zone = u[2]
        ge_proj = ccrs.Mercator()#central_longitude=dLongitude_center)
        
        pSpatial_reference_source = osr.SpatialReference()  
        pSpatial_reference_source.ImportFromEPSG(4326)

        pSpatial_reference_target = osr.SpatialReference()  
        

        iZoom = 7       
        pSpatial_reference_target.ImportFromWkt(ge_proj.to_wkt())
        uv_xcenter, uv_ycenter = reproject_coordinates(dLongitude_center, dLatitude_center, pSpatial_reference_source, pSpatial_reference_target)
        xsize_ge = 600
        ysize_ge = 600
        scale = Google_MetersPerPixel(iZoom)
        xrange = (uv_xcenter - (xsize_ge/2.0*scale), uv_xcenter + (xsize_ge/2.0*scale))
        yrange = (uv_ycenter - (ysize_ge/2.0*scale), uv_ycenter + (ysize_ge/2.0*scale))


        aImage_extent = [xrange[0], xrange[1],yrange[0], yrange[1] ]
        ax_ge = fig.add_axes([0.1, 0.75, 0.2, 0.2] , projection=ge_proj )
         
        sLongitude_center ="{:0f}".format(dLongitude_center)
        sLatitude_center= "{:0f}".format(dLatitude_center)
        sZoom="{:0d}".format(iZoom)
        
        sResolution = "{:0d}".format(xsize_ge) + 'x' + "{:0d}".format(ysize_ge)
        sMap_type="hybrid"
        sGoogleMap = "http://maps.googleapis.com/maps/api/staticmap?" + \
            "center=" + sLatitude_center + ',' + sLongitude_center + \
            "&zoom=" + sZoom + "&size=" + sResolution + \
            "&maptype="+sMap_type+"&sensor=false&format=png32" + "&key=AIzaSyCT0NPGHJqRFysOYELOWBdqov6AbphgaFY"

        r = requests.get(sGoogleMap)  
        # wb mode is stand for write binary mode
        f = open('google_map.png', 'wb')
        # r.content gives content,
        # in this case gives image
        f.write(r.content)
        # close method of file object
        # save and close the file
        f.close()
        

        img = mpimg.imread('google_map.png')
        ax_ge.imshow(img,extent=aImage_extent,       transform=ge_proj)
        gl_ge = ax_ge.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.3, linestyle='--')
        #gl_ge.xlocator = mticker.FixedLocator(aLon_range)
        #gl_ge.ylocator = mticker.FixedLocator(aLat_range)
        gl_ge.xformatter = LongitudeFormatter()
        gl_ge.yformatter = LatitudeFormatter()
        gl_ge.right_labels=False
        gl_ge.bottom_labels=False        
       
        # declare text size
        XTEXT_SIZE = 2.5
        YTEXT_SIZE = 2.5
        # to facilitate text rotation at bottom edge, ...
        # text justification: 'ha':'right' is used to avoid clashing with map's boundary
        # default of 'ha' is center, often causes trouble when text rotation is not zero
        gl_ge.xlabel_style = {'size': XTEXT_SIZE, 'color': 'k', 'rotation':0, 'ha':'right'}
        gl_ge.ylabel_style = {'size': YTEXT_SIZE, 'color': 'k', 'rotation':90,'weight': 'normal'}
        #add boundary
        sFilename_boundary = '/qfs/people/liao313/data/hexwatershed/susquehanna/vector/hydrology/boundary_ge.shp'
        reader = shpreader.Reader(sFilename_boundary)
        shape_feature = ShapelyFeature(reader.geometries(),
                                ge_proj, facecolor='none')
        ax_ge.add_feature(shape_feature)

        #histogram

        ax_histo = fig.add_axes([0.4, 0.8, 0.4, 0.15]  )
        sFilename_slope = '/qfs/people/liao313/data/hexwatershed/susquehanna/raster/dem/slope_ext.tif'
        dummy = gdal_read_geotiff_file(sFilename_slope)
        aSlope= dummy[0]
        missing_value = dummy[6]
        aSlope = aSlope[np.where(aSlope != missing_value)]
        dMax_x=60
        dMin_x=0
        dSpace_x=4
        ax_histo.hist(aSlope,  int(  (dMax_x-dMin_x)/ dSpace_x), color = "skyblue", ec="skyblue")  
        sLabel_x= 'Slope (degree)'
        sLabel_y ='Frequency'
        ax_histo.set_xlabel(sLabel_x,fontsize=4 )
        ax_histo.set_ylabel(sLabel_y,fontsize=4 )       
        ax_histo.set_xlim( dMin_x, dMax_x )

        formatter = mticker.ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)        
        ax_histo.yaxis.set_major_formatter(formatter)
        ax_histo.tick_params(axis='x', labelsize=5 )
        ax_histo.tick_params(axis='y', labelsize=5 )
        ax_histo.yaxis.offsetText.set_fontsize(5)

        sText =  'Outlet'
        ax_dem.text(0.85, 0.05, sText, \
                verticalalignment='top', horizontalalignment='left',\
                transform=ax_dem.transAxes, \
                color='black', fontsize=7)

        sFilename  = 'study_area.png'
        sFilename_out = os.path.join(sWorkspace_output_case, sFilename)
        plt.savefig(sFilename_out, bbox_inches='tight')

        
        
        return

    def plot_mesh(self, aExtent_in=None):

        sWorkspace_output_case = self.sWorkspace_output

        sFilename_json  =  self.sFilename_mesh
        sMesh_type = self.sMesh_type

        fig = plt.figure(dpi=300)
        fig.set_figwidth( 4 )
        fig.set_figheight( 4 )
        ax = fig.add_axes([0.1, 0.15, 0.75, 0.7] , projection=desired_proj )

        ax.set_global()
        pDriver = ogr.GetDriverByName('GeoJSON')
        pDataset = pDriver.Open(sFilename_json, gdal.GA_ReadOnly)
        pLayer = pDataset.GetLayer(0)
    
        pSrs = osr.SpatialReference()  
        pSrs.ImportFromEPSG(4326)    # WGS84 lat/lon

        dLat_min = 90
        dLat_max = -90
        dLon_min = 180
        dLon_max = -180
        
        for pFeature in pLayer:
            pGeometry_in = pFeature.GetGeometryRef()
            sGeometry_type = pGeometry_in.GetGeometryName()


            lID =0 
            if sGeometry_type =='POLYGON':
                dummy0 = loads( pGeometry_in.ExportToWkt() )
                aCoords_gcs = dummy0.exterior.coords
                aCoords_gcs= np.array(aCoords_gcs)
                nvertex = len(aCoords_gcs)

                for i in range(nvertex):
                    dLon = aCoords_gcs[i][0]
                    dLat = aCoords_gcs[i][1]
                    if dLon > dLon_max:
                        dLon_max = dLon

                    if dLon < dLon_min:
                        dLon_min = dLon

                    if dLat > dLat_max:
                        dLat_max = dLat

                    if dLat < dLat_min:
                        dLat_min = dLat


                polygon = mpatches.Polygon(aCoords_gcs[:,0:2], closed=True,  linewidth=0.25, \
                    alpha=0.8, edgecolor = 'black',facecolor='none', \
                        transform=ccrs.PlateCarree() )

                ax.add_patch(polygon)                   


        
        if aExtent_in is None:
            marginx  = (dLon_max - dLon_min) / 20
            marginy  = (dLat_max - dLat_min) / 20
            aExtent_in = [dLon_min - marginx , dLon_max + marginx , dLat_min - marginy , dLat_max + marginy]
       
        
        ax.set_extent( aExtent_in )  
    

        #add mesh info
        ax.set_title( sMesh_type.title() +  ' mesh')

        #ax.text(0.95, 0.95, sMesh_type, \
        #verticalalignment='top', horizontalalignment='right',\
        #        transform=ax.transAxes, \
        #        color='black', fontsize=15)

        ax.coastlines()#resolution='110m')
        ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.3, linestyle='--')


        sDirname = os.path.dirname(sFilename_json)
        sFilename  = Path(sFilename_json).stem + '.png'
        sFilename_out = os.path.join(sDirname, sFilename)
        plt.savefig(sFilename_out, bbox_inches='tight')

        pDataset = pLayer = pFeature  = None   
        #plt.show()   
        return

    def plot_mesh_with_flowline(self, aExtent_in=None):
        sWorkspace_output_case = self.sWorkspace_output

        sFilename_mesh  =  self.sFilename_mesh
        sMesh_type = self.sMesh_type

        fig = plt.figure( dpi=300 )
        fig.set_figwidth( 4 )
        fig.set_figheight( 4 )
        ax = fig.add_axes([0.1, 0.15, 0.75, 0.7] , projection=desired_proj )

        ax.set_global()
        pDriver = ogr.GetDriverByName('GeoJSON')
        pDataset = pDriver.Open(sFilename_mesh, gdal.GA_ReadOnly)
        pLayer = pDataset.GetLayer(0)
    
        pSrs = osr.SpatialReference()  
        pSrs.ImportFromEPSG(4326)    # WGS84 lat/lon

        dLat_min = 90
        dLat_max = -90
        dLon_min = 180
        dLon_max = -180
        
        for pFeature in pLayer:
            pGeometry_in = pFeature.GetGeometryRef()
            sGeometry_type = pGeometry_in.GetGeometryName()            
            if sGeometry_type =='POLYGON':
                dummy0 = loads( pGeometry_in.ExportToWkt() )
                aCoords_gcs = dummy0.exterior.coords
                aCoords_gcs= np.array(aCoords_gcs)
                nvertex = len(aCoords_gcs)

                for i in range(nvertex):
                    dLon = aCoords_gcs[i][0]
                    dLat = aCoords_gcs[i][1]
                    if dLon > dLon_max:
                        dLon_max = dLon

                    if dLon < dLon_min:
                        dLon_min = dLon

                    if dLat > dLat_max:
                        dLat_max = dLat

                    if dLat < dLat_min:
                        dLat_min = dLat

                polygon = mpatches.Polygon(aCoords_gcs[:,0:2], closed=True,   linewidth=0.25, \
                    alpha=0.8, edgecolor = 'black',facecolor='none', \
                        transform=ccrs.PlateCarree() )

                ax.add_patch(polygon)                   

        #draw base flowline first with black color
        lID = 0 
        for pBasin in self.aBasin:
            sWorkspace_output_basin=  pBasin.sWorkspace_output_basin                        
            sFilename = pBasin.sFilename_flowline_segment_order_before_intersect
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename)
            sFilename_json = os.path.join(sWorkspace_output_basin, sFilename_out)
            pDriver = ogr.GetDriverByName('GeoJSON')
            pDataset = pDriver.Open(sFilename_json, gdal.GA_ReadOnly)
            pLayer = pDataset.GetLayer(0)
            n_colors = pLayer.GetFeatureCount()     
            for pFeature in pLayer:
                pGeometry_in = pFeature.GetGeometryRef()
                sGeometry_type = pGeometry_in.GetGeometryName()
                if sGeometry_type =='LINESTRING':
                    dummy0 = loads( pGeometry_in.ExportToWkt() )
                    aCoords_gcs = dummy0.coords
                    aCoords_gcs= np.array(aCoords_gcs)
                    nvertex = len(aCoords_gcs)   
                    codes = np.full(nvertex, mpath.Path.LINETO, dtype=int )
                    codes[0] = mpath.Path.MOVETO
                    path = mpath.Path(aCoords_gcs, codes)            
                    x, y = zip(*path.vertices)
                    line, = ax.plot(x, y, color= 'black', linewidth=0.5)
                    lID = lID + 1
                pass
            pass

        #plot new flowline now
        lID = 0 
        for pBasin in self.aBasin:
            sWorkspace_output_basin=  pBasin.sWorkspace_output_basin            
            sFilename_out = pBasin.sFilename_flowline_final           
            sFilename_json = os.path.join(sWorkspace_output_basin, sFilename_out)
            pDriver = ogr.GetDriverByName('GeoJSON')
            pDataset = pDriver.Open(sFilename_json, gdal.GA_ReadOnly)
            pLayer = pDataset.GetLayer(0)
            n_colors = pLayer.GetFeatureCount()
        
            colours = cm.rainbow(np.linspace(0, 1, n_colors))
            for pFeature in pLayer:
                pGeometry_in = pFeature.GetGeometryRef()
                sGeometry_type = pGeometry_in.GetGeometryName()
                if sGeometry_type =='LINESTRING':
                    dummy0 = loads( pGeometry_in.ExportToWkt() )
                    aCoords_gcs = dummy0.coords
                    aCoords_gcs= np.array(aCoords_gcs)
                    nvertex = len(aCoords_gcs)   

                    codes = np.full(nvertex, mpath.Path.LINETO, dtype=int )
                    codes[0] = mpath.Path.MOVETO
                    path = mpath.Path(aCoords_gcs, codes)            
                    x, y = zip(*path.vertices)
                    line, = ax.plot(x, y, color= colours[lID], linewidth=1)
                    lID = lID + 1
                pass
            pass
        
        sDirname = os.path.dirname(sFilename_mesh)
        if aExtent_in is None:
            marginx  = (dLon_max - dLon_min) / 20
            marginy  = (dLat_max - dLat_min) / 20
            aExtent_in = [dLon_min - marginx , dLon_max + marginx , dLat_min - marginy , dLat_max + marginy]           
                
        
        sTitle = 'Conceptual flowline'
        sFilename  = Path(sFilename_mesh).stem + '_flowline_final.png'                      
        
        ax.set_extent( aExtent_in )        

        ax.coastlines()#resolution='110m')
        ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.2, color='gray', alpha=0.3, linestyle='--')
        
        ax.set_title( sTitle )
        sText = 'Mesh type: ' + sMesh_type.title()
        ax.text(0.05, 0.95, sText, \
        verticalalignment='top', horizontalalignment='left',\
                transform=ax.transAxes, \
                color='black', fontsize=8)
        
        sResolution =  'Resolution: ' + "{:0d}".format( int(self.dResolution_meter) ) + 'm'

        if self.sMesh_type != 'mpas':
            ax.text(0.05, 0.90, sResolution, \
                verticalalignment='top', horizontalalignment='left',\
                transform=ax.transAxes, \
                color='black', fontsize=8)

        
        
        sFilename_out = os.path.join(sDirname, sFilename)
        plt.savefig(sFilename_out, bbox_inches='tight')

        pDataset = pLayer = pFeature  = None   
        #plt.show()   
        return

    def compare_with_raster_dem_method(self, sFilename_dem_flowline, aExtent_in=None):
        sWorkspace_output_case = self.sWorkspace_output

        sFilename_mesh  =  self.sFilename_mesh
        sMesh_type = self.sMesh_type

        fig = plt.figure( dpi=300 )
        fig.set_figwidth( 4 )
        fig.set_figheight( 4 )
        ax = fig.add_axes([0.1, 0.15, 0.75, 0.7] , projection=desired_proj )

        ax.set_global()
        pDriver = ogr.GetDriverByName('GeoJSON')
        pDataset = pDriver.Open(sFilename_mesh, gdal.GA_ReadOnly)
        pLayer = pDataset.GetLayer(0)
    
        pSrs = osr.SpatialReference()  
        pSrs.ImportFromEPSG(4326)    # WGS84 lat/lon

        dLat_min = 90
        dLat_max = -90
        dLon_min = 180
        dLon_max = -180
        
        for pFeature in pLayer:
            pGeometry_in = pFeature.GetGeometryRef()
            sGeometry_type = pGeometry_in.GetGeometryName()            
            if sGeometry_type =='POLYGON':
                dummy0 = loads( pGeometry_in.ExportToWkt() )
                aCoords_gcs = dummy0.exterior.coords
                aCoords_gcs= np.array(aCoords_gcs)
                nvertex = len(aCoords_gcs)

                for i in range(nvertex):
                    dLon = aCoords_gcs[i][0]
                    dLat = aCoords_gcs[i][1]
                    if dLon > dLon_max:
                        dLon_max = dLon

                    if dLon < dLon_min:
                        dLon_min = dLon

                    if dLat > dLat_max:
                        dLat_max = dLat

                    if dLat < dLat_min:
                        dLat_min = dLat

                polygon = mpatches.Polygon(aCoords_gcs[:,0:2], closed=True,   linewidth=0.25, \
                    alpha=0.8, edgecolor = 'black',facecolor='none', \
                        transform=ccrs.PlateCarree() )

                ax.add_patch(polygon)                   

        #draw base flowline first with black color
        lID = 0 
        for pBasin in self.aBasin:
            sWorkspace_output_basin=  pBasin.sWorkspace_output_basin                        
            sFilename = pBasin.sFilename_flowline_segment_order_before_intersect
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename)
            sFilename_json = os.path.join(sWorkspace_output_basin, sFilename_out)
            pDriver = ogr.GetDriverByName('GeoJSON')
            pDataset = pDriver.Open(sFilename_json, gdal.GA_ReadOnly)
            pLayer = pDataset.GetLayer(0)
            n_colors = pLayer.GetFeatureCount()     
            for pFeature in pLayer:
                pGeometry_in = pFeature.GetGeometryRef()
                sGeometry_type = pGeometry_in.GetGeometryName()
                if sGeometry_type =='LINESTRING':
                    dummy0 = loads( pGeometry_in.ExportToWkt() )
                    aCoords_gcs = dummy0.coords
                    aCoords_gcs= np.array(aCoords_gcs)
                    nvertex = len(aCoords_gcs)   
                    codes = np.full(nvertex, mpath.Path.LINETO, dtype=int )
                    codes[0] = mpath.Path.MOVETO
                    path = mpath.Path(aCoords_gcs, codes)            
                    x, y = zip(*path.vertices)
                    line, = ax.plot(x, y, color= 'black', linewidth=0.5)
                    lID = lID + 1
                pass
            pass

        #plot new flowline now
        lID = 0 
        for pBasin in self.aBasin:
            sWorkspace_output_basin=  pBasin.sWorkspace_output_basin            
            sFilename_out = pBasin.sFilename_flowline_final           
            sFilename_json = os.path.join(sWorkspace_output_basin, sFilename_out)
            pDriver = ogr.GetDriverByName('GeoJSON')
            pDataset = pDriver.Open(sFilename_json, gdal.GA_ReadOnly)
            pLayer = pDataset.GetLayer(0)
            
            for pFeature in pLayer:
                pGeometry_in = pFeature.GetGeometryRef()
                sGeometry_type = pGeometry_in.GetGeometryName()
                if sGeometry_type =='LINESTRING':
                    dummy0 = loads( pGeometry_in.ExportToWkt() )
                    aCoords_gcs = dummy0.coords
                    aCoords_gcs= np.array(aCoords_gcs)
                    nvertex = len(aCoords_gcs)   
                    codes = np.full(nvertex, mpath.Path.LINETO, dtype=int )
                    codes[0] = mpath.Path.MOVETO
                    path = mpath.Path(aCoords_gcs, codes)            
                    x, y = zip(*path.vertices)
                    line, = ax.plot(x, y, color= 'blue', linewidth=1)
                    lID = lID + 1
                pass
            pass

        #now swat result
        #conver file first
        for pBasin in self.aBasin:
            sWorkspace_output_basin=  pBasin.sWorkspace_output_basin 
            sFilename_swat = os.path.join(sWorkspace_output_basin, 'swat.json')
            convert_shapefile_to_json_swat(1, sFilename_dem_flowline, sFilename_swat)
            pDriver = ogr.GetDriverByName('GeoJSON')
            pDataset = pDriver.Open(sFilename_swat, gdal.GA_ReadOnly)
            pLayer = pDataset.GetLayer(0)
           
            for pFeature in pLayer:
                pGeometry_in = pFeature.GetGeometryRef()
                sGeometry_type = pGeometry_in.GetGeometryName()
                if sGeometry_type =='LINESTRING':
                    dummy0 = loads( pGeometry_in.ExportToWkt() )
                    aCoords_gcs = dummy0.coords
                    aCoords_gcs= np.array(aCoords_gcs)
                    nvertex = len(aCoords_gcs)   

                    codes = np.full(nvertex, mpath.Path.LINETO, dtype=int )
                    codes[0] = mpath.Path.MOVETO
                    path = mpath.Path(aCoords_gcs, codes)            
                    x, y = zip(*path.vertices)
                    line, = ax.plot(x, y, color= 'red', linewidth=1)
                    lID = lID + 1
                pass
            pass
            #plot new

        
        sDirname = os.path.dirname(sFilename_mesh)
        if aExtent_in is None:
            marginx  = (dLon_max - dLon_min) / 20
            marginy  = (dLat_max - dLat_min) / 20
            aExtent_in = [dLon_min - marginx , dLon_max + marginx , dLat_min - marginy , dLat_max + marginy]           
                
        
        sTitle = 'Conceptual flowline'
        sFilename  = Path(sFilename_mesh).stem + '_flowline_compare_zoom.png'                      
        
        ax.set_extent( aExtent_in )        

        ax.coastlines()#resolution='110m')
        ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.2, color='gray', alpha=0.3, linestyle='--')
        
        ax.set_title( sTitle )
        sText = 'Mesh type: ' + sMesh_type.title()
        ax.text(0.05, 0.95, sText, \
        verticalalignment='top', horizontalalignment='left',\
                transform=ax.transAxes, \
                color='black', fontsize=8)
        
        sResolution =  'Resolution: ' + "{:0d}".format( int(self.dResolution_meter/1000) ) + ' km'

        if self.sMesh_type != 'mpas':
            ax.text(0.05, 0.90, sResolution, \
                verticalalignment='top', horizontalalignment='left',\
                transform=ax.transAxes, \
                color='black', fontsize=8)
        
        sText = 'Topologic relationship-based'
        ax.text(0.05, 0.85, sText, \
        verticalalignment='top', horizontalalignment='left',\
                transform=ax.transAxes, \
                color='blue', fontsize=8)
        
        sText = 'Raster DEM-based'
        ax.text(0.05, 0.80, sText, \
        verticalalignment='top', horizontalalignment='left',\
                transform=ax.transAxes, \
                color='red', fontsize=8)
        
        
        sFilename_out = os.path.join(sDirname, sFilename)
        plt.savefig(sFilename_out, bbox_inches='tight')

        pDataset = pLayer = pFeature  = None   
        return
   
    def export(self):
        
        self.export_mesh_info_to_json()
        for pBasin in self.aBasin:
            pBasin.export()

        self.tojson()

    def export_mesh_info_to_json(self):
        
        aCell_all = self.aCell
        sFilename_json = self.sFilename_mesh_info
        ncell=len(aCell_all)
        iFlag_flowline = self.iFlag_flowline
        aCellID_outlet = self.aCellID_outlet

        aFlowline = self.aFlowline_conceptual

        if iFlag_flowline == 1:
            nFlowline = len(aFlowline)
            for i in range(nFlowline):
                pFlowline = aFlowline[i]
                nEdge = pFlowline.nEdge
                nVertex = pFlowline.nVertex
                aEdge = pFlowline.aEdge

                iStream_segment = pFlowline.iStream_segment
                iStream_order = pFlowline.iStream_order

                for j in range(nEdge):
                    pEdge = aEdge[j]
                    pVertex_start = pEdge.pVertex_start
                    pVertex_end = pEdge.pVertex_end
                    for k in range(ncell):
                        pVertex_center = aCell_all[k].pVertex_center

                        if pVertex_center == pVertex_start:
                            aCell_all[k].iStream_segment_burned = iStream_segment
                            aCell_all[k].iStream_order_burned = iStream_order

                            for l in range(ncell):
                                pVertex_center2 = aCell_all[l].pVertex_center
                                lCellID = aCell_all[l].lCellID
                                if pVertex_center2 == pVertex_end:
                                    aCell_all[k].lCellID_downstream_burned = lCellID
                                    if lCellID in aCellID_outlet:
                                        aCell_all[l].iStream_segment_burned = iStream_segment
                                        aCell_all[l].iStream_order_burned = iStream_order

                                    break
                                
                                

                pass
        else:
            #only mesh, no flowline
            pass

        with open(sFilename_json, 'w', encoding='utf-8') as f:
            sJson = json.dumps([json.loads(ob.tojson()) for ob in aCell_all], indent = 4)        
            f.write(sJson)    
            f.close()

        return

    def tojson(self):
        sJson = json.dumps(self.__dict__, \
            sort_keys=True, \
                indent = 4, \
                    ensure_ascii=True, \
                        cls=CaseClassEncoder)
        return sJson

    def export_config_to_json(self, sFilename_output):
        with open(sFilename_output, 'w', encoding='utf-8') as f:
            json.dump(self.__dict__, f,sort_keys=True, \
                ensure_ascii=False, \
                indent=4, \
                cls=CaseClassEncoder)
        return
