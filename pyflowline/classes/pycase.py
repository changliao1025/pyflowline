import os
from pathlib import Path
from abc import ABCMeta, abstractmethod
import json
from json import JSONEncoder
import numpy as np
import datetime
import json
from osgeo import ogr, osr, gdal, gdalconst
import matplotlib.pyplot as plt
import matplotlib.path as mpath
from shapely.wkt import loads
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from matplotlib.collections import PatchCollection
from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.basin import pybasin
from pyflowline.classes.flowline import pyflowline


from pyflowline.algorithms.auxiliary.gdal_functions import retrieve_geotiff_metadata
from pyflowline.algorithms.auxiliary.gdal_functions import reproject_coordinates
from pyflowline.algorithms.auxiliary.gdal_functions  import degree_to_meter
from pyflowline.algorithms.auxiliary.gdal_functions  import meter_to_degree
from pyflowline.mesh.hexagon.create_hexagon_mesh import create_hexagon_mesh
from pyflowline.mesh.latlon.create_latlon_mesh import create_latlon_mesh
from pyflowline.mesh.square.create_square_mesh import create_square_mesh
from pyflowline.mesh.mpas.create_mpas_mesh import create_mpas_mesh
from pyflowline.mesh.tin.create_tin_mesh import create_tin_mesh



from pyflowline.algorithms.auxiliary.gdal_functions import reproject_coordinates
from pyflowline.formats.read_flowline import read_flowline_shapefile

from pyflowline.formats.read_mesh import read_mesh_json
from pyflowline.formats.read_flowline import read_flowline_geojson
from pyflowline.formats.export_vertex import export_vertex_to_json
from pyflowline.formats.export_flowline import export_flowline_to_json
from pyflowline.algorithms.intersect.intersect_flowline_with_mesh import intersect_flowline_with_mesh

from pyflowline.algorithms.simplification.remove_returning_flowline import remove_returning_flowline
from pyflowline.algorithms.simplification.remove_duplicate_flowline import remove_duplicate_flowline
from pyflowline.algorithms.simplification.remove_duplicate_edge import remove_duplicate_edge
from pyflowline.algorithms.direction.correct_flowline_direction import correct_flowline_direction
from pyflowline.algorithms.loop.remove_flowline_loop import remove_flowline_loop
from pyflowline.algorithms.split.find_flowline_vertex import find_flowline_vertex
from pyflowline.algorithms.split.find_flowline_confluence import find_flowline_confluence
from pyflowline.algorithms.split.split_flowline import split_flowline
from pyflowline.algorithms.split.split_flowline_to_edge import split_flowline_to_edge

from pyflowline.algorithms.merge.merge_flowline import merge_flowline

from pyflowline.algorithms.index.define_stream_order import define_stream_order
from pyflowline.algorithms.index.define_stream_segment_index import define_stream_segment_index

import cartopy.crs as ccrs

desired_proj = ccrs.Orthographic(central_longitude=-75, central_latitude=42, globe=None)

desired_proj = ccrs.PlateCarree()

pDate = datetime.datetime.today()
sDate_default = "{:04d}".format(pDate.year) + "{:02d}".format(pDate.month) + "{:02d}".format(pDate.day)

class CaseClassEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.float):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, pybasin):
            return obj.lBasinID
        
        if isinstance(obj, pyflowline):
            return obj.lFlowlineID
        
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
    dResolution=0.0
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
    aFlowline=list()
    
    def __init__(self, aParameter):
        
        if 'sFilename_model_configuration' in aParameter:
            self.sFilename_model_configuration    = aParameter[ 'sFilename_model_configuration']
        
        if 'sWorkspace_bin' in aParameter:
            self.sWorkspace_bin= aParameter[ 'sWorkspace_bin']
            
        if 'sWorkspace_data' in aParameter:
            self.sWorkspace_data = aParameter[ 'sWorkspace_data']
        
        if 'sWorkspace_project' in aParameter:
            self.sWorkspace_project= aParameter[ 'sWorkspace_project']
        
        if 'sWorkspace_output' in aParameter:
            self.sWorkspace_output = aParameter[ 'sWorkspace_output']
        
        if 'sRegion' in aParameter:
            self.sRegion               = aParameter[ 'sRegion']
        
        if 'sModel' in aParameter:
            self.sModel                = aParameter[ 'sModel']

        if 'iFlag_standalone' in aParameter:
            self.iFlag_standalone = int(aParameter['iFlag_standalone'])
    
        if 'iFlag_flowline' in aParameter:
            self.iFlag_flowline             = int(aParameter[ 'iFlag_flowline'])
        
        if 'iFlag_simplification' in aParameter:
            self.iFlag_simplification = int(aParameter['iFlag_simplification'])


        if 'iFlag_create_mesh' in aParameter:
            self.iFlag_create_mesh = int(aParameter['iFlag_create_mesh'])    
        
        if 'iFlag_save_mesh' in aParameter:
            self.iFlag_save_mesh             = int(aParameter[ 'iFlag_save_mesh'])

        if 'iFlag_rotation' in aParameter:
            self.iFlag_rotation = int(aParameter['iFlag_rotation'])

        if 'iFlag_intersect' in aParameter:
            self.iFlag_intersect = int(aParameter['iFlag_intersect'])
      
        if 'iFlag_use_mesh_dem' in aParameter:
            self.iFlag_use_mesh_dem = int(aParameter['iFlag_use_mesh_dem'])
        
        if 'nOutlet' in aParameter:
            self.nOutlet = int(aParameter['nOutlet'])

        if 'iFlag_global' in aParameter:
            self.iFlag_global             = int(aParameter[ 'iFlag_global'])
        
        if 'iFlag_multiple_outlet' in aParameter:
            self.iFlag_multiple_outlet             = int(aParameter[ 'iFlag_multiple_outlet'])    
               
        iCase_index = int(aParameter['iCase_index'])
        sCase_index = "{:03d}".format( iCase_index )
        sDate   = aParameter[ 'sDate']
        if sDate is not None:
            self.sDate= sDate
        else:
            self.sDate = sDate_default

        self.iCase_index =   iCase_index
        sCase = self.sModel  + self.sDate + sCase_index
        self.sCase = sCase

        #the model can be run as part of hexwatershed or standalone
        if self.iFlag_standalone ==1:
            sPath = str(Path(self.sWorkspace_output)  /  sCase)
            self.sWorkspace_output = sPath
        else:
            sPath = self.sWorkspace_output
        
        Path(sPath).mkdir(parents=True, exist_ok=True)

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
         
        
        self.dResolution = float(aParameter['dResolution']) 
        self.dResolution_meter = float(aParameter['dResolution_meter']) 

        
        self.dLongitude_left = float(aParameter['dLongitude_left']) 
        self.dLongitude_right = float(aParameter['dLongitude_right']) 
        self.dLatitude_bot = float(aParameter['dLatitude_bot']) 
        self.dLatitude_top = float(aParameter['dLatitude_top']) 
       
       
        self.sFilename_spatial_reference = aParameter['sFilename_spatial_reference']
        self.sFilename_dem = aParameter['sFilename_dem']

        if 'sFilename_mesh_netcdf' in aParameter:
            self.sFilename_mesh_netcdf = aParameter['sFilename_mesh_netcdf']

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
            

        self.sJob =  aParameter['sJob'] 

        

        #model generated files

   
        self.sFilename_mesh = os.path.join(str(Path(self.sWorkspace_output)  ) , sMesh_type + ".json" )
        
        

        self.sFilename_mesh_info= os.path.join(str(Path(self.sWorkspace_output)  ) , sMesh_type + "_mesh_info.json"  )
    
        self.sWorkspace_data_project = str(Path(self.sWorkspace_data ) / self.sWorkspace_project)

                
        return
        

    def convert_flowline_to_json(self):
        for pBasin in self.aBasin:            
            pBasin.convert_flowline_to_json()
            pass

    def plot(self, sVariable_in=None):
        if sVariable_in == 'mesh':
            self.plot_mesh()
        else:
            if sVariable_in == 'overlap':
                self.plot_mesh_with_flowline()
            else:
                for pBasin in self.aBasin:            
                    pBasin.plot(sVariable_in= sVariable_in)
                    pass
        
        return
    
    def plot_mesh(self):

        sWorkspace_output_case = self.sWorkspace_output

        sFilename_json  =  self.sFilename_mesh

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
        
        for pFeature_shapefile in pLayer:
            pGeometry_in = pFeature_shapefile.GetGeometryRef()
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


                polygon = mpatches.Polygon(aCoords_gcs[:,0:2], closed=True,  linewidth=1, \
                    alpha=0.8, edgecolor = 'black',facecolor='none', \
                        transform=ccrs.PlateCarree() )

                ax.add_patch(polygon)                   


        dDiff_lon = dLon_max - dLon_min
        dDiff_lat = dLat_max - dLat_min
    
        ax.set_extent([dLon_min  , dLon_max , dLat_min , dLat_max ])

        ax.coastlines()#resolution='110m')
        ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.3, linestyle='--')


        sDirname = os.path.dirname(sFilename_json)
        sFilename  = Path(sFilename_json).stem + '.png'
        sFilename_out = os.path.join(sDirname, sFilename)
        plt.savefig(sFilename_out, bbox_inches='tight')

        pDataset = pLayer = pFeature  = None   
        plt.show()   
        return

    def plot_mesh_with_flowline(self):
        sWorkspace_output_case = self.sWorkspace_output

        sFilename_mesh  =  self.sFilename_mesh

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
        
        for pFeature_shapefile in pLayer:
            pGeometry_in = pFeature_shapefile.GetGeometryRef()
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


                polygon = mpatches.Polygon(aCoords_gcs[:,0:2], closed=True,   linewidth=1, \
                    alpha=0.8, edgecolor = 'black',facecolor='none', \
                        transform=ccrs.PlateCarree() )

                ax.add_patch(polygon)                   


        dDiff_lon = dLon_max - dLon_min
        dDiff_lat = dLat_max - dLat_min

        #plot flowline now
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

    
        ax.set_extent([dLon_min  , dLon_max , dLat_min , dLat_max ])

        ax.coastlines()#resolution='110m')
        ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.3, linestyle='--')


        sDirname = os.path.dirname(sFilename_mesh)
        sFilename  = Path(sFilename_mesh).stem + '_flowline.png'
        sFilename_out = os.path.join(sDirname, sFilename)
        plt.savefig(sFilename_out, bbox_inches='tight')

        pDataset = pLayer = pFeature  = None   
        plt.show()   
        return

    def preprocess_flowline(self):
        aFlowline_out = list()   #store all the flowline
        if self.iFlag_simplification == 1: 
            for pBasin in self.aBasin:
                aFlowline_basin = pBasin.preprocess_flowline()                
                aFlowline_out = aFlowline_out + aFlowline_basin

            self.aFlowline = aFlowline_out
        return aFlowline_out
    
    def create_mesh(self):        
        iFlag_global =  self.iFlag_global
        iMesh_type = self.iMesh_type
        iFlag_save_mesh = self.iFlag_save_mesh
        iFlag_rotation = self.iFlag_rotation
        dResolution = self.dResolution
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
                dResolution = meter_to_degree(dResolution_meter, dLatitude_mean)


            dX_left = dOriginX
            dY_top = dOriginY
        else:
            pass

        if iMesh_type ==1: #hexagon

            #hexagon edge
            dResolution_meter = degree_to_meter(dLatitude_mean, dResolution )
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
                    dResolution_meter = degree_to_meter(dLatitude_mean, dResolution)
                    dArea = np.power(dResolution_meter,2.0)
                    dLatitude_top    = self.dLatitude_top   
                    dLatitude_bot    = self.dLatitude_bot   
                    dLongitude_left  = self.dLongitude_left 
                    dLongitude_right = self.dLongitude_right
                    ncolumn= int( (dLongitude_right - dLongitude_left) / dResolution )
                    nrow= int( (dLatitude_top - dLatitude_bot) / dResolution )
                    aLatlon = create_latlon_mesh(dLongitude_left, dLatitude_bot, dResolution, ncolumn, nrow, \
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
                                dLatitude_top, dLatitude_bot, dLongitude_left, dLongitude_right,\
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
    
    def intersect_flowline_with_mesh(self):
        iMesh_type = self.iMesh_type
        iFlag_intersect = self.iFlag_intersect
        sWorkspace_output = self.sWorkspace_output
        nOutlet = self.nOutlet
        sFilename_mesh=self.sFilename_mesh
        aMesh, pSpatial_reference_mesh = read_mesh_json(sFilename_mesh)
        aCell = list()
        aCell_intersect = list()
        aFlowline = list()   #store all the flowline
        aOutletID = list()
        aBasin = list()
        if iFlag_intersect == 1:
            for pBasin in self.aBasin:
                pBasin.intersect_flowline_with_mesh(iMesh_type,sFilename_mesh)
                aFlowline = aFlowline + pBasin.aFlowline_basin
                aBasin.append(pBasin)
                aOutletID.append(pBasin.lCellID_outlet)
            
            #save basin info
            sPath = os.path.dirname(self.sFilename_basins)
            sName = str(Path(self.sFilename_basins).stem ) + '_new.json'
            sFilename_configuration  =  os.path.join( sPath  , sName)
            with open(sFilename_configuration, 'w', encoding='utf-8') as f:
                sJson = json.dumps([json.loads(ob.tojson()) for ob in aBasin],\
                    sort_keys=True, \
                    indent = 4)        

                f.write(sJson)    
                f.close()

            return aCell, aCell_intersect, aFlowline, aOutletID

        return


    def tojson(self):
        sJson = json.dumps(self.__dict__, \
            sort_keys=True, \
                indent = 4, \
                    ensure_ascii=True, \
                        cls=CaseClassEncoder)
        return sJson