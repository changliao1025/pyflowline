import os
from pathlib import Path
from abc import ABCMeta
import json
from json import JSONEncoder
import datetime
import numpy as np
from osgeo import ogr, osr
from pyflowline.classes.mpas import pympas
from pyflowline.classes.hexagon import pyhexagon
from pyflowline.classes.latlon import pylatlon
from pyflowline.classes.square import pysquare
from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.basin import pybasin
from pyflowline.classes.flowline import pyflowline
from pyflowline.classes.edge import pyedge
from pyflowline.formats.read_mesh import read_mesh_json
from pyflowline.algorithms.auxiliary.gdal_functions import reproject_coordinates
from pyflowline.algorithms.auxiliary.gdal_functions  import degree_to_meter
from pyflowline.algorithms.auxiliary.gdal_functions  import meter_to_degree
from pyflowline.algorithms.auxiliary.gdal_functions import retrieve_geotiff_metadata
from pyflowline.mesh.hexagon.create_hexagon_mesh import create_hexagon_mesh
from pyflowline.mesh.latlon.create_latlon_mesh import create_latlon_mesh
from pyflowline.mesh.square.create_square_mesh import create_square_mesh
from pyflowline.mesh.mpas.create_mpas_mesh import create_mpas_mesh
from pyflowline.mesh.tin.create_tin_mesh import create_tin_mesh

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
    iFlag_use_shapefile_extent=1
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
    sWorkspace_input=''      
    sWorkspace_output=''    
    #sWorkspace_output_case=''    
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

    
    def __init__(self, aConfig_in,\
            iFlag_standalone_in= None,\
            sModel_in = None,\
            sWorkspace_output_in = None):
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
        
        if 'iFlag_global' in aConfig_in:
            self.iFlag_global             = int(aConfig_in[ 'iFlag_global'])
        
        if 'iFlag_multiple_outlet' in aConfig_in:
            self.iFlag_multiple_outlet             = int(aConfig_in[ 'iFlag_multiple_outlet'])  
        
        if 'iFlag_simplification' in aConfig_in:
            self.iFlag_simplification = int(aConfig_in['iFlag_simplification'])


        if 'iFlag_create_mesh' in aConfig_in:
            self.iFlag_create_mesh = int(aConfig_in['iFlag_create_mesh'])    
        
        if 'iFlag_save_mesh' in aConfig_in:
            self.iFlag_save_mesh             = int(aConfig_in[ 'iFlag_save_mesh'])
        
        if 'iFlag_use_mesh_dem' in aConfig_in:
            self.iFlag_use_mesh_dem = int(aConfig_in['iFlag_use_mesh_dem'])

        if 'iFlag_use_shapefile_extent' in aConfig_in:
            self.iFlag_use_shapefile_extent = int(aConfig_in['iFlag_use_shapefile_extent'])    

        if 'iFlag_rotation' in aConfig_in:
            self.iFlag_rotation = int(aConfig_in['iFlag_rotation'])
        

        if 'iFlag_intersect' in aConfig_in:
            self.iFlag_intersect = int(aConfig_in['iFlag_intersect'])

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

        if 'sFilename_spatial_reference' in aConfig_in:
            self.sFilename_spatial_reference = aConfig_in['sFilename_spatial_reference']

        if 'sFilename_dem' in aConfig_in:
            self.sFilename_dem = aConfig_in['sFilename_dem']

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
                        if sMesh_type =='tin': #tin
                            self.iMesh_type = 5
                        else:
                            print('Unsupported mesh type?')
        
        #the model can be run as part of hexwatershed or standalone
        if self.iFlag_standalone == 1:
            #in standalone case, will add case information and update output path
            sPath = str(Path(self.sWorkspace_output)  /  sCase)
            self.sWorkspace_output = sPath
        else:
            #use specified output path, also do not add output or input tag
            sPath = self.sWorkspace_output
        
        Path(sPath).mkdir(parents=True, exist_ok=True)

        self.aBasin = list()
        if self.iFlag_flowline == 1:
            if 'sFilename_basins' in aConfig_in:
                self.sFilename_basins = aConfig_in['sFilename_basins']

                if os.path.isfile(self.sFilename_basins):
                    pass
                else:
                    print('This basin configuration file does not exist: ', self.sFilename_basins )
                    exit()

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
            spatial_reference_target = osr.SpatialReference()  
            spatial_reference_target.ImportFromEPSG(4326)
            if self.iFlag_use_shapefile_extent==1:
                pDriver_shapefile = ogr.GetDriverByName('Esri Shapefile')
                pDataset_shapefile = pDriver_shapefile.Open(self.sFilename_spatial_reference, 0)
                pLayer_shapefile = pDataset_shapefile.GetLayer(0)
                pSpatial_reference = pLayer_shapefile.GetSpatialRef()
                (dOriginX, dX_right, dY_bot, dOriginY) = pLayer_shapefile.GetExtent() # not in same order as ogrinfo
                dLongitude_left,  dLatitude_bot= reproject_coordinates(dOriginX, dY_bot,pSpatial_reference,    spatial_reference_target)
                dLongitude_right, dLatitude_top= reproject_coordinates(dX_right, dOriginY,pSpatial_reference,  spatial_reference_target)
                dLatitude_mean = 0.5 * (dLatitude_top + dLatitude_bot)
                pass
            else:
                dPixelWidth, dOriginX, dOriginY, nrow, ncolumn, pSpatialRef_dem, pProjection, pGeotransform\
                     = retrieve_geotiff_metadata(sFilename_dem)

                dY_bot = dOriginY - (nrow+1) * dPixelWidth
                dLongitude_left,  dLatitude_bot= reproject_coordinates(dOriginX, dY_bot,pSpatialRef_dem,    spatial_reference_target)
                dX_right = dOriginX + (ncolumn +1) * dPixelWidth

                dLongitude_right, dLatitude_top= reproject_coordinates(dX_right, dOriginY,pSpatialRef_dem,  spatial_reference_target)
                dLatitude_mean = 0.5 * (dLatitude_top + dLatitude_bot)
                pass

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
    
    def reconstruct_topological_relationship(self, aCell_raw):
        iMesh_type = self.iMesh_type
        iFlag_intersect = self.iFlag_intersect
        sWorkspace_output = self.sWorkspace_output
        nOutlet = self.nOutlet
        sFilename_mesh=self.sFilename_mesh
        self.aCell, pSpatial_reference_mesh = read_mesh_json(iMesh_type, sFilename_mesh)

        self.aCell = self.merge_cell_info(aCell_raw)
        
        aFlowline_conceptual = list()   #store all the flowline
        aCellID_outlet = list()
        aBasin = list()
        aCell_intersect=list()
        if iFlag_intersect == 1:
            for pBasin in self.aBasin:
                aCell_intersect_basin = pBasin.reconstruct_topological_relationship(iMesh_type,sFilename_mesh)
                aFlowline_conceptual = aFlowline_conceptual + pBasin.aFlowline_basin_conceptual
                aBasin.append(pBasin)
                aCellID_outlet.append(pBasin.lCellID_outlet)
                aCell_intersect = aCell_intersect + aCell_intersect_basin

                #update length?
            for pCell in self.aCell:
                for pCell2 in aCell_intersect:
                    if pCell2.lCellID == pCell.lCellID:
                        pCell.dLength_flowline = pCell2.dLength_flowline
            
            self.aFlowline_conceptual = aFlowline_conceptual
            self.aCellID_outlet = aCellID_outlet

            return self.aCell, aFlowline_conceptual, aCellID_outlet
        else:
            return None
            

        
   
    def merge_cell_info(self, aCell_raw):

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
        
    def analyze(self):
        for pBasin in self.aBasin:
            pBasin.analyze()
        return
        
    

    def setup(self):
        self.convert_flowline_to_json()
        return

    def run(self):
        #self.flowline_simplification()
        aCell = self.mesh_generation()
        aCell_out, a, b = self.reconstruct_topological_relationship(aCell)
        
        return aCell_out

   
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
        aSkip = ['aBasin', \
                'aFlowline_simplified','aFlowline_conceptual','aCellID_outlet',
                'aCell']      

        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)
        sJson = json.dumps(obj,\
            sort_keys=True, \
                indent = 4, \
                    ensure_ascii=True, \
                        cls=CaseClassEncoder)
        return sJson

    def export_config_to_json(self, sFilename_output_in = None):

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


    def export_basin_config_to_json(self, sFilename_output_in= None):
        if self.iFlag_standalone == 1:
            if sFilename_output_in is not None:
                sFilename_output = sFilename_output_in
            else:
                #use current output path
                sName = 'configuration_basin.json'
                sFilename_output  =  os.path.join( self.sWorkspace_output  , sName)                
            
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