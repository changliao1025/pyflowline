from abc import ABCMeta, abstractmethod
import datetime
from pyearth.system.define_global_variables import *
pDate = datetime.datetime.today()
sDate_default = "{:04d}".format(pDate.year) + "{:02d}".format(pDate.month) + "{:02d}".format(pDate.day)

class flowlinecase(object):
    __metaclass__ = ABCMeta
    iCase_index= 0
    sMesh_type = 1
    iFlag_standalone=1
    iFlag_multiple = 0
    iFlag_use_mpas_dem=0
    iFlag_simplification = 1 #user can turn on/off
    iFlag_create_mesh=1
    iFlag_intersect = 1
    iFlag_disconnected =0
    iFlag_dam=0
    iFlag_rotation=0

    nOutlet = 1 #by default , there shoule ne only one ouelet

    dResolution=0.0
    dResolution_meter=0.0
    dThreshold_small_river=0.0

    dLongitude_left = -79.44374
    dLongitude_right = -74.24774 

    dLatitude_bot = 39.00 #,1399152.687,1978258.386
    dLatitude_top = 43.00334 # ,1748363.409,2424316.881


    sFilename_model_configuration=''

    sWorkspace_data=''
 
    
    sWorkspace_project=''
    
    sWorkspace_output=''
    
    
    sRegion=''
    sModel=''
    iMesh_type ='hexagon'

    sCase=''
    sDate=''
    

    sFilename_spatial_reference=''
    sFilename_dem=''
    sFilename_flowline_raw=''
    sFilename_flowline_filter=''
    sFilename_dam=''
    sFilename_flowline_topo=''
    #before intersect
    sFilename_flowline_segment_order_before_intersect=''
    sFilename_flowline_segment_index_before_intersect=''

    #intersect
    sFilename_mesh=''
    sFilename_mesh_info=''
    sFilename_mesh_netcdf=''
    sFilename_flowline_intersect = ''
    #after intersect
    sFilename_flowline_simplified_after_intersect=''
    sFilename_vertex_without_confluence_after_intersect=''
    flowline_split_by_point_after_intersect=''
    
    def __init__(self, aParameter):
        self.sFilename_model_configuration    = aParameter[ 'sFilename_model_configuration']

        self.sWorkspace_bin= aParameter[ 'sWorkspace_bin']
        self.sWorkspace_data = aParameter[ 'sWorkspace_data']
        
        self.sWorkspace_project= aParameter[ 'sWorkspace_project']
        self.sWorkspace_output = aParameter[ 'sWorkspace_output']
        
        self.sRegion               = aParameter[ 'sRegion']
        self.sModel                = aParameter[ 'sModel']


        if 'iFlag_standalone' in aParameter:
            self.iFlag_standalone = int(aParameter['iFlag_standalone'])
    
        if 'iFlag_multiple' in aParameter:
            self.iFlag_multiple = int(aParameter['iFlag_multiple'])
        
        if 'iFlag_simplification' in aParameter:
            self.iFlag_simplification = int(aParameter['iFlag_simplification'])


        if 'iFlag_create_mesh' in aParameter:
            self.iFlag_create_mesh = int(aParameter['iFlag_create_mesh'])    

        if 'iFlag_intersect' in aParameter:
            self.iFlag_intersect = int(aParameter['iFlag_intersect'])


        if 'iFlag_dam' in aParameter:
            self.iFlag_dam = int(aParameter['iFlag_dam'])

        if 'iFlag_use_mpas_dem' in aParameter:
            self.iFlag_use_mpas_dem = int(aParameter['iFlag_use_mpas_dem'])
        
        if 'nOutlet' in aParameter:
            self.nOutlet = int(aParameter['nOutlet'])
        
        
               
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
            sPath = self.sWorkspace_output + slash + sCase
            self.sWorkspace_output = sPath
        else:
            sPath = self.sWorkspace_output
        
        Path(sPath).mkdir(parents=True, exist_ok=True)

        self.sMesh_type =  aParameter['sMesh_type']
        
        sMesh_type = self.sMesh_type
        if sMesh_type =='hexagon': #hexagon
            self.iMesh_type = 1
        else:
            if sMesh_type =='square': #sqaure
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
         
     
        self.iFlag_disconnected =  int(aParameter['iFlag_disconnected'])

        self.iFlag_rotation = int(aParameter['iFlag_rotation'])
        
        self.dResolution = float(aParameter['dResolution']) 
        self.dResolution_meter = float(aParameter['dResolution_meter']) 

        self.dThreshold_small_river =  float(aParameter['dThreshold_small_river']) 

        
        
        self.dLongitude_left = float(aParameter['dLongitude_left']) 
        self.dLongitude_right = float(aParameter['dLongitude_right']) 
        self.dLatitude_bot = float(aParameter['dLatitude_bot']) 
        self.dLatitude_top = float(aParameter['dLatitude_top']) 
        

        if self.iFlag_multiple == 0:
            self.dLon_outlet = float(aParameter['dLon_outlet']) 
            self.dLat_outlet = float(aParameter['dLat_outlet']) 
        else:
            #use a file name to store outlet locations
            self.sFilename_outlet =  aParameter['sFilename_outlet']
            pass

        self.sFilename_spatial_reference = aParameter['sFilename_spatial_reference']
        self.sFilename_dem = aParameter['sFilename_dem']

        if 'sFilename_mesh_netcdf' in aParameter:
            self.sFilename_mesh_netcdf = aParameter['sFilename_mesh_netcdf']

        self.sFilename_flowline_filter = aParameter['sFilename_flowline_filter']

        if 'sFilename_dam' in aParameter:
            self.sFilename_dam = aParameter['sFilename_dam']

        if 'sFilename_flowline_topo' in aParameter:
            self.sFilename_flowline_topo = aParameter['sFilename_flowline_topo']

        if 'sFilename_flowline_raw' in aParameter:
            self.sFilename_flowline_raw = aParameter['sFilename_flowline_raw']

        #model generated files
        self.sFilename_mesh = self.sWorkspace_output + slash  + sMesh_type + ".shp"
        
        self.sFilename_flowline_segment_index_before_intersect = self.sWorkspace_output + slash + 'flowline_segment_index_before_intersect.shp'
        self.sFilename_flowline_segment_order_before_intersect = self.sWorkspace_output + slash + 'flowline_segment_order_before_intersect.shp'


        self.sFilename_mesh_info= self.sWorkspace_output + slash + sMesh_type + "_mesh_info.json"  
        
        self.sFilename_flowline_intersect  = self.sWorkspace_output + slash + 'flowline_intersect.shp'
        self.sJob =  aParameter['sJob'] 

        self.sWorkspace_data_project = self.sWorkspace_data +  slash + self.sWorkspace_project

                
        return
        