from abc import ABCMeta, abstractmethod
import datetime
from pyearth.system.define_global_variables import *
pDate = datetime.datetime.today()
sDate_default = "{:04d}".format(pDate.year) + "{:02d}".format(pDate.month) + "{:02d}".format(pDate.day)

class streamcase(object):
    __metaclass__ = ABCMeta
    iCase_index= 0
    sMesh_type = 1
    
    iFlag_mode=0    
    iFlag_use_dem = 0
    iFlag_disconnected =0
    dResolution=0.0
    dResolution_meter=0.0
    dThreshold_small_river=0.0

    dLongitude_left = -79.44374
    dLongitude_right = -74.24774 

    dLatitude_bot = 39.00 #,1399152.687,1978258.386
    dLatitude_top = 43.00334 # ,1748363.409,2424316.881


    sFilename_model_configuration=''

    sWorkspace_data=''
    sWorkspace_scratch=''
    
    sWorkspace_project=''
    
    sWorkspace_simulation=''
    sWorkspace_simulation_case=''
    
    sRegion=''
    sModel=''
    iMesh_type ='hexagon'

    sCase=''
    sDate=''
    

    sFilename_spatial_reference=''
    sFilename_dem=''
    sFilename_flowlinw_raw=''
    #before intersect
    sFilename_flowline_segment_order_before_intersect=''
    sFilename_flowline_segment_index_before_intersect=''

    #intersect
    sFilename_mesh=''
    sFilename_mesh_netcdf=''
    sFilename_flowline_intersect = ''
    #after intersect
    sFilename_flowline_simplified_after_intersect=''
    sFilename_vertex_without_confluence_after_intersect=''
    flowline_split_by_point_after_intersect=''
    
    def __init__(self, aParameter):
        self.sFilename_model_configuration    = aParameter[ 'sFilename_model_configuration']

        self.sWorkspace_data = aParameter[ 'sWorkspace_data']
        self.sWorkspace_scratch    = aParameter[ 'sWorkspace_scratch']
        self.sWorkspace_project= aParameter[ 'sWorkspace_project']
        self.sWorkspace_bin= aParameter[ 'sWorkspace_bin']
        
        self.sRegion               = aParameter[ 'sRegion']
        self.sModel                = aParameter[ 'sModel']
       
    
        self.sWorkspace_simulation = self.sWorkspace_scratch + slash + '04model' + slash \
            + self.sModel + slash + self.sRegion +  slash + 'simulation'
        sPath = self.sWorkspace_simulation
        Path(sPath).mkdir(parents=True, exist_ok=True)

       
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
        
        self.sWorkspace_simulation_case = self.sWorkspace_simulation + slash + self.sMesh_type + slash + sCase
        sPath = self.sWorkspace_simulation_case
        Path(sPath).mkdir(parents=True, exist_ok=True)


       
        self.iFlag_simulation =  int(aParameter['iFlag_simulation']) 
        self.iFlag_mode =  int(aParameter['iFlag_mode']) 
        self.iFlag_use_dem =  int(aParameter['iFlag_use_dem']) 
        self.iFlag_disconnected =  int(aParameter['iFlag_disconnected'])
        self.dResolution = float(aParameter['dResolution']) 
        self.dResolution_meter = float(aParameter['dResolution_meter']) 

        self.dThreshold_small_river =  float(aParameter['dThreshold_small_river']) 

        
        
        self.dLongitude_left = float(aParameter['dLongitude_left']) 
        self.dLongitude_right = float(aParameter['dLongitude_right']) 
        self.dLatitude_bot = float(aParameter['dLatitude_bot']) 
        self.dLatitude_top = float(aParameter['dLatitude_top']) 

        self.dx_outlet = float(aParameter['dx_outlet']) 
        self.dy_outlet = float(aParameter['dy_outlet']) 

        self.sFilename_spatial_reference = aParameter['sFilename_spatial_reference']
        self.sFilename_dem = aParameter['sFilename_dem']

        self.sFilename_mesh_netcdf = aParameter['sFilename_mesh_netcdf']
        self.sFilename_flowlinw_raw = aParameter['sFilename_flowlinw_raw']


        self.sFilename_mesh = self.sWorkspace_simulation_case + slash + aParameter['sFilename_mesh']
        
        self.sFilename_flowline_segment_index_before_intersect = self.sWorkspace_simulation_case + slash + aParameter['sFilename_flowline_segment_index_before_intersect']
        self.sFilename_flowline_segment_order_before_intersect = self.sWorkspace_simulation_case + slash + aParameter['sFilename_flowline_segment_order_before_intersect']
        
        self.sFilename_flowline_intersect  = self.sWorkspace_simulation_case + slash + aParameter['sFilename_flowline_intersect']
        self.sJob =  aParameter['sJob'] 

        self.sWorkspace_data_project = self.sWorkspace_data +  slash + self.sWorkspace_project

                
        return
        