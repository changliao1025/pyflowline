import os
from pathlib import Path
from abc import ABCMeta, abstractmethod

import datetime
import json
from pyflowline.classes.basin import pybasin
pDate = datetime.datetime.today()
sDate_default = "{:04d}".format(pDate.year) + "{:02d}".format(pDate.month) + "{:02d}".format(pDate.day)

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
        for pBasin in self.aBasin:            
            pBasin.plot(sVariable_in= sVariable_in)
            pass

    def preprocess_flowline(self):
        aFlowline_out = list()   #store all the flowline
        if self.iFlag_simplification == 1: 
            for pBasin in self.aBasin:
                aFlowline_basin = pBasin.preprocess_flowline()
                
                aFlowline_out = aFlowline_out + aFlowline_basin

            self.aFlowline = aFlowline_out
        return aFlowline_out
