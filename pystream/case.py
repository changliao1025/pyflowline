from abc import ABCMeta, abstractmethod
import datetime
from pyearth.system.define_global_variables import *
pDate = datetime.datetime.today()
sDate_default = "{:04d}".format(pDate.year) + "{:02d}".format(pDate.month) + "{:02d}".format(pDate.day)

class streamcase(object):
    __metaclass__ = ABCMeta
    iCase_index=0
    iMesh_type = 1
    
    iFlag_mode=0    


    sFilename_model_configuration=''

    sWorkspace_data=''
    sWorkspace_scratch=''
    
    sWorkspace_project=''
    
    sWorkspace_simulation=''
    sWorkspace_simulation_case=''
    
    sRegion=''
    sModel=''

    sCase=''
    sDate=''
    

    

    #before intersect

    #intersect
    sFilename_mesh=''
    sFilename_intersect = ''
    #after intersect
    sFilename_flowline_simplified_after_intersect=''
    sFilename_vertex_without_confluence_after_intersect=''
    flowline_split_by_point_after_intersect=''
    
    def __init__(self, aParameter):
        self.sFilename_model_configuration    = aParameter[ 'sFilename_model_configuration']

        self.sWorkspace_data = aParameter[ 'sWorkspace_data']
       
        self.sWorkspace_scratch    = aParameter[ 'sWorkspace_scratch']
        sWorkspace_scratch = self.sWorkspace_scratch
        
        self.sRegion               = aParameter[ 'sRegion']
        self.sModel                = aParameter[ 'sModel']
       
        self.sWorkspace_project= aParameter[ 'sWorkspace_project']
        self.sWorkspace_bin= aParameter[ 'sWorkspace_bin']

        self.sWorkspace_simulation = sWorkspace_scratch + slash + '04model' + slash \
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
        sCase = self.sModel + self.sDate + sCase_index
        self.sCase = sCase
        

        self.sWorkspace_simulation_case = self.sWorkspace_simulation + slash + sCase
        sPath = self.sWorkspace_simulation_case
        Path(sPath).mkdir(parents=True, exist_ok=True)
      
       
        self.iFlag_simulation =  int(aParameter['iFlag_simulation']) 
     

        self.iFlag_mode =  int(aParameter['iFlag_mode']) 
      

        self.sJob =  aParameter['sJob'] 

        self.sWorkspace_data_project = self.sWorkspace_data +  slash + self.sWorkspace_project

                
        return
        