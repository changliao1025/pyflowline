import os 
import sys #used to add system path
from jdcal import gcal2jd, jd2gcal
import datetime



from pyearth.system.define_global_variables import *
from pyearth.toolbox.reader.parse_xml_file import parse_xml_file
from pyearth.toolbox.reader.read_configuration_file import read_configuration_file




pDate = datetime.datetime.today()
sDate_default = "{:04d}".format(pDate.year) + "{:02d}".format(pDate.month) + "{:02d}".format(pDate.day)

def pystream_read_model_configuration_file(sFilename_configuration_in,\
     iCase_index_in=None, \
         sJob_in=None,\
          iFlag_mode_in=None, \
         aVariable_in = None, \
             aValue_in = None, \
                 sDate_in = None):


    config = parse_xml_file(sFilename_configuration_in)
    sModel = config['sModel']  
    sRegion = config['sRegion']

    sWorkspace_data=  config['sWorkspace_data']
    sWorkspace_scratch=  config['sWorkspace_scratch']
    


    
    if iFlag_mode_in is not None:
        iFlag_mode = iFlag_mode_in
    else:
        iFlag_mode = 1
    

   

    if sDate_in is not None:
        sDate = sDate_in
    else:
        sDate = config['sDate']
        pass

    if iCase_index_in is not None:        
        iCase_index = iCase_index_in
    else:       
        iCase_index = int( config['iCase_index'])

    

    config['iCase_index'] = iCase_index
    config['sDate'] = sDate
   
   
    

    #based on global variable, a few variables are calculate once
    #calculate the modflow simulation period
    #https://docs.python.org/3/library/datetime.html#datetime-objects
   
   
    
    #data
    
    #simulation
    
    
   
    
    return config