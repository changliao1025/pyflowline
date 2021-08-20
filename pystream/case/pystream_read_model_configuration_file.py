import os 
import sys #used to add system path
from jdcal import gcal2jd, jd2gcal
import datetime
import json


from pyearth.system.define_global_variables import *

from pystream.case.pycase import streamcase




pDate = datetime.datetime.today()
sDate_default = "{:04d}".format(pDate.year) + "{:02d}".format(pDate.month) + "{:02d}".format(pDate.day)

def pystream_read_model_configuration_file(sFilename_configuration_in,\
     iCase_index_in=None, \
         dResolution_in = None,\
         dResolution_meter_in = None,\
         sJob_in=None,\
         aVariable_in = None, \
             aValue_in = None, \
                 sDate_in = None,\
                     sWorkspace_output_in = None):


    
    # Opening JSON file
    with open(sFilename_configuration_in) as json_file:
        data = json.load(json_file)   
    
    if iCase_index_in is not None:        
        iCase_index = iCase_index_in
    else:       
        iCase_index = int( data['iCase_index'])

    if sDate_in is not None:
        sDate = sDate_in
    else:
        sDate = data['sDate']
        pass

    if dResolution_in is not None:
        dResolution = dResolution_in
    else:
        dResolution = data['dResolution']
        pass

    if dResolution_meter_in is not None:
        dResolution_meter = dResolution_meter_in
    else:
        dResolution_meter = data['dResolution_meter']
        pass

    if sWorkspace_output_in is not None:
        sWorkspace_output = sWorkspace_output_in
    else:
        sWorkspace_output = data['sWorkspace_output']
        pass

    

    data['iCase_index'] = iCase_index
    data['dResolution'] = dResolution
    data['dResolution_meter'] = dResolution_meter

    data['sDate'] = sDate
    data['sWorkspace_output'] = sWorkspace_output
   
    

    #based on global variable, a few variables are calculate once
    #calculate the modflow simulation period
    #https://docs.python.org/3/library/datetime.html#datetime-objects
   
   
    
    #data
    
    #simulation
    
    oPystream = streamcase(data)
   
    
    return oPystream