"""
this module will be used to define all the global variables 
"""

import platform  
from pathlib import Path
import getpass
sPlatform_os = platform.system()

sUsername = getpass.getuser()

sWorkspace_home = str(Path.home())

if sPlatform_os == 'Windows':  #windows
    slash = '\\'
    sMachine ='None'
       
    sWorkspace_scratch = 'C:'    
else:  #linux or unix
    slash = '/'    
   
    if (sPlatform_os == 'Linux'):
        sWorkspace_scratch =  sWorkspace_home + slash + 'scratch'            
        
    else:
        if (sPlatform_os == 'Darwin'):
            sMachine ='mac'
            sWorkspace_scratch =  sWorkspace_home + slash + 'scratch'  
        else:
            pass

#data file type
sExtension_txt = '.txt'
sExtension_envi = '.dat'
sExtension_tiff = '.tif'
sExtension_header ='.hdr'
sExtension_netcdf = '.nc'
sExtension_shapefile = '.shp'
sExtension_json = '.json'

#graphics format

sExtension_png = '.png'
sExtension_jpg = '.jpg'
sExtension_ps = '.ps'
sExtension_vtk = '.vtk'


#constant values
missing_value = -9999.0

nmonth = 12 #be careful with this one

#unit conversion

iMonth_start = 1
iMonth_end = 12
mms2mmd = 24 * 3600.0
feet2meter = 0.3048
inch2mm = 25.4
cms2cmd = 24 * 3600
