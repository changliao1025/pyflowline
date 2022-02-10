import sys, os, stat
import numpy as np
from pathlib import Path


from pyearth.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex

from pyearth.visual.barplot.barplot_data_with_reference import barplot_data_with_reference

sRegion = 'Susquehanna'

aResolution = ['50km', '10km', '5km']
nResolution =len(aResolution)



iFlag_outlet =1

if iFlag_outlet==1:    

    a3 = np.array([885.0, 905.2, 925.3])   
    a4 = np.array([821.7, 952, 1051.5])
    a5 = np.array([806.5, 1010.7, 1021.5])
    
    
    y10 = 1084 #nhd
    y20 = 971.9 #sim
    y60 = 859.9 #mpas
    y_label = r'River length (ratio)'
    sTitle = r'Susquehanna River Basin'
    sFilename_out= '/people/liao313/data/hexwatershed/susquehanna/river_length.png'


a1 =np.array([y10, y10, y10]) 
a2 =np.array([y20, y20, y20]) 
a6 =np.array([y60, y60, y60]) 

iSize_x = 12
iSize_y =  9
iDPI = 150
nData= 6

aLinestyle = [  'dotted',  'dotted'  , 'solid', 'solid', 'solid', 'dashdot']

aMarker=  [  '+',  '^'  , 'o', 'p', 'd', '*']
aColor= create_diverge_rgb_color_hex(nData )
aLabel_legend= ['NHDPlus HR','Simplified','Latlon','Square','Hexagon','MPAS']
aHatch = [ '.',   '*', '+', '|', '-', 'o']

sFormat_x = ''

sFormat_y = '%.1f'

#need to transpose
#a = np.array([  a3, a4, a5])
#a3, a4 , a5= np.transpose(a)
aData= np.array( [ a1, a2, a3, a4, a5, a6]) 
aData = aData / (a1)

aReference_in= [0,1]
barplot_data_with_reference(aData, \
             aResolution, \
             aLabel_legend,\
             sFilename_out,\
                 aReference_in,\
             dMax_y_in = 1.5,\
             dMin_y_in = 0.0,\
             sFormat_y_in = sFormat_y,
             sLabel_y_in= y_label,\
             ncolumn_in= 3,\
                 aLinestyle_in = aLinestyle,\
             aColor_in= aColor,\
                 aMarker_in = aMarker,\
             aHatch_in = aHatch,\
             sTitle_in = sTitle)

print('finished')
