import sys, os, stat
import numpy as np
from pathlib import Path
import json

from pyearth.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex
from pyearth.visual.barplot.barplot_data_with_reference import barplot_data_with_reference

sRegion = 'Susquehanna'

aResolution = ['50km', '10km', '5km']
nResolution =len(aResolution)

iSize_x = 12
iSize_y =  9
iDPI = 150
nData= 6

#read data
nCase = 10
sDate='20220201'
aColor= np.array(create_diverge_rgb_color_hex( 6 ))
sWorkspace_output = '/compyfs/liao313/04model/pyflowline/susquehanna/'

aLinestyle = np.array([  'dotted',  'dotted'  , 'solid', 'solid', 'solid', 'dashdot'])
aMarker=  np.array([  '+',  '^'  , 'o', 'p', 'd', '*'])

aLabel_legend= np.array(['NHDPlus HR','Simplified','Latlon','Square','Hexagon','MPAS'])
aHatch = np.array([ '.',   '*', '+', '|', '-', 'o'])

for r in range(nResolution):
    aData_full = np.full( (7,5),0.0,dtype=float )
    aCase = np.arange(3) * 3 +1 + r
    aCase= np.append(aCase,10)
    iFlag_sim =0
    for i in range( len( aCase) ):
        iCase_index = aCase[i]
        sCase_index = "{:03d}".format( iCase_index )
        sWorkspace_output_case = sWorkspace_output + '/' + 'pyflowline' + sDate + sCase_index
        if iFlag_sim ==0:
            sFilename_json = sWorkspace_output_case + '/' + '001' + '/'+ 'flowline_simplified_info.json'
            with open(sFilename_json) as json_file:
                data_sim = json.load(json_file)  
                for j in range(len(data_sim)):
                    aData_full[j,0] = float(data_sim[j]['dSinuosity'])
                iFlag_sim = 1
            
    
        sFilename_json = sWorkspace_output_case + '/' + '001' + '/'+ 'flowline_conceptual_info.json'
        with open(sFilename_json) as json_file:
            data_con = json.load(json_file) 
            for j in range(len(data_con)):
                    aData_full[j,i+1] = float(data_con[j]['dSinuosity'])
                

    sResolution = aResolution[r]
    y_label = r'River sinuosity (ratio)'
    sTitle = r'Susquehanna river basin '  
    sFilename_out= '/people/liao313/data/hexwatershed/susquehanna/river_length_by_segment_' + sResolution + '.png'

    aIndex = np.arange(7) + 1
    
    aSegment = ['Seg ' +  "{:0d}".format( x )  for x in aIndex]

    sFormat_x = ''

    sFormat_y = '%.1f'
    aSubset_index = np.arange(5) + 1

    aReference_in= []
    aData_full=np.transpose(aData_full)
    barplot_data_with_reference(aData_full, \
                 aSegment, \
                 aLabel_legend[aSubset_index],\
                 sFilename_out,\
                     aReference_in,\
                 dMax_y_in = 4,\
                 dMin_y_in = 0.0,\
                 sFormat_y_in = sFormat_y,
                 sLabel_y_in= y_label,\
                 ncolumn_in= 3,\
                     aLinestyle_in = aLinestyle[aSubset_index],\
                 aColor_in= aColor[aSubset_index],\
                     aMarker_in = aMarker[aSubset_index],\
                 aHatch_in = aHatch[aSubset_index],\
                 sTitle_in = sTitle)

print('finished')
