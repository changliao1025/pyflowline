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
nSegment=7

#read data
nCase = 10
sDate='20220201'
aColor= np.array(create_diverge_rgb_color_hex( 6 ))
sWorkspace_output = '/compyfs/liao313/04model/pyflowline/susquehanna/'

aLinestyle = np.array([  'dotted',  'dotted'  , 'solid', 'solid', 'solid', 'dashdot'])
aMarker=  np.array([  '+',  '^'  , 'o', 'p', 'd', '*'])

aLabel_legend= np.array(['NHDPlus HR','Simplified','Latlon','Square','Hexagon','MPAS'])
aHatch = np.array([ '.',   '*', '+', '|', '-', 'o'])

for iSeg in range(nSegment):
    sSegment = "{:0d}".format( iSeg+1  )
    aData_full = np.full((5,3),0.0, dtype=float )     
    
    iFlag_sim = 0
    for iRes in range( nResolution ):
        #simplified
        sCase_index = "{:03d}".format( 1 )
        sWorkspace_output_case = sWorkspace_output + '/' + 'pyflowline' + sDate + sCase_index
        if iFlag_sim ==0:
            sFilename_json = sWorkspace_output_case + '/' + '001' + '/'+ 'flowline_simplified_info.json'
            with open(sFilename_json) as json_file:
                data_sim = json.load(json_file)  
                data_seg = data_sim[iSeg]                
                aData_full[0,:] = float(data_seg['dSinuosity'])
                iFlag_sim = 1
            
        #latlon
        sCase_index = "{:03d}".format( iRes+7  )
        sWorkspace_output_case = sWorkspace_output + '/' + 'pyflowline' + sDate + sCase_index
        sFilename_json = sWorkspace_output_case + '/' + '001' + '/'+ 'flowline_conceptual_info.json'
        with open(sFilename_json) as json_file:
            data_con = json.load(json_file) 
            data_seg = data_con[iSeg]           
            aData_full[1,iRes] = float(data_seg['dSinuosity'])

        #square 
        sCase_index = "{:03d}".format( iRes +4 )
        sWorkspace_output_case = sWorkspace_output + '/' + 'pyflowline' + sDate + sCase_index
        sFilename_json = sWorkspace_output_case + '/' + '001' + '/'+ 'flowline_conceptual_info.json'
        with open(sFilename_json) as json_file:
            data_con = json.load(json_file) 
            data_seg = data_con[iSeg]           
            aData_full[2,iRes] = float(data_seg['dSinuosity'])
        
        #hexagon
        sCase_index = "{:03d}".format( iRes+1 )
        sWorkspace_output_case = sWorkspace_output + '/' + 'pyflowline' + sDate + sCase_index
        sFilename_json = sWorkspace_output_case + '/' + '001' + '/'+ 'flowline_conceptual_info.json'
        with open(sFilename_json) as json_file:
            data_con = json.load(json_file) 
            data_seg = data_con[iSeg]           
            aData_full[3,iRes] = float(data_seg['dSinuosity'])

        #mpas 
        sCase_index = "{:03d}".format( 10 )
        sWorkspace_output_case = sWorkspace_output + '/' + 'pyflowline' + sDate + sCase_index
        sFilename_json = sWorkspace_output_case + '/' + '001' + '/'+ 'flowline_conceptual_info.json'
        with open(sFilename_json) as json_file:
            data_con = json.load(json_file) 
            data_seg = data_con[iSeg]           
            aData_full[4,iRes] = float(data_seg['dSinuosity'])

    
    y_label = r'River sinuosity (ratio)'
    sTitle = r'Susquehanna river basin '  
    sFilename_out= '/people/liao313/data/hexwatershed/susquehanna/river_sinuosity_by_segment_' + sSegment + '.png'

    aIndex = np.arange(7) + 1
    
    #aSegment = ['Seg ' +  "{:0d}".format( x )  for x in aIndex]

    sFormat_x = ''

    sFormat_y = '%.1f'
    aSubset_index = np.arange(5) + 1

    aReference_in= [0]
    #aData_full=np.transpose(aData_full)
    barplot_data_with_reference(aData_full, \
                 aResolution, \
                 aLabel_legend[aSubset_index],\
                 sFilename_out,\
                     aReference_in,\
                 dMax_y_in = 4,\
                 dMin_y_in = 0.0,\
                 sFormat_y_in = sFormat_y,
                 sLabel_y_in= y_label,\
                     sLabel_info_in = 'Segment: '+sSegment, \
                 ncolumn_in= 3,\
                     aLinestyle_in = aLinestyle[aSubset_index],\
                 aColor_in= aColor[aSubset_index],\
                     aMarker_in = aMarker[aSubset_index],\
                 aHatch_in = aHatch[aSubset_index],\
                 sTitle_in = sTitle)

print('finished')
