import os, sys

import numpy as np
import osr
from pyearth.gis.gdal.gdal_function import reproject_coordinates
from pyearth.toolbox.reader.text_reader_string import text_reader_string
from pyflowline.shared.vertex import pyvertex
from pyearth.system.define_global_variables import *

from pyflowline.format.read_flowline_shapefile import read_flowline_shapefile
from pyflowline.format.read_nhdplus_flowline_shapefile import read_nhdplus_flowline_shapefile_attribute
from pyflowline.format.read_nhdplus_flowline_shapefile import extract_nhdplus_flowline_shapefile_by_attribute
from pyflowline.format.read_nhdplus_flowline_shapefile import track_nhdplus_flowline

from pyflowline.format.export_flowline_to_shapefile import export_flowline_to_shapefile
from pyflowline.format.export_vertex_to_shapefile import export_vertex_to_shapefile

from pyflowline.algorithm.connect.connect_disconnect_flowline import connect_disconnect_flowline
from pyflowline.algorithm.direction.correct_flowline_direction import correct_flowline_direction

#merge
from pyflowline.algorithm.merge.merge_flowline import merge_flowline

#split
from pyflowline.algorithm.split.split_flowline import split_flowline
from pyflowline.algorithm.split.find_flowline_confluence import find_flowline_confluence
from pyflowline.algorithm.split.find_flowline_vertex import find_flowline_vertex


from pyflowline.algorithm.loop.remove_flowline_loop import remove_flowline_loop
#
from pyflowline.algorithm.simplification.remove_small_river import remove_small_river

from pyflowline.algorithm.index.define_stream_order import define_stream_order
from pyflowline.algorithm.index.define_stream_segment_index import define_stream_segment_index


def read_list_of_files(sFilename_in):
    aFilename = text_reader_string(sFilename_in)
    return aFilename
"""
prepare the flowline using multiple step approach
"""

def preprocess_flowline_op(opyflowline_in):
    
    #read shapefile and store information in the list
    iFlag_multiple = opyflowline_in.iFlag_multiple
    
    iFlag_dam = opyflowline_in.iFlag_dam
    iFlag_disconnected = opyflowline_in.iFlag_disconnected
    dThreshold = opyflowline_in.dThreshold_small_river

    sFilename_flowline_filter = opyflowline_in.sFilename_flowline_filter

    sWorkspace_output = opyflowline_in.sWorkspace_output

    if iFlag_multiple == 0:
        #only one outlet in this case
        aFlowline, pSpatialRef_pcs = read_flowline_shapefile(sFilename_flowline_filter)
        #we also need to save the spatial reference information for the output purpose
        pass
    else:
        nOutlet = opyflowline_in.nOutlet
        aFilename_flowline_filter = read_list_of_files(sFilename_flowline_filter)
        aFlowline = list()
        for i in range(nOutlet):
            #in this case, the sFilename_flowline_filter is a list of files
            sFilename_flowline_filter_dummy = aFilename_flowline_filter(i)

            aFlowline_dummy, pSpatialRef_pcs = read_flowline_shapefile( sFilename_flowline_filter_dummy )
            aFlowline = aFlowline + aFlowline_dummy
            pass


    


    if iFlag_dam ==1:
        sFilename_dam = opyflowline_in.sFilename_dam
        aData_dam = text_reader_string(sFilename_dam, iSkipline_in =1,cDelimiter_in=',' )
        sFilename_flowline_topo = opyflowline_in.sFilename_flowline_topo
        aData_flowline_topo = text_reader_string(sFilename_flowline_topo, iSkipline_in =1,cDelimiter_in=',' )

        aFromFlowline = aData_flowline_topo[:,1].astype(int).ravel()
        aToFlowline = aData_flowline_topo[:,2].astype(int).ravel()

        sFilename_flowline_raw = opyflowline_in.sFilename_flowline_raw
        aNHDPlusID_filter = read_nhdplus_flowline_shapefile_attribute(sFilename_flowline_filter)
        aNHDPlusID_raw = read_nhdplus_flowline_shapefile_attribute(sFilename_flowline_raw)
        ndam = len(aData_dam)
        aNHDPlusID_dams_headwater = list()
        aNHDPlusID_dams_nonheadwater = list()
        
        for i in range(0, ndam):
            print(i)
            dLon = float(aData_dam[i][1])
            dLat = float(aData_dam[i][0])
            sDam = aData_dam[i][4]            
            lNHDPlusID = int(aData_dam[i][5])
            aNHDPlusID_dams_headwater.append(lNHDPlusID)
            if i==2:
                print('debug')
            if lNHDPlusID in aNHDPlusID_filter:
                #remove by id
                for j in range(len(aFlowline)):
                    if aFlowline[j].lNHDPlusID == lNHDPlusID:
                        aFlowline.pop(j)
                        break
                pass
            else:                                
                aNHDPlusID_dam_nonheadwater = track_nhdplus_flowline(aNHDPlusID_filter, aFromFlowline, aToFlowline, lNHDPlusID)
                aNHDPlusID_filter = aNHDPlusID_filter + aNHDPlusID_dams_headwater+ aNHDPlusID_dam_nonheadwater  
                aNHDPlusID_dams_nonheadwater = aNHDPlusID_dams_nonheadwater + aNHDPlusID_dam_nonheadwater
        


        print('finished')
        aFlowline_dams_headwater = extract_nhdplus_flowline_shapefile_by_attribute(sFilename_flowline_raw, aNHDPlusID_dams_headwater )
        for i in range(len(aFlowline_dams_headwater)):
            aFlowline_dams_headwater[i].iFlag_dam = 1

        aFlowline_dams_nonheadwater = extract_nhdplus_flowline_shapefile_by_attribute(sFilename_flowline_raw, aNHDPlusID_dams_nonheadwater )

        aFlowline = aFlowline + aFlowline_dams_headwater + aFlowline_dams_nonheadwater

    #the flowline should not be in GCS because it cannot be used for distance directly
    pSpatialRef_gcs = osr.SpatialReference()
    pSpatialRef_gcs.ImportFromEPSG(4326)
    pSpatialRef_gcs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    sFilename_out = 'flowline_before_intersect.shp'
    sFilename_out = os.path.join(sWorkspace_output, sFilename_out)

    iFlag_projected = 0
    export_flowline_to_shapefile(iFlag_projected, aFlowline, pSpatialRef_gcs, sFilename_out)
    

    if iFlag_disconnected ==1:
        #need a better way to include this capability
         
        aThreshold = np.full(2, 300.0, dtype=float)
        #aFlowline = connect_disconnect_flowline(aFlowline, aVertex, aThreshold)
        sFilename_out = 'flowline_connect.json'
        sFilename_out = os.path.join(sWorkspace_output, sFilename_out)
        
        export_flowline_to_shapefile(iFlag_projected, aFlowline,pSpatialRef_gcs, sFilename_out)
    else:
        pass


    aVertex = find_flowline_vertex(aFlowline)
    sFilename_out = 'flowline_vertex_without_confluence_before_intersect.shp'
    sFilename_out = os.path.join(sWorkspace_output, sFilename_out)
    export_vertex_to_shapefile(iFlag_projected, aVertex,pSpatialRef_gcs, sFilename_out)

    aFlowline = split_flowline(aFlowline, aVertex)
    sFilename_out = 'flowline_split_by_point_before_intersect.shp'
    sFilename_out = os.path.join(sWorkspace_output, sFilename_out)
    export_flowline_to_shapefile(iFlag_projected, aFlowline,pSpatialRef_gcs, sFilename_out)

    #ues location to find outlet
  
    point= dict()   
    point['lon'] = opyflowline_in.dLon_outlet
    point['lat'] = opyflowline_in.dLat_outlet
    pVertex_outlet=pyvertex(point)

    aFlowline= correct_flowline_direction(aFlowline,  pVertex_outlet )

    pVertex_outlet = aFlowline[0].pVertex_end

    sFilename_out = 'flowline_direction_before_intersect.shp'
    sFilename_out = os.path.join(sWorkspace_output, sFilename_out)
    export_flowline_to_shapefile(iFlag_projected, aFlowline, pSpatialRef_gcs, sFilename_out)

    #step 4: remove loops

    aFlowline = remove_flowline_loop(aFlowline)    
    sFilename_out = 'flowline_loop_before_intersect.shp'
    sFilename_out = os.path.join(sWorkspace_output, sFilename_out)
    export_flowline_to_shapefile(iFlag_projected, aFlowline,pSpatialRef_gcs, sFilename_out)

    #using loop to remove small river, here we use 5 steps

    for i in range(3):
        sStep = "{:02d}".format(i+1)
        aFlowline = remove_small_river(aFlowline, dThreshold)
        sFilename_out = 'flowline_large_'+ sStep +'_before_intersect.shp'
        sFilename_out =os.path.join(sWorkspace_output, sFilename_out)
        export_flowline_to_shapefile(iFlag_projected, aFlowline, pSpatialRef_gcs, sFilename_out)
        

        aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity = find_flowline_confluence(aFlowline,  pVertex_outlet)
        sFilename_out = 'flowline_vertex_with_confluence_'+ sStep +'_before_intersect.shp'
        sFilename_out = os.path.join(sWorkspace_output, sFilename_out)
        export_vertex_to_shapefile(iFlag_projected, aVertex, pSpatialRef_gcs, sFilename_out, aAttribute_data=aConnectivity)

        aFlowline = merge_flowline( aFlowline,aVertex, pVertex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence  )  
        sFilename_out = 'flowline_merge_'+ sStep +'_before_intersect.shp'
        sFilename_out = os.path.join(sWorkspace_output, sFilename_out)
        export_flowline_to_shapefile(iFlag_projected, aFlowline, pSpatialRef_gcs, sFilename_out)

        if len(aFlowline) ==1:
            break


    #build segment index
    aFlowline, aStream_segment = define_stream_segment_index(aFlowline)
    sFilename_out = opyflowline_in.sFilename_flowline_segment_index_before_intersect
    sFilename_out = os.path.join(sWorkspace_output, sFilename_out)
    export_flowline_to_shapefile(iFlag_projected, aFlowline, pSpatialRef_gcs, sFilename_out, \
        aAttribute_data=[aStream_segment], aAttribute_field=['iseg'], aAttribute_dtype=['int'])

    #build stream order 
    aFlowline, aStream_order = define_stream_order(aFlowline)
    sFilename_out = opyflowline_in.sFilename_flowline_segment_order_before_intersect
    
    export_flowline_to_shapefile(iFlag_projected, aFlowline, pSpatialRef_gcs, sFilename_out, \
        aAttribute_data=[aStream_segment, aStream_order], aAttribute_field=['iseg','iord'], aAttribute_dtype=['int','int'])

    

    print('Finished')