import os, sys

import numpy as np
import osr

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


def preprocess_flowline_op(oPyflowline_in):
    
    #read shapefile and store information in the list
    iFlag_simplification = oPyflowline_in.iFlag_simplification 
    #the flowline should not be in GCS because it cannot be used for distance directly
    pSpatialRef_gcs = osr.SpatialReference()
    pSpatialRef_gcs.ImportFromEPSG(4326)
    pSpatialRef_gcs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

    nOutlet = oPyflowline_in.nOutlet
    if iFlag_simplification == 1: 
        
        aFlowline = list()        #store all the flowline
        for i in range(nOutlet):
            sBasin =  "{:03d}".format(i+1)
            sWorkspace_output_basin = oPyflowline_in.sWorkspace_output + slash + sBasin
            Path(sWorkspace_output_basin).mkdir(parents=True, exist_ok=True)
            #in this case, the sFilename_flowline_filter is a list of files            
            pBasin = oPyflowline_in.aBasin[i]
            iFlag_dam  = pBasin.iFlag_dam
            iFlag_disconnected = pBasin.iFlag_disconnected
            dThreshold = pBasin.dThreshold_small_river
            sFilename_flowline_filter = pBasin.sFilename_flowline_filter
            aFlowline_basin, pSpatialRef_pcs = read_flowline_shapefile( sFilename_flowline_filter )                
            if iFlag_dam ==1:
                sFilename_dam = pBasin.sFilename_dam
                aData_dam = text_reader_string(sFilename_dam, iSkipline_in =1,cDelimiter_in=',' )
                sFilename_flowline_topo = pBasin.sFilename_flowline_topo
                aData_flowline_topo = text_reader_string(sFilename_flowline_topo, iSkipline_in =1,cDelimiter_in=',' )

                aFromFlowline = aData_flowline_topo[:,1].astype(int).ravel()
                aToFlowline = aData_flowline_topo[:,2].astype(int).ravel()

                sFilename_flowline_raw = pBasin.sFilename_flowline_raw
                aNHDPlusID_filter = read_nhdplus_flowline_shapefile_attribute(sFilename_flowline_filter)
                aNHDPlusID_raw = read_nhdplus_flowline_shapefile_attribute(sFilename_flowline_raw)
                ndam = len(aData_dam)
                aNHDPlusID_dams_headwater = list()
                aNHDPlusID_dams_nonheadwater = list()
                for j in range(0, ndam):
                    #print(j)
                    dLon = float(aData_dam[j][1])
                    dLat = float(aData_dam[j][0])
                    sDam = aData_dam[j][4]            
                    lNHDPlusID = int(aData_dam[j][5])
                    aNHDPlusID_dams_headwater.append(lNHDPlusID)
                    
                    if lNHDPlusID in aNHDPlusID_filter:
                        #remove by id
                        for k in range(len(aFlowline_basin)):
                            if aFlowline_basin[k].lNHDPlusID == lNHDPlusID:
                                aFlowline_basin.pop(k)
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

                aFlowline_basin = aFlowline_basin + aFlowline_dams_headwater + aFlowline_dams_nonheadwater
            else:
                pass

            if iFlag_disconnected == 1:
                #need a better way to include this capability    
                #aThreshold = np.full(2, 300.0, dtype=float)
                #aFlowline_basin = connect_disconnect_flowline(aFlowline_basin, aVertex, aThreshold)
                #sFilename_out = 'flowline_connect.json'
                #sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)    
                #export_flowline_to_shapefile(iFlag_projected, aFlowline_basin,pSpatialRef_gcs, sFilename_out)
                pass
            else:
                pass

            
            sFilename_out = 'flowline_before_intersect_' + sBasin + '.shp'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)

            iFlag_projected = 0
            export_flowline_to_shapefile(iFlag_projected, aFlowline_basin, pSpatialRef_gcs, sFilename_out)


            aVertex = find_flowline_vertex(aFlowline_basin)
            sFilename_out = 'flowline_vertex_without_confluence_before_intersect_' + sBasin + '.shp'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_vertex_to_shapefile(iFlag_projected, aVertex,pSpatialRef_gcs, sFilename_out)

            aFlowline_basin = split_flowline(aFlowline_basin, aVertex)
            sFilename_out = 'flowline_split_by_point_before_intersect_' + sBasin + '.shp'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_flowline_to_shapefile(iFlag_projected, aFlowline_basin,pSpatialRef_gcs, sFilename_out)

            #ues location to find outlet

            point= dict()   
            point['lon'] = pBasin.dLon_outlet
            point['lat'] = pBasin.dLat_outlet
            pVertex_outlet=pyvertex(point)

            aFlowline_basin = correct_flowline_direction(aFlowline_basin,  pVertex_outlet )

            pVertex_outlet = aFlowline_basin[0].pVertex_end

            sFilename_out = 'flowline_direction_before_intersect_' + sBasin + '.shp'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_flowline_to_shapefile(iFlag_projected, aFlowline_basin, pSpatialRef_gcs, sFilename_out)

            #step 4: remove loops

            aFlowline_basin = remove_flowline_loop(aFlowline_basin)    
            sFilename_out = 'flowline_loop_before_intersect_' + sBasin + '.shp'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_flowline_to_shapefile(iFlag_projected, aFlowline_basin,pSpatialRef_gcs, sFilename_out)

            #using loop to remove small river, here we use 5 steps

            for i in range(3):
                sStep = "{:02d}".format(i+1)
                aFlowline_basin = remove_small_river(aFlowline_basin, dThreshold)
                sFilename_out = 'flowline_large_'+ sStep +'_before_intersect_' + sBasin + '.shp'
                sFilename_out =os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_shapefile(iFlag_projected, aFlowline_basin, pSpatialRef_gcs, sFilename_out)


                aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity = find_flowline_confluence(aFlowline_basin,  pVertex_outlet)
                sFilename_out = 'flowline_vertex_with_confluence_'+ sStep +'_before_intersect_' + sBasin + '.shp'
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_vertex_to_shapefile(iFlag_projected, aVertex, pSpatialRef_gcs, sFilename_out, aAttribute_data=aConnectivity)

                aFlowline_basin = merge_flowline( aFlowline_basin,aVertex, pVertex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence  )  
                sFilename_out = 'flowline_merge_'+ sStep +'_before_intersect_' + sBasin + '.shp'
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_shapefile(iFlag_projected, aFlowline_basin, pSpatialRef_gcs, sFilename_out)

                if len(aFlowline_basin) ==1:
                    break


            #build segment index
            aFlowline_basin, aStream_segment = define_stream_segment_index(aFlowline_basin)
            sFilename_out = pBasin.sFilename_flowline_segment_index_before_intersect
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_flowline_to_shapefile(iFlag_projected, \
                aFlowline_basin, pSpatialRef_gcs, sFilename_out, \
                aAttribute_data=[aStream_segment], aAttribute_field=['iseg'], aAttribute_dtype=['int'])

            #build stream order 
            aFlowline_basin, aStream_order = define_stream_order(aFlowline_basin)
            sFilename_out = pBasin.sFilename_flowline_segment_order_before_intersect
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_flowline_to_shapefile(iFlag_projected, \
                aFlowline_basin, pSpatialRef_gcs, sFilename_out, \
                aAttribute_data=[aStream_segment, aStream_order], aAttribute_field=['iseg','iord'], aAttribute_dtype=['int','int'])

            
            #add this basin into all flowline

            aFlowline = aFlowline + aFlowline_basin
    else:
        return

    print('Finished')