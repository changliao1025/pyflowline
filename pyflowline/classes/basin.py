import os
from pathlib import Path
from abc import ABCMeta, abstractmethod
import json
from json import JSONEncoder
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
from matplotlib import cm
import cartopy.crs as ccrs

from osgeo import ogr, osr, gdal
from shapely.wkt import loads

from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.edge import pyedge
from pyflowline.classes.flowline import pyflowline
from pyflowline.classes.confluence import pyconfluence
from pyflowline.formats.read_flowline import read_flowline_geojson
from pyflowline.formats.read_nhdplus_flowline_shapefile import read_nhdplus_flowline_shapefile_attribute
from pyflowline.formats.read_nhdplus_flowline_shapefile import extract_nhdplus_flowline_shapefile_by_attribute
from pyflowline.formats.read_nhdplus_flowline_shapefile import track_nhdplus_flowline
from pyflowline.formats.convert_shapefile_to_json import convert_shapefile_to_json
from pyflowline.formats.export_flowline import export_flowline_to_json
from pyflowline.formats.export_vertex import export_vertex_to_json

from pyflowline.algorithms.auxiliary.text_reader_string import text_reader_string
from pyflowline.algorithms.auxiliary.find_index_in_list import find_vertex_in_list
from pyflowline.algorithms.auxiliary.calculate_area_of_difference import calculate_area_of_difference_simplified


from pyflowline.algorithms.split.find_flowline_vertex import find_flowline_vertex
from pyflowline.algorithms.split.find_flowline_confluence import find_flowline_confluence
from pyflowline.algorithms.split.split_flowline import split_flowline
from pyflowline.algorithms.split.split_flowline_to_edge import split_flowline_to_edge

from pyflowline.algorithms.merge.merge_flowline import merge_flowline

from pyflowline.algorithms.direction.correct_flowline_direction import correct_flowline_direction

from pyflowline.algorithms.loop.remove_flowline_loop import remove_flowline_loop

from pyflowline.algorithms.simplification.remove_small_river import remove_small_river
from pyflowline.algorithms.simplification.remove_returning_flowline import remove_returning_flowline
from pyflowline.algorithms.simplification.remove_duplicate_flowline import remove_duplicate_flowline

from pyflowline.algorithms.index.define_stream_order import define_stream_order
from pyflowline.algorithms.index.define_stream_segment_index import define_stream_segment_index

from pyflowline.algorithms.intersect.intersect_flowline_with_mesh import intersect_flowline_with_mesh
from pyflowline.algorithms.intersect.intersect_flowline_with_flowline import intersect_flowline_with_flowline

#desired_proj = ccrs.Orthographic(central_longitude=-75, central_latitude=42, globe=None)
desired_proj = ccrs.PlateCarree()

class BasinClassEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.float32):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, list):
            pass  
        if isinstance(obj, pyvertex):
            return json.loads(obj.tojson()) #lVertexID
        if isinstance(obj, pyedge):
            return obj.lEdgeID        
        if isinstance(obj, pyflowline):
            return obj.lFlowlineID
        if isinstance(obj, pyconfluence):
            return obj.dAngle_upstream
        
        #if isinstance(obj, pybasin):
        #    return json.loads(obj.tojson())    
            
        return JSONEncoder.default(self, obj)



class pybasin(object):
    lBasinID =1 
    sBasinID=''
    lCellID_outlet=-1
    iFlag_debug = 0
    iFlag_disconnected =0
    iFlag_dam=0
    dLongitude_outlet_degree = -9999.
    dLatitude_outlet_degree = -9999.
    dAccumulation_threshold= 100000.0
    dThreshold_small_river = 10000
    dLength_flowline_before_simplification = 0.0
    dLength_flowline_after_simplification = 0.0
    dLength_flowline_conceptual = 0.0

    dArea_of_difference=0.0
    dDistance_displace = 0.0
    sWorkspace_output_basin=''
    sFilename_flowline_raw=''    
    sFilename_flowline_filter=''
    sFilename_flowline_filter_json=''
    sFilename_dam=''
    sFilename_flowline_topo=''
    #before intersect
    sFilename_flowline_segment_order_before_intersect=''
    sFilename_flowline_segment_index_before_intersect=''
    sFilename_flowline_final=''
    sFilename_flowline_edge=''
    sFilename_basin_info=''
    sFilename_flowline_simplified_info=''
    sFilename_flowline_conceptual_info=''
    #sFilename_basin_configuration=''
    aFlowline_basin_simplified=None
    aFlowline_basin_conceptual=None
    #aVertex_confluence=None
    pVertex_outlet=None
    aConfluence_basin_simplified= None
    aConfluence_basin_conceptual= None
    
    def __init__(self, aParameter):

        if 'lBasinID' in aParameter:            
            self.lBasinID             = int(aParameter['lBasinID'])
        else:
            self.lBasinID   = 1
        
        
        if 'lCellID_outlet' in aParameter:            
            self.lCellID_outlet             = int(aParameter['lCellID_outlet'])
        else:
            self.lCellID_outlet   = -1

        if 'iFlag_disconnected' in aParameter:            
            self.iFlag_disconnected             = int(aParameter['iFlag_disconnected'])
        else:
            self.iFlag_disconnected   = 0
        
        if 'iFlag_dam' in aParameter:            
            self.iFlag_dam             = int(aParameter['iFlag_dam'])
        else:
            self.iFlag_dam   = 0
        
        if 'dLongitude_outlet_degree' in aParameter:            
            self.dLongitude_outlet_degree             = float(aParameter['dLongitude_outlet_degree'])
        else:
            self.dLongitude_outlet_degree   = -9999.
        
        if 'dLatitude_outlet_degree' in aParameter:            
            self.dLatitude_outlet_degree             = float(aParameter['dLatitude_outlet_degree'])
        else:
            self.dLatitude_outlet_degree   = -9999.
        
        if 'dThreshold_small_river' in aParameter:            
            self.dThreshold_small_river             = float(aParameter['dThreshold_small_river'])
        else:
            self.dThreshold_small_river   =10000.0

        if 'dAccumulation_threshold' in aParameter:            
            self.dAccumulation_threshold             = float(aParameter['dAccumulation_threshold'])
        else:
            self.dAccumulation_threshold = 100000.0
        


        if 'sFilename_flowline_raw' in aParameter:
            self.sFilename_flowline_raw = aParameter['sFilename_flowline_raw']
        else:
            self.sFilename_flowline_raw   =''
       
        if 'sFilename_flowline_filter' in aParameter:
            self.sFilename_flowline_filter = aParameter['sFilename_flowline_filter']
        else:
            self.sFilename_flowline_filter   = ''

        if 'sWorkspace_output_basin' in aParameter:
            self.sWorkspace_output_basin = aParameter['sWorkspace_output_basin']
        else:
            self.sWorkspace_output_basin   = '.'
            print('The basin output path is not specified!')
        Path(self.sWorkspace_output_basin).mkdir(parents=True, exist_ok=True)

        self.sFilename_flowline_filter_json = os.path.join(str(self.sWorkspace_output_basin ), "flowline_filter.json"  )

        if 'sFilename_dam' in aParameter:
            self.sFilename_dam = aParameter['sFilename_dam']
        else:
            self.sFilename_dam   = ''

        if 'sFilename_flowline_topo' in aParameter:
            self.sFilename_flowline_topo = aParameter['sFilename_flowline_topo']
        else:
            self.sFilename_flowline_topo   =''

        self.sBasinID = sBasinID = "{:03d}".format(self.lBasinID)

        self.sFilename_flowline_segment_index_before_intersect = 'flowline_segment_index_before_intersect.json'
        self.sFilename_flowline_segment_order_before_intersect = 'flowline_segment_order_before_intersect.json'
        self.sFilename_flowline_intersect  = 'flowline_intersect_mesh.json'
        self.sFilename_flowline_final = 'flowline_final.json'
        self.sFilename_flowline_edge = 'flowline_edge.json'
        self.sFilename_area_of_difference = 'area_of_difference.json'
        self.sFilename_basin_info = 'basin_info.json'
        self.sFilename_flowline_conceptual_info = 'flowline_conceptual_info.json'
        self.sFilename_flowline_simplified_info = 'flowline_simplified_info.json'
        return
        
    def flowline_simplification(self):

        
        sFilename_flowline_filter = self.sFilename_flowline_filter
        sFilename_flowline_filter_json = self.sFilename_flowline_filter_json
        aFlowline_basin_simplified, pSpatial_reference = read_flowline_geojson( sFilename_flowline_filter_json )   
        sWorkspace_output_basin = self.sWorkspace_output_basin
        if self.iFlag_dam ==1:
            sFilename_dam = self.sFilename_dam
            aData_dam = text_reader_string(sFilename_dam, iSkipline_in =1,cDelimiter_in=',' )
            sFilename_flowline_topo = self.sFilename_flowline_topo
            aData_flowline_topo = text_reader_string(sFilename_flowline_topo, iSkipline_in =1,cDelimiter_in=',' )
            aFromFlowline = aData_flowline_topo[:,1].astype(int).ravel()
            aToFlowline = aData_flowline_topo[:,2].astype(int).ravel()
            sFilename_flowline_raw = self.sFilename_flowline_raw
            aNHDPlusID_filter = read_nhdplus_flowline_shapefile_attribute(sFilename_flowline_filter)
            aNHDPlusID_raw = read_nhdplus_flowline_shapefile_attribute(sFilename_flowline_raw)
            ndam = len(aData_dam)
            aNHDPlusID_dams_headwater = list()
            aNHDPlusID_dams_nonheadwater = list()
            for j in range(0, ndam):
                dLon = float(aData_dam[j][1])
                dLat = float(aData_dam[j][0])
                sDam = aData_dam[j][4]            
                lNHDPlusID = int(aData_dam[j][5])
                aNHDPlusID_dams_headwater.append(lNHDPlusID)
                if lNHDPlusID in aNHDPlusID_filter:
                    #remove by id
                    for k in range(len(aFlowline_basin_simplified)):
                        if aFlowline_basin_simplified[k].lNHDPlusID == lNHDPlusID:
                            aFlowline_basin_simplified.pop(k)
                            break
                    pass
                else:                                
                    aNHDPlusID_dam_nonheadwater = track_nhdplus_flowline(aNHDPlusID_filter, aFromFlowline, aToFlowline, lNHDPlusID)
                    aNHDPlusID_filter = aNHDPlusID_filter + aNHDPlusID_dams_headwater+ aNHDPlusID_dam_nonheadwater  
                    aNHDPlusID_dams_nonheadwater = aNHDPlusID_dams_nonheadwater + aNHDPlusID_dam_nonheadwater
            aFlowline_dams_headwater = extract_nhdplus_flowline_shapefile_by_attribute(sFilename_flowline_raw, aNHDPlusID_dams_headwater )
            for i in range(len(aFlowline_dams_headwater)):
                aFlowline_dams_headwater[i].iFlag_dam = 1
            aFlowline_dams_nonheadwater = extract_nhdplus_flowline_shapefile_by_attribute(sFilename_flowline_raw, aNHDPlusID_dams_nonheadwater )
            aFlowline_basin_simplified = aFlowline_basin_simplified + aFlowline_dams_headwater + aFlowline_dams_nonheadwater
        else:
            pass
        if self.iFlag_disconnected == 1:                
            #aThreshold = np.full(2, 300.0, dtype=float)
            #aFlowline_basin_simplified = connect_disconnect_flowline(aFlowline_basin_simplified, aVertex, aThreshold)
            #sFilename_out = 'flowline_connect.json'
            #sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)    
            #export_flowline_to_json(iFlag_projected, aFlowline_basin_simplified,pSpatial_reference_gcs, sFilename_out)
            pass
        else:
            pass
            
        if self.iFlag_debug ==1:
            sFilename_out = 'flowline_before_intersect.json'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)            
            self.export_flowline(aFlowline_basin_simplified, sFilename_out)
        #calculate length
        self.dLength_flowline_before_simplification = self.calculate_flowline_length(aFlowline_basin_simplified)
        
        aVertex = find_flowline_vertex(aFlowline_basin_simplified)
        if self.iFlag_debug ==1:
            sFilename_out = 'flowline_vertex_without_confluence_before_intersect.json'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_vertex_to_json( aVertex, sFilename_out)
        aFlowline_basin_simplified = split_flowline(aFlowline_basin_simplified, aVertex)
        if self.iFlag_debug ==1:
            sFilename_out = 'flowline_split_by_point_before_intersect.json'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_flowline_to_json(aFlowline_basin_simplified, sFilename_out)
        #ues location to find outlet
        point= dict()   
        point['dLongitude_degree'] = self.dLongitude_outlet_degree
        point['dLatitude_degree'] = self.dLatitude_outlet_degree
        pVertex_outlet=pyvertex(point)
        aFlowline_basin_simplified = correct_flowline_direction(aFlowline_basin_simplified,  pVertex_outlet )
        pVertex_outlet = aFlowline_basin_simplified[0].pVertex_end
        self.pVertex_outlet = pVertex_outlet
        if self.iFlag_debug ==1:
            sFilename_out = 'flowline_direction_before_intersect.json'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_flowline_to_json( aFlowline_basin_simplified,  sFilename_out)
        #step 4: remove loops
        aFlowline_basin_simplified = remove_flowline_loop(aFlowline_basin_simplified)    
        if self.iFlag_debug ==1:
            sFilename_out = 'flowline_loop_before_intersect.json'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_flowline_to_json( aFlowline_basin_simplified, sFilename_out)
        #using loop to remove small river, here we use 5 steps
        for i in range(3):
            sStep = "{:02d}".format(i+1)
            aFlowline_basin_simplified = remove_small_river(aFlowline_basin_simplified, self.dThreshold_small_river)
            if self.iFlag_debug ==1:
                sFilename_out = 'flowline_large_'+ sStep +'_before_intersect.json'
                sFilename_out =os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_json( aFlowline_basin_simplified,  sFilename_out)
            aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity = find_flowline_confluence(aFlowline_basin_simplified,  pVertex_outlet)
            if self.iFlag_debug ==1:
                sFilename_out = 'flowline_vertex_with_confluence_'+ sStep +'_before_intersect.json'
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_vertex_to_json( aVertex,  sFilename_out, aAttribute_data=aConnectivity)
            aFlowline_basin_simplified = merge_flowline( aFlowline_basin_simplified,aVertex, pVertex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence  )  
            if self.iFlag_debug ==1:
                sFilename_out = 'flowline_merge_'+ sStep +'_before_intersect.json'
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_json( aFlowline_basin_simplified,  sFilename_out)
            if len(aFlowline_basin_simplified) == 1:
                break
        
        #the final vertex info
        aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity = find_flowline_confluence(aFlowline_basin_simplified,  pVertex_outlet)
        sFilename_out = 'flowline_vertex_with_confluence_before_intersect.json'
        sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
        export_vertex_to_json( aVertex,  sFilename_out, aAttribute_data=aConnectivity)
        
        #self.dLength_flowline_after_simplification = self.calculate_flowline_length(aFlowline_basin_simplified)
        aVertex = np.array(aVertex)
        aIndex_confluence = np.array(aIndex_confluence)
        if aIndex_confluence.size > 0:        
            aVertex_confluence = aVertex[aIndex_confluence]
            self.aConfluence_basin_simplified = self.build_confluence(aFlowline_basin_simplified, aVertex_confluence) 
        

        #build segment index
        aFlowline_basin_simplified, aStream_segment = define_stream_segment_index(aFlowline_basin_simplified)
        if self.iFlag_debug ==1:
            sFilename_out = self.sFilename_flowline_segment_index_before_intersect
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_flowline_to_json(  aFlowline_basin_simplified, sFilename_out, \
                aAttribute_data=[aStream_segment], aAttribute_field=['iseg'], aAttribute_dtype=['int'])
        #build stream order 
        aFlowline_basin_simplified, aStream_order = define_stream_order(aFlowline_basin_simplified)
        sFilename_out = self.sFilename_flowline_segment_order_before_intersect
        sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
        export_flowline_to_json(  aFlowline_basin_simplified, sFilename_out, \
                aAttribute_data=[aStream_segment, aStream_order], aAttribute_field=['iseg','iord'], aAttribute_dtype=['int','int'])
        

        self.aFlowline_basin_simplified= aFlowline_basin_simplified
        return aFlowline_basin_simplified

    def reconstruct_topological_relationship(self, iMesh_type, sFilename_mesh):
        
        sWorkspace_output_basin = self.sWorkspace_output_basin
        sFilename_flowline = self.sFilename_flowline_segment_order_before_intersect
        sFilename_flowline_in = os.path.join(sWorkspace_output_basin, sFilename_flowline)
        aFlowline_basin_simplified, pSpatial_reference = read_flowline_geojson( sFilename_flowline_in )   
                
        sFilename_flowline_intersect = self.sFilename_flowline_intersect
        sFilename_flowline_intersect_out = os.path.join(sWorkspace_output_basin, sFilename_flowline_intersect)
        aCell, aCell_intersect_basin, aFlowline_intersect_all = intersect_flowline_with_mesh(iMesh_type, sFilename_mesh, \
            sFilename_flowline_in, sFilename_flowline_intersect_out)
        sFilename_flowline_filter_json = self.sFilename_flowline_filter
        
        point= dict()
        point['dLongitude_degree'] = self.dLongitude_outlet_degree
        point['dLatitude_degree'] = self.dLatitude_outlet_degree
        pVertex_outlet_initial=pyvertex(point)

        #from this point, aFlowline_basin is conceptual

        #segment based
        aFlowline_basin_conceptual, aFlowline_no_parallel, lCellID_outlet, pVertex_outlet \
            = remove_returning_flowline(iMesh_type, aCell_intersect_basin, pVertex_outlet_initial)
        if self.iFlag_debug ==1:
            sFilename_out = 'flowline_simplified_after_intersect.json'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)  
            export_flowline_to_json(aFlowline_basin_conceptual,  sFilename_out)

        #edge based
        aFlowline_basin_conceptual, aEdge = split_flowline_to_edge(aFlowline_basin_conceptual)
        aFlowline_basin_conceptual = remove_duplicate_flowline(aFlowline_basin_conceptual)
        aFlowline_basin_conceptual = correct_flowline_direction(aFlowline_basin_conceptual,  pVertex_outlet )
        aFlowline_basin_conceptual = remove_flowline_loop(  aFlowline_basin_conceptual )  
  
        aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity\
            = find_flowline_confluence(aFlowline_basin_conceptual,  pVertex_outlet)
        if self.iFlag_debug ==1:
            sFilename_out = 'flowline_vertex_with_confluence_after_intersect.json'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_vertex_to_json( aVertex,  sFilename_out, aAttribute_data=aConnectivity)
        
        #segment based
        aFlowline_basin_conceptual = merge_flowline( aFlowline_basin_conceptual,aVertex, pVertex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence  )                          
        aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity\
            = find_flowline_confluence(aFlowline_basin_conceptual,  pVertex_outlet)        
        aFlowline_basin_conceptual, aStream_segment = define_stream_segment_index(aFlowline_basin_conceptual)
        aFlowline_basin_conceptual, aStream_order = define_stream_order(aFlowline_basin_conceptual)

        #save confluence
        aVertex = np.array(aVertex)
        aIndex_confluence = np.array(aIndex_confluence)
        if aIndex_confluence.size > 0:        
            aVertex_confluence = aVertex[aIndex_confluence] 
            self.aConfluence_basin_conceptual = self.build_confluence(aFlowline_basin_conceptual, aVertex_confluence) 
          

        #edge based
        aFlowline_basin_edge, aEdge = split_flowline_to_edge(aFlowline_basin_conceptual)
        sFilename_out = self.sFilename_flowline_edge
        sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
        export_flowline_to_json(  aFlowline_basin_edge, sFilename_out)

        sFilename_out = self.sFilename_flowline_final
        sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
        export_flowline_to_json(  aFlowline_basin_conceptual, sFilename_out, \
            aAttribute_data=[aStream_segment, aStream_order], aAttribute_field=['iseg','iord'], aAttribute_dtype=['int','int'])

        self.aFlowline_basin_conceptual = aFlowline_basin_conceptual     
        
        self.lCellID_outlet = lCellID_outlet
        self.dLongitude_outlet_degree = pVertex_outlet.dLongitude_degree
        self.dLatitude_outlet_degree = pVertex_outlet.dLatitude_degree
        

        return aCell_intersect_basin

    def build_confluence(self, aFlowline_basin_in, aVertex_confluence_in):    
        #this can only be calculated for confluence
        aConfluence_basin=list()
        for pVertex in aVertex_confluence_in:
            #aEdge=list()
            aFlowline_upstream =list()
            for pFlowline in aFlowline_basin_in:
                pVertex_start = pFlowline.pVertex_start
                pVertex_end = pFlowline.pVertex_end
                if pVertex_end == pVertex:
                    #pEdge = pFlowline.aEdge[pFlowline.nEdge-1]
                    aFlowline_upstream.append(pFlowline)
                    pass
                if pVertex_start == pVertex:
                    pFlowline_downstream=pFlowline

            pConfluence = pyconfluence(pVertex, aFlowline_upstream, pFlowline_downstream)
            aConfluence_basin.append(pConfluence)   
        return aConfluence_basin

    def analyze(self):        
        self.dLength_flowline_after_simplification = self.calculate_flowline_length(self.aFlowline_basin_simplified)
        self.dLength_flowline_conceptual = self.calculate_flowline_length(self.aFlowline_basin_conceptual)
        self.calculate_river_sinuosity()
        self.calculate_confluence_branching_angle()
        return    

    def evaluate(self, iMesh_type, sMesh_type):        
        self.evaluate_area_of_difference(iMesh_type, sMesh_type)
        return

    def evaluate_area_of_difference(self, iMesh_type, sMesh_type):

        sFilename_simplified =  self.sFilename_flowline_segment_order_before_intersect
        sFilename_simplified= os.path.join(self.sWorkspace_output_basin, sFilename_simplified)

        sFilename_flowline_edge = self.sFilename_flowline_edge
        sFilename_flowline_edge= os.path.join(self.sWorkspace_output_basin, sFilename_flowline_edge)

        
        #intersect first
        sFilename_output= os.path.join(self.sWorkspace_output_basin, 'flowline_intersect_flowline.json')
        aVertex_intersect = intersect_flowline_with_flowline(sFilename_simplified, sFilename_flowline_edge, sFilename_output)

        #get confluence simple
        point= dict()   
        point['dLongitude_degree'] = self.dLongitude_outlet_degree
        point['dLatitude_degree'] = self.dLatitude_outlet_degree
        pVertex_outlet=pyvertex(point)

        aFlowline_simplified,pSpatial_reference = read_flowline_geojson( sFilename_simplified )   
        aVertex_simplified, lIndex_outlet_simplified, \
            aIndex_headwater_simplified, aIndex_middle, \
                aIndex_confluence_simplified, aConnectivity\
                = find_flowline_confluence(aFlowline_simplified,  pVertex_outlet)

        
        aFlowline_conceptual,pSpatial_reference = read_flowline_geojson( sFilename_flowline_edge ) 
        aVertex_conceptual, lIndex_outlet_conceptual, \
            aIndex_headwater_conceptual, aIndex_middle_conceptual, \
            aIndex_confluence_conceptual,  aConnectivity \
                = find_flowline_confluence(aFlowline_conceptual,  pVertex_outlet)

        #merge intersect with confluence
        a=np.array(aIndex_confluence_simplified)
        b=np.array(aIndex_confluence_conceptual)
        c=np.array(aIndex_middle_conceptual)
        d = list(np.array(aVertex_simplified)[a] )
        e = list(np.array(aVertex_conceptual)[b] )
        f = list(np.array(aVertex_conceptual)[c] )

        g= aVertex_intersect  + d + e + f
        
        aVertex_all =list()
        for i in g:
            iFlag_exist, lIndex = find_vertex_in_list( aVertex_all,  i)
            if iFlag_exist ==1:
                pass
            else:
                aVertex_all.append(i)

        h = aVertex_intersect  + d         
        aVertex_all_simplified =list()
        for i in h:
            iFlag_exist, lIndex = find_vertex_in_list( aVertex_all_simplified,  i)
            if iFlag_exist ==1:
                pass
            else:
                aVertex_all_simplified.append(i)

        j = aVertex_intersect  + e + f       
        aVertex_all_conceptual =list()
        for i in j:
            iFlag_exist, lIndex = find_vertex_in_list( aVertex_all_conceptual,  i)
            if iFlag_exist ==1:
                pass
            else:
                aVertex_all_conceptual.append(i)
        
        #export 
        self.iFlag_debug =1
        if self.iFlag_debug ==1:
            sFilename_output= os.path.join(self.sWorkspace_output_basin, 'vertex_split_all.json')
            export_vertex_to_json( aVertex_all, sFilename_output)
        
        #split 
        aFlowline_simplified_split = split_flowline(aFlowline_simplified, aVertex_all_simplified,iFlag_intersect =1)
        self.iFlag_debug =1
        if self.iFlag_debug ==1:
            sFilename_out = 'flowline_split_simplified.json'
            sFilename_out = os.path.join(self.sWorkspace_output_basin, sFilename_out)
            export_flowline_to_json(aFlowline_simplified_split, sFilename_out)     
        
        aFlowline_conceptual_split = split_flowline(aFlowline_conceptual, aVertex_all_conceptual,\
            iFlag_intersect =1, iFlag_use_id=1)
        self.iFlag_debug =1
        if self.iFlag_debug ==1:
            sFilename_out = 'flowline_split_conceptual.json'
            sFilename_out = os.path.join(self.sWorkspace_output_basin, sFilename_out)  
            export_flowline_to_json(aFlowline_conceptual_split, sFilename_out)
            #aFlowline_conceptual_split, dummy = read_flowline_geojson(sFilename_out)

        aFlowline_all = aFlowline_simplified_split + aFlowline_conceptual_split

        sFilename_area_of_difference = self.sFilename_area_of_difference
        sFilename_output = os.path.join(self.sWorkspace_output_basin, sFilename_area_of_difference)
        #remove headwater not needed here

        aPolygon_out, dArea = calculate_area_of_difference_simplified(aFlowline_all, aVertex_all, sFilename_output)
        print('Area of difference: ', dArea)
        self.dArea_of_difference = dArea
        self.dDistance_displace = dArea / self.dLength_flowline_after_simplification
        self.plot_area_of_difference(iMesh_type, sMesh_type)
        return

    def plot(self, sVariable_in=None):
        sWorkspace_output_basin = self.sWorkspace_output_basin
        if sVariable_in is not None:
            if sVariable_in == 'flowline_filter_json':
                sFilename_json = self.sFilename_flowline_filter_json
                sTitle = 'Original flowline'
            else:
                if sVariable_in == 'flowline_simplified':
                    sFilename_out = self.sFilename_flowline_segment_index_before_intersect
                    sFilename_json = os.path.join(sWorkspace_output_basin, sFilename_out)
                    sTitle = 'Simplified flowline'
                else:
                    sFilename_out = self.sFilename_flowline_final
                    sFilename_json = os.path.join(sWorkspace_output_basin, sFilename_out)
                    sTitle = 'Conceptual flowline'
                pass
        else:
            #default 
            sFilename_json = self.sFilename_flowline_filter_json
        
        #request = cimgt.OSM()
        fig = plt.figure( dpi=300)
        fig.set_figwidth( 4 )
        fig.set_figheight( 4 )
        ax = fig.add_axes([0.1, 0.15, 0.75, 0.8] , projection=desired_proj ) #request.crs
        
        pDriver = ogr.GetDriverByName('GeoJSON')
        pDataset = pDriver.Open(sFilename_json, gdal.GA_ReadOnly)
        pLayer = pDataset.GetLayer(0)
    
        pSrs = osr.SpatialReference()  
        pSrs.ImportFromEPSG(4326)    # WGS84 lat/lon
    
        lID = 0
        dLat_min = 90
        dLat_max = -90
        dLon_min = 180
        dLon_max = -180          
        

        #ax.add_image(request, 6)    # 5 = zoom level

        n_colors = pLayer.GetFeatureCount()
        
        colours = cm.rainbow(np.linspace(0, 1, n_colors))
        for pFeature in pLayer:
            pGeometry_in = pFeature.GetGeometryRef()
            sGeometry_type = pGeometry_in.GetGeometryName()
            if sGeometry_type =='LINESTRING':
                dummy0 = loads( pGeometry_in.ExportToWkt() )
                aCoords_gcs = dummy0.coords
                aCoords_gcs= np.array(aCoords_gcs)
                nvertex = len(aCoords_gcs)
                for i in range(nvertex):
                    dLon = aCoords_gcs[i][0]
                    dLat = aCoords_gcs[i][1]
                    if dLon > dLon_max:
                        dLon_max = dLon
                    
                    if dLon < dLon_min:
                        dLon_min = dLon
                    
                    if dLat > dLat_max:
                        dLat_max = dLat
    
                    if dLat < dLat_min:
                        dLat_min = dLat
    
                codes = np.full(nvertex, mpath.Path.LINETO, dtype=int )
                codes[0] = mpath.Path.MOVETO
                path = mpath.Path(aCoords_gcs, codes)            
                x, y = zip(*path.vertices)
                line, = ax.plot(x, y, color= colours[lID],linewidth=1)
                lID = lID + 1
                
    
        pDataset = pLayer = pFeature  = None    
        sDirname = os.path.dirname(sFilename_json)
        marginx  = (dLon_max - dLon_min) / 20
        marginy  = (dLat_max - dLat_min) / 20
        aExtent_in = [dLon_min - marginx , dLon_max + marginx , dLat_min - marginy , dLat_max + marginy]
        #aExtent_in = [-76.5,-76.2, 41.6,41.9]
        #aExtent_in = [-76.95,-76.75, 40.7,40.9]
        sFilename  = Path(sFilename_json).stem + '.png'
        #sFilename  = Path(sFilename_json).stem + '_meander.png'       
        #sFilename  = Path(sFilename_json).stem + '_loop.png'  
        ax.set_extent(aExtent_in)       
    
        ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.3, linestyle='--')
        ax.set_title( sTitle)       
        
        sFilename_out = os.path.join(sDirname, sFilename)
        plt.savefig(sFilename_out, bbox_inches='tight')
        #plt.show()
    
        return

    def plot_area_of_difference(self, iMesh_type, sMesh_type):
        #request = cimgt.OSM()
        sFilename_json = self.sFilename_area_of_difference
        sFilename_json = os.path.join(self.sWorkspace_output_basin, sFilename_json)
       
        fig = plt.figure( dpi=300)
        fig.set_figwidth( 4 )
        fig.set_figheight( 4 )
        ax = fig.add_axes([0.1, 0.15, 0.75, 0.8] , projection=desired_proj ) #request.crs
        
        pDriver = ogr.GetDriverByName('GeoJSON')
        pDataset = pDriver.Open(sFilename_json, gdal.GA_ReadOnly)
        pLayer = pDataset.GetLayer(0)
    
        pSrs = osr.SpatialReference()  
        pSrs.ImportFromEPSG(4326)    # WGS84 lat/lon
    
        lID = 0
        dLat_min = 90
        dLat_max = -90
        dLon_min = 180
        dLon_max = -180          
        

        #ax.add_image(request, 6)    # 5 = zoom level

        n_colors = pLayer.GetFeatureCount()
        
        colours = cm.rainbow(np.linspace(0, 1, n_colors))
        lID=0
        for pFeature in pLayer:
            pGeometry_in = pFeature.GetGeometryRef()
            sGeometry_type = pGeometry_in.GetGeometryName()
            if sGeometry_type =='POLYGON':
                dummy0 = loads( pGeometry_in.ExportToWkt() )
                aCoords_gcs = dummy0.exterior.coords
                aCoords_gcs= np.array(aCoords_gcs)
                nvertex = len(aCoords_gcs)
                for i in range(nvertex):
                    dLon = aCoords_gcs[i][0]
                    dLat = aCoords_gcs[i][1]
                    if dLon > dLon_max:
                        dLon_max = dLon
                    
                    if dLon < dLon_min:
                        dLon_min = dLon
                    
                    if dLat > dLat_max:
                        dLat_max = dLat
    
                    if dLat < dLat_min:
                        dLat_min = dLat
    
                polygon = mpatches.Polygon(aCoords_gcs[:,0:2], closed=True,  linewidth=0.25, \
                    alpha=0.8, edgecolor = 'red',facecolor='red', \
                        transform=ccrs.PlateCarree() )

                ax.add_patch(polygon)   
                lID = lID + 1
                
    
        pDataset = pLayer = pFeature  = None    
        sDirname = os.path.dirname(sFilename_json)
       
        marginx  = (dLon_max - dLon_min) / 20
        marginy  = (dLat_max - dLat_min) / 20
        aExtent_in = [dLon_min - marginx , dLon_max + marginx , dLat_min - marginy , dLat_max + marginy]
       
        
        ax.set_extent( aExtent_in )  
        
        #aExtent_in = [-76.5,-76.2, 41.6,41.9]
        #aExtent_in = [-76.95,-76.75, 40.7,40.9]
        sFilename  = Path(sFilename_json).stem + '.png'
        #sFilename  = Path(sFilename_json).stem + '_meander.png'       
        #sFilename  = Path(sFilename_json).stem + '_loop.png'  
              
        sTitle = 'Area of difference'
        ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.3, linestyle='--')
        ax.set_title( sTitle)  

        sText = 'Mesh type: ' + sMesh_type.title()
        ax.text(0.05, 0.95, sText, \
        verticalalignment='top', horizontalalignment='left',\
                transform=ax.transAxes, \
                color='black', fontsize=8)
     
        sText = 'Total area: ' + "{:4.1f}".format( int(self.dArea_of_difference/1.0E6) ) + ' km^2'
        ax.text(0.05, 0.90, sText, \
        verticalalignment='top', horizontalalignment='left',\
                transform=ax.transAxes, \
                color='blue', fontsize=8)
        
        sFilename_out = os.path.join(sDirname, sFilename)
        plt.savefig(sFilename_out, bbox_inches='tight')
        #plt.show()
        return
    
    def export(self):
        self.export_basin_info_to_json()
        self.export_flowline_info_to_json()
        self.tojson()    
        return

    def export_flowline(self, aFlowline_in, sFilename_json_in,iFlag_projected_in = None,  pSpatial_reference_in = None):
        export_flowline_to_json(aFlowline_in, sFilename_json_in,\
            iFlag_projected_in= iFlag_projected_in, \
            pSpatial_reference_in = pSpatial_reference_in)

    def export_basin_info_to_json(self):
        sFilename_json = self.sFilename_basin_info
        sFilename_json = os.path.join(str(Path(self.sWorkspace_output_basin)  ) , sFilename_json  )

        with open(sFilename_json, 'w', encoding='utf-8') as f:
            sJson = json.dumps(self.__dict__, \
            sort_keys=True, \
                indent = 4, \
                    ensure_ascii=True, \
                        cls=BasinClassEncoder)      
            f.write(sJson)    
            f.close()
        return

    def export_flowline_info_to_json(self):
        sFilename_json = self.sFilename_flowline_simplified_info
        sFilename_json = os.path.join(str(Path(self.sWorkspace_output_basin)  ) , sFilename_json  )
        with open(sFilename_json, 'w', encoding='utf-8') as f:
            sJson = json.dumps([json.loads(ob.tojson()) for ob in self.aFlowline_basin_simplified], indent = 4)        
            f.write(sJson)    
            f.close()

        sFilename_json = self.sFilename_flowline_conceptual_info
        sFilename_json = os.path.join(str(Path(self.sWorkspace_output_basin)  ) , sFilename_json  )
        with open(sFilename_json, 'w', encoding='utf-8') as f:
            sJson = json.dumps([json.loads(ob.tojson()) for ob in self.aFlowline_basin_conceptual], indent = 4)        
            f.write(sJson)    
            f.close()
        return

    def tojson(self):
        sJson = json.dumps(self.__dict__, \
            sort_keys=True, \
                indent = 4, \
                    ensure_ascii=True, \
                        cls=BasinClassEncoder)
        return sJson
    
    def convert_flowline_to_json(self):
        sFilename_raw = self.sFilename_flowline_filter            
        sFilename_out = self.sFilename_flowline_filter_json
        convert_shapefile_to_json(1, sFilename_raw, sFilename_out)
        
    def calculate_flowline_length(self, aFlowline_in):

        dLength = 0.0
        nflowline = len(aFlowline_in)
        for i in range(nflowline):
            pFlowline= aFlowline_in[i]
            pFlowline.calculate_length()
            dLength = dLength + pFlowline.dLength        
        return dLength

    def calculate_river_sinuosity(self):
        #the numner of segment
        for pFlowline in self.aFlowline_basin_simplified:
            pFlowline.calculate_flowline_sinuosity() 

        for pFlowline in self.aFlowline_basin_conceptual:
            pFlowline.calculate_flowline_sinuosity()    

        return

    def calculate_confluence_branching_angle(self):
        for pConfluence in self.aConfluence_basin_simplified:
            pConfluence.calculate_branching_angle()
        for pConfluence in self.aConfluence_basin_conceptual:
            pConfluence.calculate_branching_angle()    
        return
    