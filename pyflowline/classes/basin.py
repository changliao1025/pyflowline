import os
from pathlib import Path
from abc import ABCMeta, abstractmethod
import json
from json import JSONEncoder
import numpy as np
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
    dLength_flowline_filtered = 0.0
    dLength_flowline_simplified = 0.0
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
    sFilename_flowline_simplified=''
    sFilename_flowline_segment_index_before_intersect=''
    sFilename_flowline_conceptual=''
    sFilename_flowline_edge=''
    sFilename_basin_info=''
    sFilename_flowline_simplified_info=''
    sFilename_flowline_conceptual_info=''
    sFilename_confluence_simplified_info=''
    sFilename_confluence_conceptual_info=''    
    aFlowline_basin_filtered=None
    aFlowline_basin_simplified=None
    aFlowline_basin_conceptual=None    
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
            self.dThreshold_small_river   = 10000.0

        if 'dAccumulation_threshold' in aParameter:            
            self.dAccumulation_threshold             = float(aParameter['dAccumulation_threshold'])
        else:
            self.dAccumulation_threshold = 100000.0   

        if 'sFilename_flowline_raw' in aParameter:
            self.sFilename_flowline_raw = aParameter['sFilename_flowline_raw']
        else:
            self.sFilename_flowline_raw   = ''
       
        if 'sFilename_flowline_filter' in aParameter:
            self.sFilename_flowline_filter = aParameter['sFilename_flowline_filter']
        else:
            self.sFilename_flowline_filter   = ''

        if 'sWorkspace_output_basin' in aParameter:
            self.sWorkspace_output_basin = aParameter['sWorkspace_output_basin']
        else:
            self.sWorkspace_output_basin   = '.'
            
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
        self.sFilename_flowline_simplified = 'flowline_simplified.json'
        self.sFilename_flowline_intersect  = 'flowline_intersect_mesh.json'
        self.sFilename_flowline_conceptual = 'flowline_conceptual.json'
        self.sFilename_flowline_edge = 'flowline_edge.json'
        self.sFilename_area_of_difference = 'area_of_difference.json'
        self.sFilename_basin_info = 'basin_info.json'
        self.sFilename_flowline_conceptual_info = 'flowline_conceptual_info.json'
        self.sFilename_flowline_simplified_info = 'flowline_simplified_info.json'
        self.sFilename_confluence_conceptual_info = 'confluence_conceptual_info.json'
        self.sFilename_confluence_simplified_info = 'confluence_simplified_info.json'
        return
        
    def flowline_simplification(self):

        
        sFilename_flowline_filter = self.sFilename_flowline_filter
        sFilename_flowline_filter_json = self.sFilename_flowline_filter_json
        aFlowline_basin_filtered, pSpatial_reference = read_flowline_geojson( sFilename_flowline_filter_json )   
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
                    for k in range(len(aFlowline_basin_filtered)):
                        if aFlowline_basin_filtered[k].lNHDPlusID == lNHDPlusID:
                            aFlowline_basin_filtered.pop(k)
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
            aFlowline_basin_filtered = aFlowline_basin_filtered + aFlowline_dams_headwater + aFlowline_dams_nonheadwater
        else:
            pass
        if self.iFlag_disconnected == 1:
            #not used anymore                
            #aThreshold = np.full(2, 300.0, dtype=float)
            #aFlowline_basin_filtered = connect_disconnect_flowline(aFlowline_basin_filtered, aVertex, aThreshold)
            #sFilename_out = 'flowline_connect.json'
            #sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)    
            #export_flowline_to_json(iFlag_projected, aFlowline_basin_filtered,pSpatial_reference_gcs, sFilename_out)
            pass
        else:
            pass
            
        if self.iFlag_debug ==1:
            sFilename_out = 'flowline_before_intersect.json'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)            
            self.export_flowline(aFlowline_basin_filtered, sFilename_out)
        #calculate length
        self.aFlowline_basin_filtered = aFlowline_basin_filtered
        self.dLength_flowline_filtered = self.calculate_flowline_length(aFlowline_basin_filtered)

        #simplification started
        aVertex = find_flowline_vertex(aFlowline_basin_filtered)
        if self.iFlag_debug ==1:
            sFilename_out = 'flowline_vertex_without_confluence_before_intersect.json'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_vertex_to_json( aVertex, sFilename_out)
        aFlowline_basin_simplified = split_flowline(aFlowline_basin_filtered, aVertex)
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
        sFilename_out = 'vertex_simplified.json'
        sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
        export_vertex_to_json( aVertex,  sFilename_out, aAttribute_data=aConnectivity)
        
        #self.dLength_flowline_simplified = self.calculate_flowline_length(aFlowline_basin_simplified)
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
        sFilename_out = self.sFilename_flowline_simplified
        sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
        export_flowline_to_json(  aFlowline_basin_simplified, sFilename_out, \
                aAttribute_data=[aStream_segment, aStream_order], aAttribute_field=['iseg','iord'], aAttribute_dtype=['int','int'])
        

        self.aFlowline_basin_simplified= aFlowline_basin_simplified
        return aFlowline_basin_simplified

    def reconstruct_topological_relationship(self, iMesh_type, sFilename_mesh):
        
        sWorkspace_output_basin = self.sWorkspace_output_basin
        sFilename_flowline = self.sFilename_flowline_simplified
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
        aFlowline_basin_conceptual, lCellID_outlet, pVertex_outlet \
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

        sFilename_out = self.sFilename_flowline_conceptual
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
            aFlowline_upstream =list()
            for pFlowline in aFlowline_basin_in:
                pVertex_start = pFlowline.pVertex_start
                pVertex_end = pFlowline.pVertex_end
                if pVertex_end == pVertex:                 
                    aFlowline_upstream.append(pFlowline)
                    pass
                if pVertex_start == pVertex:
                    pFlowline_downstream=pFlowline

            pConfluence = pyconfluence(pVertex, aFlowline_upstream, pFlowline_downstream)
            aConfluence_basin.append(pConfluence)   
        return aConfluence_basin

    def analyze(self):      
        if self.aFlowline_basin_filtered is None:
            sFilename_flowline_filter = self.sFilename_flowline_filter
            sFilename_flowline_filter_json = self.sFilename_flowline_filter_json
            self.aFlowline_basin_filtered, pSpatial_reference = read_flowline_geojson( sFilename_flowline_filter_json )   
            self.dLength_flowline_filtered = self.calculate_flowline_length(self.aFlowline_basin_filtered)
        
        point= dict()
        point['dLongitude_degree'] = self.dLongitude_outlet_degree
        point['dLatitude_degree'] = self.dLatitude_outlet_degree
        pVertex_outlet_initial=pyvertex(point)
        if self.aFlowline_basin_simplified is None:
            sFilename_flowline = self.sFilename_flowline_simplified
            sFilename_flowline_in = os.path.join(self.sWorkspace_output_basin, sFilename_flowline)
            aFlowline_simplified,pSpatial_reference = read_flowline_geojson( sFilename_flowline_in )   
            read_flowline_geojson
            self.aFlowline_basin_simplified = aFlowline_simplified
            aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity\
            = find_flowline_confluence(self.aFlowline_basin_simplified,  pVertex_outlet_initial)  
            aVertex = np.array(aVertex)
            aIndex_confluence = np.array(aIndex_confluence)
            if aIndex_confluence.size > 0:        
                aVertex_confluence = aVertex[aIndex_confluence] 
                self.aConfluence_basin_simplified = self.build_confluence(self.aFlowline_basin_simplified, aVertex_confluence)

        self.dLength_flowline_simplified = self.calculate_flowline_length(self.aFlowline_basin_simplified)

        if self.aFlowline_basin_conceptual is None:
            sFilename_flowline = self.sFilename_flowline_conceptual
            sFilename_flowline_in = os.path.join(self.sWorkspace_output_basin, sFilename_flowline)
            aFlowline_conceptual, pSpatial_reference = read_flowline_geojson( sFilename_flowline_in )   
            self.aFlowline_basin_conceptual = aFlowline_conceptual
            
            aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity\
            = find_flowline_confluence(self.aFlowline_basin_conceptual,  pVertex_outlet_initial)  
            aVertex = np.array(aVertex)
            aIndex_confluence = np.array(aIndex_confluence)
            if aIndex_confluence.size > 0:        
                aVertex_confluence = aVertex[aIndex_confluence] 
                self.aConfluence_basin_conceptual = self.build_confluence(self.aFlowline_basin_conceptual, aVertex_confluence)

        self.dLength_flowline_conceptual = self.calculate_flowline_length(self.aFlowline_basin_conceptual)
        self.calculate_river_sinuosity()
        self.calculate_confluence_branching_angle()
        return    
    
    def export(self):
        self.export_basin_info_to_json()
        self.export_flowline_info_to_json()
        self.export_confluence_info_to_json()        
        return

    def export_flowline(self, aFlowline_in, sFilename_json_in,iFlag_projected_in = None,  pSpatial_reference_in = None):
        export_flowline_to_json(aFlowline_in, sFilename_json_in,\
            iFlag_projected_in= iFlag_projected_in, \
            pSpatial_reference_in = pSpatial_reference_in)

    def export_basin_info_to_json(self):
        sFilename_json = self.sFilename_basin_info
        sFilename_json = os.path.join(str(Path(self.sWorkspace_output_basin)  ) , sFilename_json  )

        aSkip = ['aFlowline_basin_filtered', \
                'aFlowline_basin_simplified','aFlowline_basin_conceptual','aConfluence_basin_simplified',
                'aConfluence_basin_conceptual']
        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)
        with open(sFilename_json, 'w', encoding='utf-8') as f:
            sJson = json.dumps(obj, default=lambda o: o.__dict__,\
            sort_keys=True, \
                indent = 4, \
                    ensure_ascii=True, \
                        cls=BasinClassEncoder)      
            f.write(sJson)    
            f.close()
        return

    def export_flowline_info_to_json(self):
        iFlag_export_simplified=0
        if iFlag_export_simplified==1:
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

    def export_confluence_info_to_json(self):
        iFlag_export_confluence =0
        if iFlag_export_confluence==1:
            sFilename_json = self.sFilename_confluence_simplified_info
            sFilename_json = os.path.join(str(Path(self.sWorkspace_output_basin)  ) , sFilename_json  )
            with open(sFilename_json, 'w', encoding='utf-8') as f:
                sJson = json.dumps([json.loads(ob.tojson()) for ob in self.aConfluence_basin_simplified], indent = 4)        
                f.write(sJson)    
                f.close()

        sFilename_json = self.sFilename_confluence_conceptual_info
        sFilename_json = os.path.join(str(Path(self.sWorkspace_output_basin)  ) , sFilename_json  )
        
        with open(sFilename_json, 'w', encoding='utf-8') as f:
            sJson = json.dumps([json.loads(ob.tojson()) for ob in self.aConfluence_basin_conceptual], indent = 4)        
            f.write(sJson)    
            f.close()
        return   

    def tojson(self):
        aSkip = ['aFlowline_basin_filtered', \
                'aFlowline_basin_simplified','aFlowline_basin_conceptual','aConfluence_basin_simplified',
                'aConfluence_basin_conceptual']

        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)
 
    
        sJson = json.dumps(obj, \
            sort_keys=True, \
                indent = 4, \
                    ensure_ascii=True, \
                        cls=BasinClassEncoder)
        return sJson
    
    def export_config_to_json(self, sFilename_output_in = None):
        #single basin
        if sFilename_output_in is not None:
            sFilename_output = sFilename_output_in
        else:
            sFilename_output = os.path.join(self.sWorkspace_output_basin, 'configuration_basin.json' )

        aSkip = ['aFlowline_basin_filtered', \
                'aFlowline_basin_simplified','aFlowline_basin_conceptual','aConfluence_basin_simplified',
                'aConfluence_basin_conceptual']
        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)
        with open(sFilename_output, 'w', encoding='utf-8') as f:
            json.dump(obj, f,sort_keys=True, \
                ensure_ascii=False, \
                indent=4, \
                cls=BasinClassEncoder)
        return

    def convert_flowline_to_json(self):
        sFilename_raw = self.sFilename_flowline_filter            
        sFilename_out = self.sFilename_flowline_filter_json
        print('This is the filtered flowline:', sFilename_raw )
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
    