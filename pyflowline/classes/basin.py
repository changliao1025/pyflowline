import os
from abc import ABCMeta, abstractmethod
import json
from json import JSONEncoder
from pathlib import Path
import numpy as np
from osgeo import ogr, osr, gdal, gdalconst
from shapely.wkt import loads
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
from matplotlib import cm

import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt

from pyflowline.classes.vertex import pyvertex

from pyflowline.classes.edge import pyedge
from pyflowline.classes.cell import pycell
from pyflowline.classes.flowline import pyflowline
from pyflowline.algorithms.auxiliary.text_reader_string import text_reader_string
from pyflowline.formats.read_flowline import read_flowline_geojson

from pyflowline.formats.read_nhdplus_flowline_shapefile import read_nhdplus_flowline_shapefile_attribute
from pyflowline.formats.read_nhdplus_flowline_shapefile import extract_nhdplus_flowline_shapefile_by_attribute
from pyflowline.formats.read_nhdplus_flowline_shapefile import track_nhdplus_flowline
from pyflowline.formats.convert_shapefile_to_json import convert_shapefile_to_json
from pyflowline.formats.export_flowline import export_flowline_to_json
from pyflowline.algorithms.split.find_flowline_vertex import find_flowline_vertex
from pyflowline.formats.export_vertex import export_vertex_to_json
from pyflowline.algorithms.merge.merge_flowline import merge_flowline
from pyflowline.algorithms.split.split_flowline import split_flowline
from pyflowline.algorithms.split.find_flowline_confluence import find_flowline_confluence
from pyflowline.algorithms.split.find_flowline_vertex import find_flowline_vertex
from pyflowline.algorithms.direction.correct_flowline_direction import correct_flowline_direction
from pyflowline.algorithms.loop.remove_flowline_loop import remove_flowline_loop
from pyflowline.algorithms.simplification.remove_small_river import remove_small_river
from pyflowline.algorithms.index.define_stream_order import define_stream_order
from pyflowline.algorithms.index.define_stream_segment_index import define_stream_segment_index

from pyflowline.formats.read_mesh import read_mesh_json

from pyflowline.formats.export_vertex import export_vertex_to_json
from pyflowline.formats.export_flowline import export_flowline_to_json
from pyflowline.algorithms.intersect.intersect_flowline_with_mesh import intersect_flowline_with_mesh

from pyflowline.algorithms.simplification.remove_returning_flowline import remove_returning_flowline
from pyflowline.algorithms.simplification.remove_duplicate_flowline import remove_duplicate_flowline
from pyflowline.algorithms.simplification.remove_duplicate_edge import remove_duplicate_edge
from pyflowline.algorithms.direction.correct_flowline_direction import correct_flowline_direction
from pyflowline.algorithms.loop.remove_flowline_loop import remove_flowline_loop
from pyflowline.algorithms.split.find_flowline_vertex import find_flowline_vertex
from pyflowline.algorithms.split.find_flowline_confluence import find_flowline_confluence
from pyflowline.algorithms.split.split_flowline import split_flowline
from pyflowline.algorithms.split.split_flowline_to_edge import split_flowline_to_edge
from pyflowline.formats.export_vertex import export_vertex_to_shapefile
from pyflowline.algorithms.merge.merge_flowline import merge_flowline

from pyflowline.algorithms.index.define_stream_order import define_stream_order
from pyflowline.algorithms.index.define_stream_segment_index import define_stream_segment_index
from pyflowline.algorithms.area.calculate_area_of_difference import calculate_area_of_difference_simplified

from pyflowline.algorithms.intersect.intersect_flowline_with_flowline import intersect_flowline_with_flowline


desired_proj = ccrs.Orthographic(central_longitude=-75, central_latitude=42, globe=None)
desired_proj = ccrs.PlateCarree()

class BasinClassEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.float):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, pyedge):
            return obj.lEdgeID
        if isinstance(obj, pyvertex):
            return json.loads(obj.tojson()) #lVertexID
        if isinstance(obj, pyflowline):
            return obj.lFlowlineID
        
        return JSONEncoder.default(self, obj)


class pybasin(object):
    lBasinID =1 
    sBasinID=''
    lCellID_outlet=-1
    iFlag_disconnected =0
    iFlag_dam=0
    dLongitude_outlet_degree = -9999.
    dLatitude_outlet_degree = -9999.
    dAccumulation_threshold= 100000.0
    dThreshold_small_river = 10000
    dLength_flowline_total = 0.0
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

    aFlowline_basin=None
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

        self.sFilename_flowline_segment_index_before_intersect = 'flowline_segment_index_before_intersect_' + sBasinID + '.json'
        self.sFilename_flowline_segment_order_before_intersect = 'flowline_segment_order_before_intersect_' + sBasinID + '.json'
        self.sFilename_flowline_intersect  = 'flowline_intersect_' + sBasinID + '.json'
        self.sFilename_flowline_final = 'flowline_final_' + sBasinID + '.json'
        
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
    
    def export_flowline(self, aFlowline_in, sFilename_json_in,iFlag_projected_in = None,  pSpatial_reference_in = None):

        export_flowline_to_json(aFlowline_in, sFilename_json_in,\
            iFlag_projected_in= iFlag_projected_in, \
            pSpatial_reference_in = pSpatial_reference_in)



    def calculate_flowline_length(self, aFlowline_in):

        dLength = 0.0

        nflowline = len(aFlowline_in)

        for i in range(nflowline):

            pFlowline= aFlowline_in[i]

            pFlowline.calculate_length()

            dLength = dLength + pFlowline.dLength

        self.dLength_flowline_total = dLength

        return dLength
    
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

    def preprocess_flowline(self):

        try:
            sFilename_flowline_filter = self.sFilename_flowline_filter
            sFilename_flowline_filter_json = self.sFilename_flowline_filter_json
            aFlowline_basin, pSpatial_reference = read_flowline_geojson( sFilename_flowline_filter_json )   
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
                        for k in range(len(aFlowline_basin)):
                            if aFlowline_basin[k].lNHDPlusID == lNHDPlusID:
                                aFlowline_basin.pop(k)
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
                aFlowline_basin = aFlowline_basin + aFlowline_dams_headwater + aFlowline_dams_nonheadwater
            else:
                pass

            if self.iFlag_disconnected == 1:                
                #aThreshold = np.full(2, 300.0, dtype=float)
                #aFlowline_basin = connect_disconnect_flowline(aFlowline_basin, aVertex, aThreshold)
                #sFilename_out = 'flowline_connect.json'
                #sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)    
                #export_flowline_to_json(iFlag_projected, aFlowline_basin,pSpatial_reference_gcs, sFilename_out)
                pass
            else:
                pass

            sFilename_out = 'flowline_before_intersect_' + self.sBasinID + '.json'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)            
            self.export_flowline(aFlowline_basin, sFilename_out)

            #calculate length
            self.calculate_flowline_length(aFlowline_basin)
            

            aVertex = find_flowline_vertex(aFlowline_basin)
            sFilename_out = 'flowline_vertex_without_confluence_before_intersect_' + self.sBasinID + '.json'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_vertex_to_json( aVertex, sFilename_out)

            aFlowline_basin = split_flowline(aFlowline_basin, aVertex)
            sFilename_out = 'flowline_split_by_point_before_intersect_' + self.sBasinID + '.json'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_flowline_to_json(aFlowline_basin, sFilename_out)

            #ues location to find outlet

            point= dict()   
            point['dLongitude_degree'] = self.dLongitude_outlet_degree
            point['dLatitude_degree'] = self.dLatitude_outlet_degree
            pVertex_outlet=pyvertex(point)

            aFlowline_basin = correct_flowline_direction(aFlowline_basin,  pVertex_outlet )

            pVertex_outlet = aFlowline_basin[0].pVertex_end

            sFilename_out = 'flowline_direction_before_intersect_' + self.sBasinID + '.json'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_flowline_to_json( aFlowline_basin,  sFilename_out)

            #step 4: remove loops

            aFlowline_basin = remove_flowline_loop(aFlowline_basin)    
            sFilename_out = 'flowline_loop_before_intersect_' + self.sBasinID + '.json'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_flowline_to_json( aFlowline_basin, sFilename_out)

            #using loop to remove small river, here we use 5 steps

            for i in range(3):
                sStep = "{:02d}".format(i+1)
                aFlowline_basin = remove_small_river(aFlowline_basin, self.dThreshold_small_river)
                sFilename_out = 'flowline_large_'+ sStep +'_before_intersect_' + self.sBasinID + '.json'
                sFilename_out =os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_json( aFlowline_basin,  sFilename_out)


                aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity = find_flowline_confluence(aFlowline_basin,  pVertex_outlet)
                sFilename_out = 'flowline_vertex_with_confluence_'+ sStep +'_before_intersect_' + self.sBasinID + '.json'
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_vertex_to_json( aVertex,  sFilename_out, aAttribute_data=aConnectivity)

                aFlowline_basin = merge_flowline( aFlowline_basin,aVertex, pVertex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence  )  
                sFilename_out = 'flowline_merge_'+ sStep +'_before_intersect_' + self.sBasinID + '.json'
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_json( aFlowline_basin,  sFilename_out)

                if len(aFlowline_basin) ==1:
                    break
            
            #the final vertex info
            aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity = find_flowline_confluence(aFlowline_basin,  pVertex_outlet)
            sFilename_out = 'flowline_vertex_with_confluence_before_intersect_final_' + self.sBasinID + '.json'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_vertex_to_json( aVertex,  sFilename_out, aAttribute_data=aConnectivity)

            
            dLength_total_new = self.calculate_flowline_length(aFlowline_basin)
            print(dLength_total_new)

            #build segment index
            aFlowline_basin, aStream_segment = define_stream_segment_index(aFlowline_basin)
            sFilename_out = self.sFilename_flowline_segment_index_before_intersect
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_flowline_to_json(  aFlowline_basin, sFilename_out, \
                aAttribute_data=[aStream_segment], aAttribute_field=['iseg'], aAttribute_dtype=['int'])

            #build stream order 
            aFlowline_basin, aStream_order = define_stream_order(aFlowline_basin)
            sFilename_out = self.sFilename_flowline_segment_order_before_intersect
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_flowline_to_json(  aFlowline_basin, sFilename_out, \
                aAttribute_data=[aStream_segment, aStream_order], aAttribute_field=['iseg','iord'], aAttribute_dtype=['int','int'])
        except:
            print('Flowline preprocess failed')

        self.aFlowline_basin= aFlowline_basin
        return aFlowline_basin

    def intersect_flowline_with_mesh(self, iMesh_type, sFilename_mesh):

        try:
            sWorkspace_output_basin = self.sWorkspace_output_basin
            sFilename_flowline = self.sFilename_flowline_segment_order_before_intersect
            sFilename_flowline_in = os.path.join(sWorkspace_output_basin, sFilename_flowline)
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
            aFlowline_basin, aFlowline_no_parallel, lCellID_outlet, pVertex_outlet \
                = remove_returning_flowline(iMesh_type, aCell_intersect_basin, pVertex_outlet_initial)
            sFilename_out = 'flowline_simplified_after_intersect_' + self.sBasinID + '.json'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)  

            export_flowline_to_json(aFlowline_basin,  sFilename_out)

            #added start
            aFlowline_basin, aEdge = split_flowline_to_edge(aFlowline_basin)

            aFlowline_basin = remove_duplicate_flowline(aFlowline_basin)
            aFlowline_basin = correct_flowline_direction(aFlowline_basin,  pVertex_outlet )
            aFlowline_basin = remove_flowline_loop(  aFlowline_basin )  

            sFilename_out = 'flowline_debug_' + self.sBasinID + '.json'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_flowline_to_json( aFlowline_basin,  sFilename_out)

            aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity\
                = find_flowline_confluence(aFlowline_basin,  pVertex_outlet)

            sFilename_out = 'flowline_vertex_with_confluence_01_after_intersect_' + self.sBasinID + '.json'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_vertex_to_json( aVertex,  sFilename_out, aAttribute_data=aConnectivity)


            aFlowline_basin = merge_flowline( aFlowline_basin,aVertex, pVertex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence  )  

            aFlowline_basin = remove_flowline_loop(  aFlowline_basin )    

            aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity\
                = find_flowline_confluence(aFlowline_basin,  pVertex_outlet)

            aFlowline_basin = merge_flowline( aFlowline_basin,aVertex, pVertex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence  ) 

            aFlowline_basin, aStream_segment = define_stream_segment_index(aFlowline_basin)
            aFlowline_basin, aStream_order = define_stream_order(aFlowline_basin)

            sFilename_out = self.sFilename_flowline_final
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_flowline_to_json( aFlowline_basin,  sFilename_out)

            self.aFlowline_basin = aFlowline_basin

            self.lCellID_outlet = lCellID_outlet
            self.dLongitude_outlet_degree = pVertex_outlet.dLongitude_degree
            self.dLatitude_outlet_degree = pVertex_outlet.dLatitude_degree
        except:
            print('Intersection failed')
        return

    def evaluate(self):

        sFilename_simplified =  self.sFilename_flowline_segment_order_before_intersect
        sFilename_simplified= os.path.join(self.sWorkspace_output_basin, sFilename_simplified)

        sFilename_final = self.sFilename_flowline_final
        sFilename_final= os.path.join(self.sWorkspace_output_basin, sFilename_final)

        
        #intersect first
        sFilename_output= os.path.join(self.sWorkspace_output_basin, 'intersect_flowline.json')
        intersect_flowline_with_flowline(sFilename_simplified, sFilename_final, sFilename_output)

        #plot diff
        sFilename_output = os.path.join(self.sWorkspace_output_basin, 'area_of_diff.png')

        #calculate_area_of_difference_simplified(sFilename_simplified , sFilename_final, sFilename_output)
        return