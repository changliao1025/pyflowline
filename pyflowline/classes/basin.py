import os, sys
from pathlib import Path
import json
from json import JSONEncoder
import importlib
import numpy as np

from pyflowline.classes.timer import pytimer
from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.edge import pyedge
from pyflowline.classes.flowline import pyflowline
from pyflowline.classes.confluence import pyconfluence
from pyflowline.formats.read_flowline import read_flowline_geojson
from pyflowline.formats.read_nhdplus_flowline_shapefile import  read_nhdplus_flowline_geojson_attribute
from pyflowline.formats.read_nhdplus_flowline_shapefile import extract_nhdplus_flowline_shapefile_by_attribute
from pyflowline.formats.read_nhdplus_flowline_shapefile import track_nhdplus_flowline
from pyflowline.formats.convert_flowline_to_geojson import convert_flowline_to_geojson
from pyflowline.formats.export_flowline import export_flowline_to_geojson
from pyflowline.formats.export_vertex import export_vertex_to_geojson
from pyearth.toolbox.reader.text_reader_string import text_reader_string

from pyflowline.algorithms.split.find_flowline_vertex import find_flowline_vertex
from pyflowline.algorithms.split.find_flowline_confluence import find_flowline_confluence
from pyflowline.algorithms.split.split_flowline import split_flowline
from pyflowline.algorithms.split.split_flowline_to_edge import split_flowline_to_edge
from pyflowline.algorithms.split.split_by_length import split_flowline_by_length
from pyflowline.algorithms.merge.merge_flowline import merge_flowline
from pyflowline.algorithms.direction.correct_flowline_direction import correct_flowline_direction
from pyflowline.algorithms.loop.remove_flowline_loop import remove_flowline_loop
from pyflowline.algorithms.simplification.remove_small_river import remove_small_river
from pyflowline.algorithms.simplification.remove_returning_flowline import remove_returning_flowline
from pyflowline.algorithms.simplification.remove_duplicate_flowline import remove_duplicate_flowline
from pyflowline.algorithms.index.define_stream_segment_index import define_stream_segment_index
from pyflowline.algorithms.index.define_stream_topology import define_stream_topology
from pyflowline.algorithms.index.define_stream_order import define_stream_order, update_head_water_stream_order
from pyflowline.algorithms.intersect.intersect_flowline_with_mesh import intersect_flowline_with_mesh
from pyflowline.algorithms.intersect.intersect_flowline_with_flowline import intersect_flowline_with_flowline
from pyflowline.algorithms.auxiliary.calculate_area_of_difference import calculate_area_of_difference_simplified

iFlag_cython = importlib.util.find_spec("cython")
if iFlag_cython is not None:

    from pyflowline.algorithms.cython.kernel import find_vertex_in_list
    from pyflowline.external.tinyr.tinyr.tinyr import RTree
    iFlag_use_rtree = 1
else:
    from pyflowline.algorithms.auxiliary.find_vertex_in_list import find_vertex_in_list
    iFlag_use_rtree =0

iFlag_kml = importlib.util.find_spec("simplekml")
if iFlag_kml is not None:
    #from pyearth.gis.kml.convert_geojson_to_kml import convert_geojson_to_kml
    pass


sys.setrecursionlimit(100000)

class BasinClassEncoder(JSONEncoder):
    """Basin class encoder

    Args:
        JSONEncoder (_type_): _description_
    """
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
            return json.loads(obj.tojson())
        if isinstance(obj, pyedge):
            return obj.lEdgeID
        if isinstance(obj, pyflowline):
            return obj.lFlowlineID
        if isinstance(obj, pyconfluence):
            return obj.dAngle_upstream

        return JSONEncoder.default(self, obj)

class pybasin(object):
    """Basin class

    Args:
        (object): None

    Returns:
        None: A basin object
    """
    lBasinID =1
    sBasinID=''
    lCellID_outlet=-1
    iFlag_debug = 0
    iFlag_disconnected =0
    iFlag_remove_small_river = 0
    iFlag_remove_low_order_river = 0
    iFlag_dam=0

    iFlag_break_by_distance = 0
    dLongitude_outlet_degree = -9999.
    dLatitude_outlet_degree = -9999.
    dAccumulation_threshold= 100000.0
    dThreshold_small_river = 10000
    dLength_flowline_filtered = 0.0
    dLength_flowline_simplified = 0.0
    dLength_flowline_conceptual = 0.0

    dArea_of_difference=0.0
    dDistance_displace = 0.0
    dThreshold_break_by_distance = 5000.0
    sWorkspace_output_basin=''
    sFilename_flowline_raw=''
    sFilename_flowline_filter=''
    sFilename_flowline_filter_geojson=''
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
    aFlowline_basin_edge = None
    pVertex_outlet=None
    aConfluence_basin_simplified= None
    aConfluence_basin_conceptual= None


    #
    sFilename_mesh = ''
    iMesh_type = 0
    sMesh_type = ''
    #json
    sFilename_watershed_json=''
    sFilename_stream_edge_json =''

    #geojson for hexwatershed compatibility
    sFilename_elevation=''
    sFilename_slope=''
    sFilename_drainage_area=''
    sFilename_flow_direction =''
    sFilename_distance_to_outlet = ''
    sFilename_stream_segment=''
    sFilename_stream_edge=''

    sFilename_variable_polygon=''
    sFilename_variable_polyline=''

    pRTree_flowline = None
    pRTree_edge = None


    iFlag_visual = importlib.util.find_spec("cartopy")
    if iFlag_visual is not None:
        from ._visual_basin import basin_plot
        from ._visual_basin import _plot_polyline_variable
        from ._visual_basin import _plot_polygon_variable

        from ._visual_basin import _plot_area_of_difference
    else:
        pass

    def __init__(self, aParameter):
        """
        Initialize the basin class object

        Args:
            aParameter (dict): Dictionary for parameters
        """

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

        if 'iFlag_remove_small_river' in aParameter:
            self.iFlag_remove_small_river             = int(aParameter['iFlag_remove_small_river'])
        else:
            self.iFlag_remove_small_river   = 0

        if 'iFlag_remove_low_order_river' in aParameter:
            self.iFlag_remove_low_order_river = int(aParameter['iFlag_remove_low_order_river'])
        else:
            self.iFlag_remove_low_order_river = 0

        if 'iFlag_dam' in aParameter:
            self.iFlag_dam             = int(aParameter['iFlag_dam'])
        else:
            self.iFlag_dam   = 0

        if 'iFlag_debug' in aParameter:
            self.iFlag_debug             = int(aParameter['iFlag_debug'])
        else:
            self.iFlag_debug   = 0

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
            self.sFilename_flowline_filter = ''

        if 'sWorkspace_output_basin' in aParameter:
            self.sWorkspace_output_basin = aParameter['sWorkspace_output_basin']
        else:
            self.sWorkspace_output_basin = '.'

        Path(self.sWorkspace_output_basin).mkdir(parents=True, exist_ok=True)

        self.sFilename_flowline_filter_geojson = os.path.join(str(self.sWorkspace_output_basin ), "flowline_filter.geojson"  )

        if 'sFilename_dam' in aParameter:
            self.sFilename_dam = aParameter['sFilename_dam']
        else:
            self.sFilename_dam  = ''

        if 'sFilename_flowline_topo' in aParameter:
            self.sFilename_flowline_topo = aParameter['sFilename_flowline_topo']
        else:
            self.sFilename_flowline_topo = ''

        if 'sMesh_type' in aParameter:
            self.sMesh_type =  aParameter['sMesh_type']
        else:
            self.sMesh_type = 'hexagon'

        sMesh_type = self.sMesh_type
        if sMesh_type =='hexagon': #hexagon
            self.iMesh_type = 1
        else:
            if sMesh_type =='square': #square
                self.iMesh_type = 2
            else:
                if sMesh_type =='latlon': #latlon
                    self.iMesh_type = 3
                else:
                    if sMesh_type =='mpas': #mpas
                        self.iMesh_type = 4
                    else:
                        if sMesh_type =='dggrid':
                            self.iMesh_type = 5
                        else:
                            if sMesh_type =='tin': #
                                self.iMesh_type = 6
                            else:
                                print('Unsupported mesh type?')

        self.sBasinID  = "{:08d}".format(self.lBasinID)

        if not os.path.isfile(self.sFilename_flowline_filter):
            print("The filtered flowline file does not exist!")
            pass
        if self.iFlag_dam==1:
            if not os.path.isfile(self.sFilename_flowline_raw):
                print("The raw flowline file does not exist!")
                exit
            if not os.path.isfile(self.sFilename_flowline_topo):
                print("The flowline topology file does not exist!")
                exit
            pass

        #json
        self.sFilename_watershed_json = os.path.join(str(self.sWorkspace_output_basin ), "watershed.json" )
        self.sFilename_stream_edge_json = os.path.join(str(self.sWorkspace_output_basin ), 'stream_edge.json')
        self.sFilename_basin_info = os.path.join(str(self.sWorkspace_output_basin ), 'basin_info.json')
        self.sFilename_flowline_conceptual_info = os.path.join(str(self.sWorkspace_output_basin ), 'flowline_conceptual_info.json')
        self.sFilename_flowline_simplified_info = os.path.join(str(self.sWorkspace_output_basin ), 'flowline_simplified_info.json')
        self.sFilename_confluence_conceptual_info = os.path.join(str(self.sWorkspace_output_basin ),'confluence_conceptual_info.json')
        self.sFilename_confluence_simplified_info = os.path.join(str(self.sWorkspace_output_basin ),'confluence_simplified_info.json')

        #geojson, full path of the file
        #full paths are required for the following files
        #for hexwatershed compatibility, the geojson files will be generated by the pyhexwatershed front end
        self.sFilename_flowline_segment_index_before_intersect = os.path.join(str(self.sWorkspace_output_basin ),'flowline_segment_index_before_intersect.geojson')
        self.sFilename_flowline_simplified = os.path.join(str(self.sWorkspace_output_basin ),'flowline_simplified.geojson')
        self.sFilename_flowline_split = os.path.join(str(self.sWorkspace_output_basin ),'flowline_split.geojson')
        self.sFilename_flowline_intersect  = os.path.join(str(self.sWorkspace_output_basin ),'flowline_intersect_mesh.geojson')
        self.sFilename_flowline_conceptual = os.path.join(str(self.sWorkspace_output_basin ),'flowline_conceptual.geojson')
        self.sFilename_flowline_edge = os.path.join(str(self.sWorkspace_output_basin ),'flowline_edge.geojson')
        self.sFilename_area_of_difference = os.path.join(str(self.sWorkspace_output_basin ),'area_of_difference.geojson')

        self.sFilename_elevation = os.path.join(str(self.sWorkspace_output_basin ), "elevation.geojson" )
        self.sFilename_slope = os.path.join(str(self.sWorkspace_output_basin ), "slope.geojson" )
        self.sFilename_drainage_area =  os.path.join(str(self.sWorkspace_output_basin ), "drainage_area.geojson" )
        self.sFilename_flow_direction = os.path.join(str(self.sWorkspace_output_basin ), "flow_direction.geojson" )
        self.sFilename_distance_to_outlet = os.path.join(str(self.sWorkspace_output_basin ), "distance_to_outlet.geojson" )
        self.sFilename_stream_edge = os.path.join(str(self.sWorkspace_output_basin ), "stream_edge.geojson" )
        self.sFilename_stream_segment = os.path.join(str(self.sWorkspace_output_basin ), "stream_segment.geojson" )
        self.sFilename_variable_polygon = os.path.join(str(self.sWorkspace_output_basin ), "variable_polygon.geojson" )
        self.sFilename_variable_polyline = os.path.join(str(self.sWorkspace_output_basin ), "variable_polyline.geojson" )

        #kml
        self.sFilename_flowline_conceptual_kml = os.path.join(str(self.sWorkspace_output_basin ),'flowline_conceptual.kml')
        return

    def basin_flowline_simplification(self):
        """
        Run the basin flowline simplification

        Returns:
            list [pyflowline]: A list of simplified flowline
        """
        print('Start flowline simplification:',  self.sBasinID)
        sWorkspace_output_basin = self.sWorkspace_output_basin

        ptimer = pytimer()

        if self.iFlag_dam == 1:
            sFilename_flowline_filter = self.sFilename_flowline_filter
            aFlowline_basin_filtered_raw, pSpatial_reference = read_flowline_geojson( sFilename_flowline_filter )
            aVertex_filtered = find_flowline_vertex(aFlowline_basin_filtered_raw)

            ptimer.start()
            nFlowline_before = len(aFlowline_basin_filtered_raw)
            sFilename_dam = self.sFilename_dam
            #obtain dam lookup table C
            aData_dam = text_reader_string(sFilename_dam, iSkipline_in =1,cDelimiter_in=',' )
            sFilename_flowline_topo = self.sFilename_flowline_topo
            #obtain whole topology B
            aData_flowline_topo = text_reader_string(sFilename_flowline_topo, iSkipline_in =1,cDelimiter_in=',' )
            aFromFlowline = aData_flowline_topo[:,1].astype(int).ravel()
            aToFlowline = aData_flowline_topo[:,2].astype(int).ravel()
            sFilename_flowline_raw = self.sFilename_flowline_raw
            #find A lookup table
            aNHDPlusID_filter = read_nhdplus_flowline_geojson_attribute(sFilename_flowline_filter)
            ndam = len(aData_dam)
            aNHDPlusID_dams_headwater = list()
            aFlowline_dams_nonheadwater = list()
            aVertex_dams_nonheadwater=list()
            n=0
            for j in range(0, ndam):
                dLon = float(aData_dam[j][1])
                dLat = float(aData_dam[j][0])
                sDam = aData_dam[j][4]
                #individual ID
                lNHDPlusID = int(aData_dam[j][5])
                #if lNHDPlusID in aNHDPlusID_filter:
                if lNHDPlusID in aNHDPlusID_filter: #this is already included in A
                    #change flag by id
                    for k in range(len(aFlowline_basin_filtered_raw)):
                        #remove this flowline by ID
                        if aFlowline_basin_filtered_raw[k].lNHDPlusID == lNHDPlusID:
                            aFlowline_basin_filtered_raw[k].iFlag_dam =1
                            break
                    pass
                else:
                    aNHDPlusID_dams_headwater.append(lNHDPlusID)
                    #not in A, so we need to trace it down

                    aNHDPlusID_dam_nonheadwater = track_nhdplus_flowline(aNHDPlusID_filter, aFromFlowline, aToFlowline, lNHDPlusID)
                    aFlowline_dam_nonheadwater = extract_nhdplus_flowline_shapefile_by_attribute(sFilename_flowline_raw, aNHDPlusID_dam_nonheadwater )
                    #clean up
                    aVertex_dam_nonheadwater = find_flowline_vertex(aFlowline_dam_nonheadwater)
                    dThreshold = 1.0
                    iFlag_found=0
                    n=n+len(aFlowline_dam_nonheadwater)
                    #print(j,n)
                    for pVertex in aVertex_filtered:
                        if iFlag_found ==1:
                            break

                        for i in range( len(aVertex_dam_nonheadwater) ):
                            pVertex_dam = aVertex_dam_nonheadwater[i]
                            dDistance = pVertex.calculate_distance(pVertex_dam)
                            if  dDistance<= dThreshold:
                                aVertex_dam_nonheadwater[i] = pVertex
                                for k in range(len(aFlowline_dam_nonheadwater)):
                                    if aFlowline_dam_nonheadwater[k].pVertex_end == pVertex_dam:
                                        #update
                                        aEdge = aFlowline_dam_nonheadwater[k].aEdge
                                        aEdge[-1] = pyedge( aEdge[-1].pVertex_start, pVertex)
                                        aFlowline_dam_nonheadwater[k] = pyflowline(aEdge)

                                iFlag_found = 1
                                break
                            else:
                                #print(dDistance)
                                pass

                    aVertex_dams_nonheadwater.append(aVertex_dam_nonheadwater)

                    aFlowline_dams_nonheadwater.append(aFlowline_dam_nonheadwater)


            aFlowline_dams_headwater = extract_nhdplus_flowline_shapefile_by_attribute(sFilename_flowline_raw, aNHDPlusID_dams_headwater )
            for i in range(len(aFlowline_dams_headwater)):
                aFlowline_dams_headwater[i].iFlag_dam = 1

            aFlowline_dams_nonheadwater_all = [item for sublist in aFlowline_dams_nonheadwater for item in sublist]
            aVertex_dam_nonheadwater_all = [item for sublist in aVertex_dams_nonheadwater for item in sublist]

            aFlowline_basin_filtered = aFlowline_basin_filtered_raw + aFlowline_dams_headwater + aFlowline_dams_nonheadwater_all
            aVertex_dam = find_flowline_vertex(aFlowline_dams_headwater)

            if self.iFlag_debug ==1:
                sFilename_out = 'flowline_vertex_filtered.geojson'
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_vertex_to_geojson( aVertex_filtered, sFilename_out)
                sFilename_out = 'flowline_vertex_dam.geojson'
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_vertex_to_geojson( aVertex_dam, sFilename_out)
                sFilename_out = 'flowline_vertex_dam_nonheadwater.geojson'
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_vertex_to_geojson( aVertex_dam_nonheadwater_all, sFilename_out)

            nFlowline_after = len(aFlowline_basin_filtered)
            print('Basin ',  self.sBasinID, ' has dam', nFlowline_before, nFlowline_after)
            ptimer.stop()
        else:
            print('Basin ',  self.sBasinID, ' has no dam')
            sFilename_flowline_filter = self.sFilename_flowline_filter
            aFlowline_basin_filtered, pSpatial_reference = read_flowline_geojson( sFilename_flowline_filter )
            #aVertex_filtered = find_flowline_vertex(aFlowline_basin_filtered)

            pass
        sys.stdout.flush()
        if self.iFlag_disconnected == 1:
            #not used anymore
            #aThreshold = np.full(2, 300.0, dtype=float)
            #aFlowline_basin_filtered = connect_disconnect_flowline(aFlowline_basin_filtered, aVertex, aThreshold)
            #sFilename_out = 'flowline_connect.geojson'
            #sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            #export_flowline_to_geojson(iFlag_projected, aFlowline_basin_filtered,pSpatial_reference_gcs, sFilename_out)
            pass


        if self.iFlag_debug ==1:
            sFilename_out = 'flowline_before_intersect.geojson'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            self.basin_export_flowline(aFlowline_basin_filtered, sFilename_out)
        #calculate length
        self.aFlowline_basin_filtered = aFlowline_basin_filtered
        self.dLength_flowline_filtered = self.basin_calculate_flowline_length(aFlowline_basin_filtered)

        #assign vertex id


        #simplification started
        print('Basin ',  self.sBasinID, 'find flowline vertex')
        sys.stdout.flush()
        ptimer.start()
        aVertex = find_flowline_vertex(aFlowline_basin_filtered)
        ptimer.stop()
        if self.iFlag_debug ==1:
            sFilename_out = 'flowline_vertex_without_confluence_before_intersect.geojson'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_vertex_to_geojson( aVertex, sFilename_out)

        try:
            print('Basin ', self.sBasinID, 'split flowline')
            sys.stdout.flush()
            ptimer.start()
       
            nFlowline_before = len(aFlowline_basin_filtered)
            aFlowline_basin_simplified = split_flowline(aFlowline_basin_filtered, aVertex)
            nFlowline_after = len(aFlowline_basin_simplified)
            ptimer.stop()
        except:
            print(nFlowline_before)


        if self.iFlag_debug ==1:
            sFilename_out = 'flowline_split_by_point_before_intersect.geojson'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_flowline_to_geojson(aFlowline_basin_simplified, sFilename_out)

        #use location to find outlet
        point= dict()
        point['dLongitude_degree'] = self.dLongitude_outlet_degree
        point['dLatitude_degree'] = self.dLatitude_outlet_degree
        pVertex_outlet=pyvertex(point)


        try:
            print('Basin ',  self.sBasinID, 'started correction flow direction')
            sys.stdout.flush()
            ptimer.start()

            nFlowline_before = len(aFlowline_basin_simplified)
            #this flowline is ordered from downstream to upstream
            aFlowline_basin_simplified = correct_flowline_direction(aFlowline_basin_simplified,  pVertex_outlet )
            nFlowline_after = len(aFlowline_basin_simplified)
            ptimer.stop()

            pVertex_outlet = aFlowline_basin_simplified[0].pVertex_end
            self.pVertex_outlet = pVertex_outlet
            if self.iFlag_debug ==1:
                sFilename_out = 'flowline_direction_before_intersect.geojson'
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_geojson( aFlowline_basin_simplified,  sFilename_out)
        except:
            print('Error in flow direction correction')



        #step 4: remove loops
        try:
            print('Basin ',  self.sBasinID, 'started loop removal')
            sys.stdout.flush()
            ptimer.start()
            aFlowline_basin_simplified = remove_flowline_loop(aFlowline_basin_simplified)
            ptimer.stop()

            if self.iFlag_debug ==1:
                sFilename_out = 'flowline_loop_before_intersect.geojson'
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_geojson( aFlowline_basin_simplified, sFilename_out)
        except:
            print('Error in loop removal')


        #at this stage, the stream order information is not available, but could be easily rebuild for speed up the small river removal algorithm
        print('Basin ',  self.sBasinID, 'started update stream order initial')
        sys.stdout.flush()
        ptimer.start()
        aFlowline_basin_simplified = update_head_water_stream_order(aFlowline_basin_simplified )
        ptimer.stop()

        #using loop to remove small river, here we use 3 steps
        if self.iFlag_remove_small_river ==1:
            for i in np.arange(0, 3, 1):
                sStep = "{:02d}".format(i+1)

                dThreshold = self.dThreshold_small_river
                print('Basin ',  self.sBasinID, 'started small river removal', sStep, dThreshold)
                sys.stdout.flush()
                ptimer.start()

                aFlowline_basin_simplified = remove_small_river(aFlowline_basin_simplified, dThreshold )
                if self.iFlag_debug ==1:
                    sFilename_out = 'flowline_large_'+ sStep +'_before_intersect.geojson'
                    sFilename_out =os.path.join(sWorkspace_output_basin, sFilename_out)
                    export_flowline_to_geojson( aFlowline_basin_simplified,  sFilename_out)
                    #print(len(aFlowline_basin_simplified))

                #update stream order
                aFlowline_basin_simplified = update_head_water_stream_order(aFlowline_basin_simplified )

                aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity, pVertex_outlet = find_flowline_confluence(aFlowline_basin_simplified,  pVertex_outlet)
                if self.iFlag_debug ==1:
                    sFilename_out = 'flowline_vertex_with_confluence_'+ sStep +'_before_intersect.geojson'
                    sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                    export_vertex_to_geojson( aVertex,  sFilename_out, aAttribute_data=aConnectivity)

                aFlowline_basin_simplified = merge_flowline( aFlowline_basin_simplified, aVertex, pVertex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence  )
                if self.iFlag_debug ==1:
                    sFilename_out = 'flowline_merge_'+ sStep +'_before_intersect.geojson'
                    sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                    export_flowline_to_geojson( aFlowline_basin_simplified,  sFilename_out)

                aFlowline_basin_simplified = update_head_water_stream_order(aFlowline_basin_simplified )

                ptimer.stop()
                if len(aFlowline_basin_simplified) == 1:
                    break

            sys.stdout.flush()
        else: #if we dont remove small river, we still need to merge the flowline
            aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity, pVertex_outlet = find_flowline_confluence(aFlowline_basin_simplified,  pVertex_outlet)
            if self.iFlag_debug ==1:
                sFilename_out = 'flowline_vertex_with_confluence_before_intersect.geojson'
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_vertex_to_geojson( aVertex,  sFilename_out, aAttribute_data=aConnectivity)
            aFlowline_basin_simplified = merge_flowline( aFlowline_basin_simplified, aVertex, pVertex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence  )
            if self.iFlag_debug ==1:
                sFilename_out = 'flowline_merge_before_intersect.geojson'
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_geojson( aFlowline_basin_simplified,  sFilename_out)
            aFlowline_basin_simplified = update_head_water_stream_order(aFlowline_basin_simplified )
            pass

        #the final vertex info
        print('Basin ',  self.sBasinID, 'find flowline confluence')
        sys.stdout.flush()
        ptimer.start()
        aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity, pVertex_outlet = find_flowline_confluence(aFlowline_basin_simplified,  pVertex_outlet)
        ptimer.stop()
        sFilename_out = 'vertex_simplified.geojson'
        sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
        export_vertex_to_geojson( aVertex,  sFilename_out, aAttribute_data=aConnectivity)

        #starting here, the self object can be updated because no more edit (add/remove) will be made to the flowline
        #but attributes will be updated

        #build segment index, they are reversed
        print('Basin ',  self.sBasinID, 'started stream segment definition')
        sys.stdout.flush()
        ptimer.start()
        aFlowline_basin_simplified, aStream_segment = define_stream_segment_index(aFlowline_basin_simplified)
        ptimer.stop()
        if self.iFlag_debug ==1:
            sFilename_out = self.sFilename_flowline_segment_index_before_intersect
            export_flowline_to_geojson(aFlowline_basin_simplified,
                                       sFilename_out,
                aAttribute_data=[aStream_segment],
                aAttribute_field=['stream_segment'],
                aAttribute_dtype=['int'])

        aVertex = np.array(aVertex)
        aIndex_confluence = np.array(aIndex_confluence)
        if aIndex_confluence.size > 0:
            print('Basin ',  self.sBasinID, 'started confluence definition')
            ptimer.start()
            aVertex_confluence = aVertex[aIndex_confluence]
            aConfluence_basin_simplified = self.basin_build_confluence(aFlowline_basin_simplified, aVertex_confluence)
            ptimer.stop()

        #change added a new function to build stream topology
        print('Basin ',  self.sBasinID, 'started stream topology definition')
        ptimer.start()
        aFlowline_basin_simplified = define_stream_topology(aFlowline_basin_simplified, aConfluence_basin_simplified)
        ptimer.stop()

        #build stream order
        print('Basin ',  self.sBasinID, 'started stream order definition')
        ptimer.start()
        aFlowline_basin_simplified, aStream_order = define_stream_order(aFlowline_basin_simplified, aConfluence_basin_simplified)
        ptimer.stop()
        sFilename_out = self.sFilename_flowline_simplified
        export_flowline_to_geojson(aFlowline_basin_simplified,
                                   sFilename_out,
                aAttribute_data=[aStream_segment, aStream_order],
                aAttribute_field=['stream_segment','stream_order'],
                aAttribute_dtype=['int','int'])

        if self.iFlag_break_by_distance==1:
            print('Basin ',  self.sBasinID, 'started flowline split by length')
            ptimer.start()
            aFlowline_basin_simplified_split = split_flowline_by_length(aFlowline_basin_simplified, self.dThreshold_break_by_distance)
            ptimer.stop()

            sFilename_out = self.sFilename_flowline_split
            export_flowline_to_geojson(  aFlowline_basin_simplified_split, sFilename_out  )

        self.aConfluence_basin_simplified = aConfluence_basin_simplified
        self.aFlowline_basin_simplified= aFlowline_basin_simplified
        print('Finish flowline simplification:',  self.sBasinID)
        sys.stdout.flush()
        return aFlowline_basin_simplified

    def basin_reconstruct_topological_relationship(self, iMesh_type, sFilename_mesh):
        """
        Run the basin topologic relationship reconstruction

        Args:
            iMesh_type (int): Mesh type
            sFilename_mesh (str): Filename of the geojson mesh

        Returns:
            list [pyflowline]: A list of intersected cells
        """
        print('Basin ',  self.sBasinID, 'Start topology reconstruction')

        ptimer = pytimer()

        sWorkspace_output_basin = self.sWorkspace_output_basin
        sFilename_flowline_in = self.sFilename_flowline_simplified
        sFilename_flowline_intersect_out = self.sFilename_flowline_intersect

        try:
            print('Basin ',  self.sBasinID, 'Start flowline and mesh intersection')
            ptimer.start()
            aCell, aCell_intersect_basin, aFlowline_intersect_all = intersect_flowline_with_mesh(iMesh_type, sFilename_mesh, \
                sFilename_flowline_in, sFilename_flowline_intersect_out)
            ptimer.stop()
            sys.stdout.flush()
            if self.iFlag_debug ==1:
                sFilename_out = 'flowline_intersect_flowline_with_mesh.geojson'
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_geojson(aFlowline_intersect_all,  sFilename_out)
        except:
            print('Error in flowline and mesh intersection.')


        point= dict()
        point['dLongitude_degree'] = self.dLongitude_outlet_degree
        point['dLatitude_degree'] = self.dLatitude_outlet_degree
        pVertex_outlet_initial=pyvertex(point)

        #from this point, aFlowline_basin is conceptual
        #segment based
        try:
            print('Basin ',  self.sBasinID, 'Start return flowline removal')
            ptimer.start()
            aFlowline_basin_conceptual, lCellID_outlet, pVertex_outlet = remove_returning_flowline(iMesh_type, aCell_intersect_basin, pVertex_outlet_initial)
            ptimer.stop()
            sys.stdout.flush()
            if self.iFlag_debug ==1:
                sFilename_out = 'flowline_simplified_after_intersect.geojson'
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_geojson(aFlowline_basin_conceptual, sFilename_out)
        except:
            print('Error in remove_returning_flowline.')


        #edge based
        try:
            print('Basin ',  self.sBasinID, 'Start split flowline to edge')
            ptimer.start()
            aFlowline_basin_conceptual, aEdge = split_flowline_to_edge(aFlowline_basin_conceptual)
            ptimer.stop()
            sys.stdout.flush()
            if self.iFlag_debug ==1:
                sFilename_out = 'flowline_edge_split_flowline_to_edge.geojson'
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_geojson( aFlowline_basin_conceptual,  sFilename_out)
        except:
            print('Error in split_flowline_to_edge.')

        try:
            print('Basin ',  self.sBasinID, 'Start remove duplicate flowline')
            ptimer.start()
            aFlowline_basin_conceptual = remove_duplicate_flowline(aFlowline_basin_conceptual)
            ptimer.stop()
            sys.stdout.flush()
            if self.iFlag_debug ==1:
                sFilename_out = 'flowline_edge_remove_duplicate_flowline.geojson'
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_geojson( aFlowline_basin_conceptual,  sFilename_out)
        except:
            print('Error in remove_duplicate_flowline.')

        try:
            print('Basin ',  self.sBasinID, 'Start flowline direction correction')
            ptimer.start()
            aFlowline_basin_conceptual = correct_flowline_direction(aFlowline_basin_conceptual,  pVertex_outlet )
            ptimer.stop()
            sys.stdout.flush()
            if self.iFlag_debug ==1:
                sFilename_out = 'flowline_edge_correct_flowline_direction.geojson'
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_geojson( aFlowline_basin_conceptual,  sFilename_out)
        except:
            print('Error in correct_flowline_direction.')

        try:
            print('Basin ',  self.sBasinID, 'Start flowline direction correction')
            ptimer.start()
            aFlowline_basin_conceptual = remove_flowline_loop(aFlowline_basin_conceptual )
            ptimer.stop()
            sys.stdout.flush()
            if self.iFlag_debug ==1:
                sFilename_out = 'flowline_edge_remove_flowline_loop.geojson'
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_geojson( aFlowline_basin_conceptual,  sFilename_out)
        except:
            print('Error in remove_flowline_loop.')

        try:
            aFlowline_basin_conceptual = update_head_water_stream_order(aFlowline_basin_conceptual )
        except:
            print('Error in update_head_water_stream_order.')

        try:
            print('Basin ',  self.sBasinID, 'Start find flowline confluence')
            ptimer.start()
            aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity, pVertex_outlet\
                = find_flowline_confluence(aFlowline_basin_conceptual,  pVertex_outlet)
            ptimer.stop()
            sys.stdout.flush()
            if self.iFlag_debug ==1:
                sFilename_out = 'flowline_vertex_with_confluence_after_intersect.geojson'
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_vertex_to_geojson( aVertex,  sFilename_out, aAttribute_data=aConnectivity)
        except:
            print('Error in find_flowline_confluence.')

        #segment based
        try:
            print('Basin ',  self.sBasinID, 'Start merge flowline')
            ptimer.start()
            aFlowline_basin_conceptual = merge_flowline( aFlowline_basin_conceptual,aVertex, pVertex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence  )
            ptimer.stop()
            sys.stdout.flush()
        except:
            print('Error in merge_flowline.')

        aFlowline_basin_conceptual = update_head_water_stream_order(aFlowline_basin_conceptual )

        try:
            print('Basin ',  self.sBasinID, 'Start find flowline confluence')
            ptimer.start()
            aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity, pVertex_outlet\
                = find_flowline_confluence(aFlowline_basin_conceptual,  pVertex_outlet)
            ptimer.stop()
            sys.stdout.flush()
        except:
            print('Error in find_flowline_confluence.')

        aFlowline_basin_conceptual, aStream_segment = define_stream_segment_index(aFlowline_basin_conceptual)

        #save confluence
        aVertex = np.array(aVertex)
        aIndex_confluence = np.array(aIndex_confluence)
        if aIndex_confluence.size > 0:
            aVertex_confluence = aVertex[aIndex_confluence]
            aConfluence_basin_conceptual = self.basin_build_confluence(aFlowline_basin_conceptual, aVertex_confluence)
        else:
            #there is no confluence
            pass

        #change added a new function to build stream topology
        print('Basin ',  self.sBasinID, 'started stream topology definition')
        ptimer.start()
        aFlowline_basin_conceptual = define_stream_topology(aFlowline_basin_conceptual, aConfluence_basin_conceptual)
        ptimer.stop()
        sys.stdout.flush()

        aFlowline_basin_conceptual, aStream_order = define_stream_order(aFlowline_basin_conceptual, aConfluence_basin_conceptual)

        #edge based
        aFlowline_basin_edge, aEdge = split_flowline_to_edge(aFlowline_basin_conceptual)
        sFilename_out = self.sFilename_flowline_edge
        export_flowline_to_geojson( aFlowline_basin_edge, sFilename_out)

        sFilename_out = self.sFilename_flowline_conceptual
        export_flowline_to_geojson( aFlowline_basin_conceptual,
                                   sFilename_out,
            aAttribute_data=[aStream_segment, aStream_order],
            aAttribute_field=['stream_segment','stream_order'],
            aAttribute_dtype=['int','int'])

        self.aConfluence_basin_conceptual = aConfluence_basin_conceptual
        self.aFlowline_basin_conceptual = aFlowline_basin_conceptual
        self.aFlowline_basin_edge = aFlowline_basin_edge

        self.lCellID_outlet = lCellID_outlet
        self.dLongitude_outlet_degree = pVertex_outlet.dLongitude_degree
        self.dLatitude_outlet_degree = pVertex_outlet.dLatitude_degree

        print('Finish topology reconstruction:',  self.sBasinID)
        sys.stdout.flush()
        if iFlag_use_rtree ==1:
            interleaved = True
            self.pRTree_flowline = RTree(interleaved=interleaved, max_cap=5, min_cap=2)
            for lFlowlineIndex in range(len(self.aFlowline_basin_conceptual)):
                pBound = self.aFlowline_basin_conceptual[lFlowlineIndex].pBound
                self.pRTree_flowline.insert(lFlowlineIndex, pBound)

            self.pRTree_edge = RTree(interleaved=interleaved, max_cap=5, min_cap=2)
            for lEdgeIndex in range(len(self.aFlowline_basin_edge)):
                pBound = self.aFlowline_basin_edge[lEdgeIndex].pBound
                self.pRTree_edge.insert(lEdgeIndex, pBound)
            pass
        return aCell_intersect_basin

    def basin_build_confluence(self, aFlowline_basin_in, aVertex_confluence_in):
        """
        Build the conflence

        Args:
            aFlowline_basin_in (list [pyflowline]): A list of flowlines in this basin
            aVertex_confluence_in (list [pyconfluence]): A list of vertices in this basin

        Returns:
            list [pyconfluence]: A list of confluences in this basin
        """
        #this can only be calculated for confluence
        # Create a dictionary to map each vertex to its upstream and downstream flowlines
        vertex_to_flowlines = {}
        for pFlowline in aFlowline_basin_in:
            pVertex_start = pFlowline.pVertex_start
            pVertex_end = pFlowline.pVertex_end

            if pVertex_end not in vertex_to_flowlines:
                vertex_to_flowlines[pVertex_end] = {'upstream': [], 'downstream': None}
            vertex_to_flowlines[pVertex_end]['upstream'].append(pFlowline)

            if pVertex_start not in vertex_to_flowlines:
                vertex_to_flowlines[pVertex_start] = {'upstream': [], 'downstream': None}
            vertex_to_flowlines[pVertex_start]['downstream'] = pFlowline

        # Build the confluence for each vertex
        aConfluence_basin = []
        for pVertex in aVertex_confluence_in:
            aFlowline_upstream = vertex_to_flowlines[pVertex]['upstream']
            pFlowline_downstream = vertex_to_flowlines[pVertex]['downstream']
            pConfluence = pyconfluence(pVertex, aFlowline_upstream, pFlowline_downstream)
            aConfluence_basin.append(pConfluence)

        #aConfluence_basin=list()
        #for pVertex in aVertex_confluence_in:
        #    aFlowline_upstream =list()
        #    for pFlowline in aFlowline_basin_in:
        #        pVertex_start = pFlowline.pVertex_start
        #        pVertex_end = pFlowline.pVertex_end
        #        if pVertex_end == pVertex:
        #            aFlowline_upstream.append(pFlowline)
        #            pass
        #        if pVertex_start == pVertex:
        #            pFlowline_downstream=pFlowline
        #    pConfluence = pyconfluence(pVertex, aFlowline_upstream, pFlowline_downstream)
        #    aConfluence_basin.append(pConfluence)

        return aConfluence_basin


    def basin_analyze(self):
        """
        Analyze the basin results including length, sinuosity, and breaching angle
        """

        if self.aFlowline_basin_filtered is None:
            sFilename_flowline_filter = self.sFilename_flowline_filter
            sFilename_flowline_filter_geojson = self.sFilename_flowline_filter_geojson
            self.aFlowline_basin_filtered, pSpatial_reference = read_flowline_geojson( sFilename_flowline_filter_geojson )
            self.dLength_flowline_filtered = self.calculate_flowline_length(self.aFlowline_basin_filtered)

        point= dict()
        point['dLongitude_degree'] = self.dLongitude_outlet_degree
        point['dLatitude_degree'] = self.dLatitude_outlet_degree
        pVertex_outlet_initial=pyvertex(point)
        if self.aFlowline_basin_simplified is None:
            sFilename_flowline_in = self.sFilename_flowline_simplified
            aFlowline_simplified,pSpatial_reference = read_flowline_geojson( sFilename_flowline_in )

            self.aFlowline_basin_simplified = aFlowline_simplified
            aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity, pVertex_outlet\
            = find_flowline_confluence(self.aFlowline_basin_simplified,  pVertex_outlet_initial)
            aVertex = np.array(aVertex)
            aIndex_confluence = np.array(aIndex_confluence)
            if aIndex_confluence.size > 0:
                aVertex_confluence = aVertex[aIndex_confluence]
                self.aConfluence_basin_simplified = self.build_confluence(self.aFlowline_basin_simplified, aVertex_confluence)

        self.dLength_flowline_simplified = self.calculate_flowline_length(self.aFlowline_basin_simplified)

        if self.aFlowline_basin_conceptual is None:
            sFilename_flowline_in = self.sFilename_flowline_conceptual
            aFlowline_conceptual, pSpatial_reference = read_flowline_geojson( sFilename_flowline_in )
            self.aFlowline_basin_conceptual = aFlowline_conceptual

            aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity, pVertex_outlet\
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

    def basin_export(self):
        """
        Export the basin outputs in json format
        """
        self.basin_export_basin_info_to_json()
        self.basin_export_flowline_info_to_json()
        self.basin_export_confluence_info_to_json()

        if iFlag_kml is not None:
            #only convert final conceptual flowline to kml
            sFilename_conceptual = self.sFilename_flowline_conceptual
            sFilename_conceptual_kml = self.sFilename_flowline_conceptual_kml

            #convert_geojson_to_kml(sFilename_conceptual, sFilename_conceptual_kml)


        return

    def basin_export_flowline(self, aFlowline_in, sFilename_json_in,iFlag_projected_in = None,  pSpatial_reference_in = None):
        """
        Export the basin flowline to geojson

        Args:
            aFlowline_in (list [pyflowline]): A list of flowlines
            sFilename_json_in (str): The output json filename
            iFlag_projected_in (int, optional): Flag if re-projection is needed. Defaults to None.
            pSpatial_reference_in (object, optional): The spatial reference if re-projection is needed. Defaults to None.
        """
        export_flowline_to_geojson(aFlowline_in, sFilename_json_in,
            iFlag_projected_in= iFlag_projected_in,
            pSpatial_reference_in = pSpatial_reference_in)

        return

    def basin_export_basin_info_to_json(self):
        """
        Export the basin basin object to json
        """
        sFilename_json = self.sFilename_basin_info

        aSkip = ['aFlowline_basin_filtered', \
                'aFlowline_basin_simplified','aFlowline_basin_conceptual','aConfluence_basin_simplified',
                'aConfluence_basin_conceptual', 'pRTree_flowline', 'pRTree_edge']
        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)
            pass

        with open(sFilename_json, 'w', encoding='utf-8') as f:
            sJson = json.dumps(obj, default=lambda o: o.__dict__,\
            sort_keys=True, \
                indent = 4, \
                    ensure_ascii=True, \
                        cls=BasinClassEncoder)
            f.write(sJson)
            f.close()
        return

    def basin_export_flowline_info_to_json(self):
        """
        Export the flowline object to json
        """
        iFlag_export_simplified=0
        if iFlag_export_simplified==1:
            sFilename_json = self.sFilename_flowline_simplified_info
            with open(sFilename_json, 'w', encoding='utf-8') as f:
                sJson = json.dumps([json.loads(ob.tojson()) for ob in self.aFlowline_basin_simplified], indent = 4)
                f.write(sJson)
                f.close()

        sFilename_json = self.sFilename_flowline_conceptual_info

        with open(sFilename_json, 'w', encoding='utf-8') as f:
            sJson = json.dumps([json.loads(ob.tojson()) for ob in self.aFlowline_basin_conceptual], indent = 4)
            f.write(sJson)
            f.close()
        return

    def basin_export_confluence_info_to_json(self):
        """
        Export the confluence object to json
        """
        #iFlag_export_confluence =0
        if self.aConfluence_basin_simplified is not None:
            sFilename_json = self.sFilename_confluence_simplified_info

            with open(sFilename_json, 'w', encoding='utf-8') as f:
                sJson = json.dumps([json.loads(ob.tojson()) for ob in self.aConfluence_basin_simplified], indent = 4)
                f.write(sJson)
                f.close()

        if self.aConfluence_basin_conceptual is not None:
            sFilename_json = self.sFilename_confluence_conceptual_info


            with open(sFilename_json, 'w', encoding='utf-8') as f:
                sJson = json.dumps([json.loads(ob.tojson()) for ob in self.aConfluence_basin_conceptual], indent = 4)
                f.write(sJson)
                f.close()

        return

    def tojson(self):
        """
        Export the basin object to json

        Returns:
            json str: A json string
        """
        aSkip = ['aFlowline_basin_filtered', \
                'aFlowline_basin_simplified','aFlowline_basin_conceptual','aConfluence_basin_simplified',
                'aConfluence_basin_conceptual','pRTree_flowline', 'pRTree_edge']

        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)
            pass


        sJson = json.dumps(obj,
            sort_keys=True,
                indent = 4,
                    ensure_ascii=True,
                        cls=BasinClassEncoder)
        return sJson

    def basin_export_config_to_json(self, sFilename_output_in = None):
        """
        Export the basin object to json using the encoder

        Args:
            sFilename_output_in (str, optional): The json filename. Defaults to None.
        """
        #single basin
        if sFilename_output_in is not None:
            sFilename_output = sFilename_output_in
        else:
            sFilename_output = os.path.join(self.sWorkspace_output_basin, 'configuration_basin.json' )

        aSkip = ['aFlowline_basin_filtered', \
                'aFlowline_basin_simplified','aFlowline_basin_conceptual','aConfluence_basin_simplified',
                'aConfluence_basin_conceptual','aFlowline_basin_edge']
        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)
        with open(sFilename_output, 'w', encoding='utf-8') as f:
            json.dump(obj, f,sort_keys=True, \
                ensure_ascii=False, \
                indent=4, \
                cls=BasinClassEncoder)
        return

    def basin_convert_flowline_to_geojson(self):
        """
        Convert the flowline to geojson
        """
        sFilename_raw = self.sFilename_flowline_filter
        sFilename_out = self.sFilename_flowline_filter_geojson
        print('Basin '+ self.sBasinID + ': initial flowline:', sFilename_raw )
        convert_flowline_to_geojson(1, sFilename_raw, sFilename_out)
        return

    def basin_calculate_flowline_length(self, aFlowline_in):
        """
        Calculate the length of flowlines

        Args:
            aFlowline_in (list [pyflowline]): A list of flowlines

        Returns:
            float: The total length of all flowlines
        """
        #dLength = 0.0
        #nflowline = len(aFlowline_in)
        #for i in range(nflowline):
        #    pFlowline= aFlowline_in[i]
        #    pFlowline.calculate_length()
        #    dLength = dLength + pFlowline.dLength
        #return dLength
        return sum(pFlowline.dLength for pFlowline in aFlowline_in)

    def basin_calculate_river_sinuosity(self):
        """
        Calcualte the the river sinuosity
        """
        for pFlowline in self.aFlowline_basin_simplified:
            pFlowline.calculate_flowline_sinuosity()

        for pFlowline in self.aFlowline_basin_conceptual:
            pFlowline.calculate_flowline_sinuosity()

        return

    def basin_calculate_confluence_branching_angle(self):
        """
        Calcualte the the river confluence branching angle
        """
        for pConfluence in self.aConfluence_basin_simplified:
            pConfluence.calculate_branching_angle()
        for pConfluence in self.aConfluence_basin_conceptual:
            pConfluence.calculate_branching_angle()
        return

    def basin_evaluate(self, iMesh_type, sMesh_type):
        """
        Evaluate the model performance

        Args:
            iMesh_type (int): The mesh type
            sMesh_type (str): The mesh type
        """

        self.evaluate_area_of_difference(iMesh_type, sMesh_type)
        return

    def basin_evaluate_area_of_difference(self, iMesh_type, sMesh_type):
        """
        Evaluate the model performance using area of difference

        Args:
            iMesh_type (int): The mesh type
            sMesh_type (str): The mesh type
        """


        sFilename_simplified = self.sFilename_flowline_simplified


        sFilename_flowline_edge = self.sFilename_flowline_edge

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
                aIndex_confluence_simplified, aConnectivity, pVertex_outlet\
                = find_flowline_confluence(aFlowline_simplified,  pVertex_outlet)


        aFlowline_conceptual,pSpatial_reference = read_flowline_geojson( sFilename_flowline_edge )
        aVertex_conceptual, lIndex_outlet_conceptual, \
            aIndex_headwater_conceptual, aIndex_middle_conceptual, \
            aIndex_confluence_conceptual,  aConnectivity, pVertex_outlet \
                = find_flowline_confluence(aFlowline_conceptual,  pVertex_outlet)

        #merge intersect with confluence
        a=np.array(aIndex_confluence_simplified)
        b=np.array(aIndex_confluence_conceptual)
        c=np.array(aIndex_middle_conceptual)
        d = list(np.array(aVertex_simplified)[a] )
        e = list(np.array(aVertex_conceptual)[b] )
        f = list(np.array(aVertex_conceptual)[c] )

        g = aVertex_intersect  + d + e + f

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
            export_vertex_to_geojson( aVertex_all, sFilename_output)

        #split
        aFlowline_simplified_split = split_flowline(aFlowline_simplified, aVertex_all_simplified,iFlag_intersect =1)
        self.iFlag_debug =1
        if self.iFlag_debug ==1:
            sFilename_out = 'flowline_split_simplified.json'
            sFilename_out = os.path.join(self.sWorkspace_output_basin, sFilename_out)
            export_flowline_to_geojson(aFlowline_simplified_split, sFilename_out)

        aFlowline_conceptual_split = split_flowline(aFlowline_conceptual, aVertex_all_conceptual,\
            iFlag_intersect =1, iFlag_use_id=1)
        self.iFlag_debug =1
        if self.iFlag_debug ==1:
            sFilename_out = 'flowline_split_conceptual.json'
            sFilename_out = os.path.join(self.sWorkspace_output_basin, sFilename_out)
            export_flowline_to_geojson(aFlowline_conceptual_split, sFilename_out)
            #aFlowline_conceptual_split, dummy = read_flowline_geojson(sFilename_out)

        aFlowline_all = aFlowline_simplified_split + aFlowline_conceptual_split

        sFilename_output = self.sFilename_area_of_difference
        #remove headwater not needed here

        aPolygon_out, dArea = calculate_area_of_difference_simplified(aFlowline_all, aVertex_all, sFilename_output)
        print('Area of difference: ', dArea)
        self.dArea_of_difference = dArea
        self.dDistance_displace = dArea / self.dLength_flowline_simplified

        return
