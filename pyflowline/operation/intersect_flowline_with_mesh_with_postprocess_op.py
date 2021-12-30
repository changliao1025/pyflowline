import os
import json
from pathlib import Path
from pyflowline.classes.vertex import pyvertex
from pyflowline.algorithms.auxiliary.reproject_coordinates import reproject_coordinates
from pyflowline.format.read_flowline_shapefile import read_flowline_shapefile
#from pyflowline.format.read_mesh_shapefile import read_mesh_shapefile
from pyflowline.format.read_mesh_json import read_mesh_json
from pyflowline.format.read_flowline_geojson import read_flowline_geojson

#from pyflowline.format.export_flowline_to_shapefile import export_flowline_to_shapefile
from pyflowline.format.export_vertex_to_json import export_vertex_to_json
from pyflowline.format.export_flowline_to_json import export_flowline_to_json
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
from pyflowline.format.export_vertex_to_shapefile import export_vertex_to_shapefile
from pyflowline.algorithms.merge.merge_flowline import merge_flowline

from pyflowline.algorithms.index.define_stream_order import define_stream_order
from pyflowline.algorithms.index.define_stream_segment_index import define_stream_segment_index

def intersect_flowline_with_mesh_with_postprocess_op(oPyflowline_in):

    #important
    

    iMesh_type = oPyflowline_in.iMesh_type
    iFlag_intersect = oPyflowline_in.iFlag_intersect
    sWorkspace_output = oPyflowline_in.sWorkspace_output
    nOutlet = oPyflowline_in.nOutlet

    sFilename_mesh=oPyflowline_in.sFilename_mesh
    
    aMesh, pSpatial_reference_mesh = read_mesh_json(sFilename_mesh)

    
    aCell = list()
    aCell_intersect = list()
    aFlowline = list()   #store all the flowline
    aOutletID = list()
    aBasin = list()
    if iFlag_intersect == 1:
        for i in range(0, nOutlet, 1):
            sBasin =  "{:03d}".format(i+1)    
            print('Flowline interset with postprocess ',sBasin)         
            
            pBasin = oPyflowline_in.aBasin[i]
            sWorkspace_output_basin = pBasin.sWorkspace_output_basin
            sFilename_flowline = pBasin.sFilename_flowline_segment_order_before_intersect
            sFilename_flowline_in = os.path.join(sWorkspace_output_basin, sFilename_flowline)
            sFilename_flowline_intersect = pBasin.sFilename_flowline_intersect
            sFilename_flowline_intersect_out = os.path.join(sWorkspace_output_basin, sFilename_flowline_intersect)

            aCell, aCell_intersect_basin, aFlowline_intersect_all = intersect_flowline_with_mesh(iMesh_type, sFilename_mesh, \
                sFilename_flowline_in, sFilename_flowline_intersect_out)

            sFilename_flowline_filter = pBasin.sFilename_flowline_filter

            aFlowline_basin, pSpatialRef_flowline = read_flowline_shapefile(sFilename_flowline_filter)

            iFlag_projected = 0

            point= dict()

            point['lon'] = pBasin.dLongitude_outlet_degree
            point['lat'] = pBasin.dLatitude_outlet_degree
            pVertex_outlet_initial=pyvertex(point)

            #from this point, aFlowline_basin is conceptual
            aFlowline_basin, aFlowline_no_parallel, lCellID_outlet, pVertex_outlet \
                = remove_returning_flowline(iMesh_type, aCell_intersect_basin, pVertex_outlet_initial)
            sFilename_out = 'flowline_simplified_after_intersect_' + sBasin + '.json'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)  

            pSpatial_reference =  pSpatial_reference_mesh

            export_flowline_to_json(iFlag_projected, aFlowline_basin, pSpatial_reference, sFilename_out)

            #added start
            aFlowline_basin, aEdge = split_flowline_to_edge(aFlowline_basin)

            aFlowline_basin = remove_duplicate_flowline(aFlowline_basin)
            aFlowline_basin = correct_flowline_direction(aFlowline_basin,  pVertex_outlet )
            aFlowline_basin = remove_flowline_loop(  aFlowline_basin )  

            sFilename_out = 'flowline_debug_' + sBasin + '.json'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_flowline_to_json(iFlag_projected, aFlowline_basin, pSpatial_reference, sFilename_out)

            aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity\
                = find_flowline_confluence(aFlowline_basin,  pVertex_outlet)

            sFilename_out = 'flowline_vertex_with_confluence_01_after_intersect_' + sBasin + '.json'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_vertex_to_json(iFlag_projected, aVertex, pSpatial_reference, sFilename_out, aAttribute_data=aConnectivity)


            aFlowline_basin = merge_flowline( aFlowline_basin,aVertex, pVertex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence  )  

            aFlowline_basin = remove_flowline_loop(  aFlowline_basin )    

            aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity\
                = find_flowline_confluence(aFlowline_basin,  pVertex_outlet)

            aFlowline_basin = merge_flowline( aFlowline_basin,aVertex, pVertex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence  ) 

            aFlowline_basin, aStream_segment = define_stream_segment_index(aFlowline_basin)
            aFlowline_basin, aStream_order = define_stream_order(aFlowline_basin)

            sFilename_out = pBasin.sFilename_flowline_final
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_flowline_to_json(iFlag_projected, aFlowline_basin, pSpatial_reference, sFilename_out)

            aFlowline = aFlowline + aFlowline_basin
            
            aCell_intersect = aCell_intersect + aCell_intersect_basin
            aOutletID.append(lCellID_outlet)

            pBasin.lCellID_outlet = lCellID_outlet
            pBasin.dLongitude_outlet_degree = pVertex_outlet.dLongitude_degree
            pBasin.dLatitude_outlet_degree = pVertex_outlet.dLatitude_degree

            aBasin.append(pBasin)

        print('Finished flowline intersect')
    else:
        print('Flowline intersect was skiped')
    
    #save basin info
    sPath = os.path.dirname(oPyflowline_in.sFilename_basins)
    sName = str(Path(oPyflowline_in.sFilename_basins).stem ) + '_new.json'
    sFilename_configuration  =  os.path.join( sPath  , sName)
    with open(sFilename_configuration, 'w', encoding='utf-8') as f:
        sJson = json.dumps([json.loads(ob.tojson()) for ob in aBasin],\
            sort_keys=True, \
            indent = 4)        
            
        f.write(sJson)    
        f.close()

    return aCell, aCell_intersect, aFlowline, aOutletID

