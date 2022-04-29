import os
import numpy as np
import json

def export_mesh_info_to_json(iFlag_flowline_in, aCell_in, aFlowline_in, aCellID_outlet_iin, sFilename_json_in):
    """_summary_

    Args:
        iFlag_flowline_in (_type_): _description_
        aCell_in (_type_): _description_
        aFlowline_in (_type_): _description_
        aCellID_outlet_iin (_type_): _description_
        sFilename_json_in (_type_): _description_
    """
    if os.path.exists(sFilename_json_in): 
        os.remove(sFilename_json_in)
        pass
    ncell=len(aCell_in)
    if iFlag_flowline_in == 1:
        nFlowline = len(aFlowline_in)
        for i in range(nFlowline):
            pFlowline = aFlowline_in[i]
            nEdge = pFlowline.nEdge
            nVertex = pFlowline.nVertex
            aEdge = pFlowline.aEdge
            iStream_segment = pFlowline.iStream_segment
            iStream_order = pFlowline.iStream_order
            for j in range(nEdge):
                pEdge = aEdge[j]
                pVertex_start = pEdge.pVertex_start
                pVertex_end = pEdge.pVertex_end
                for k in range(ncell):
                    pVertex_center = aCell_in[k].pVertex_center
                    if pVertex_center == pVertex_start:
                        aCell_in[k].iStream_segment_burned = iStream_segment
                        aCell_in[k].iStream_order_burned = iStream_order
                        for l in range(ncell):
                            pVertex_center2 = aCell_in[l].pVertex_center
                            lCellID = aCell_in[l].lCellID
                            if pVertex_center2 == pVertex_end:
                                aCell_in[k].lCellID_downstream_burned = lCellID
                                if lCellID in aCellID_outlet_iin:
                                    aCell_in[l].iStream_segment_burned = iStream_segment
                                    aCell_in[l].iStream_order_burned = iStream_order
                                break
                            
                            

            pass
    else:
        #only mesh, no flowline
        pass

    with open(sFilename_json_in, 'w', encoding='utf-8') as f:
        sJson = json.dumps([json.loads(ob.tojson()) for ob in aCell_in], indent = 4)        
        f.write(sJson)    
        f.close()
    
    return
