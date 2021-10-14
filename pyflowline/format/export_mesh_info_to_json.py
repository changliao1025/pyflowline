import numpy as np
import json

def export_mesh_info_to_json(aCell_in, aFlowline_in, sFilename_json_out):
    ncell=len(aCell_in)

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
                        if pVertex_center2 == pVertex_end:
                            aCell_in[k].lCellID_downstream_burned = aCell_in[l].lCellID
                            break
                    
                    

        pass

    with open(sFilename_json_out, 'w', encoding='utf-8') as f:
        sJson = json.dumps([json.loads(ob.tojson()) for ob in aCell_in], indent = 4)        
        f.write(sJson)    
        f.close()
    
    return
