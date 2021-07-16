from pystream.shared.flowline import pyflowline
def export_flowline_info_to_json(aCell_in, aFlowline_in, sFilename_json_out):
    
    #export the flowline topology to json

    ncell= len(aCell_in)
    nflowline= len(aFlowline_in)

    


    with open(sFilename_json_out, 'w', encoding='utf-8') as f:

        for i in range(1, nflowline+1):

            pFlowline = aFlowline_in[i]

            nVertex = pFlowline.nVertex
            nEdge = pFlowline.nEdge
            for j in range(1, nEdge+1):
                pEdge = pFlowline.aEdge[j]
                pVertex_start = pEdge.pVertex_start
                pVertex_end = pEdge.pVertex_end

                

            sJson = aCell_in[i-1].export_to_json()
            #print(sJson)
            f.write(sJson)

        f.close()
    
    return

    
   