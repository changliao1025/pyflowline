from pyflowline.classes.flowline import pyflowline
def split_flowline_by_length(aFlowline_in, dDistance):

    aFlowline_out=list()
    nFlowline = len(aFlowline_in)
    for i in range(nFlowline):
        pFlowline = aFlowline_in[i]
        pFlowline_out = pFlowline.split_by_length(dDistance)
            
        aFlowline_out.append(pFlowline_out)    
          

    return aFlowline_out

def split_edge_by_length(pEdge_in, dLength_in):
    aEdge_out=list()
    return aEdge_out

