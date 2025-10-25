from pyflowline.classes.flowline import pyflowline
def split_flowline_to_edge(aFlowline_in):
    aFlowline_out = list()
    aEdge_out=list()
    for pFlowline in aFlowline_in:
        for pEdge in pFlowline.aEdge:
            aEdge_out.append(pEdge)
            aFlowline_out.append(pyflowline([pEdge]))

    return aFlowline_out, aEdge_out

