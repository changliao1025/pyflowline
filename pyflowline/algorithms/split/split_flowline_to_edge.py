from pyflowline.classes.flowline import pyflowline
def split_flowline_to_edge(aFlowline_in):
    aFlowline_out = list()
    aEdge_out=list()
    #nFlowline = len(aFlowline_in)
    #for i in range(nFlowline):
    #    pFlowline = aFlowline_in[i]     
    #    nEdge = pFlowline.nEdge     
    #    for j in range(nEdge):
    #        pEdge = pFlowline.aEdge[j]
    #        aEdge = list()
    #        aEdge.append(pEdge)
    #        aEdge_out.append(pEdge)            
    #        pFlowline1 = pyflowline(aEdge)      
    #        aFlowline_out.append(pFlowline1)    
    #        pass
    #    pass   
    for pFlowline in aFlowline_in:
        for pEdge in pFlowline.aEdge:
            aEdge_out.append(pEdge)
            aFlowline_out.append(pyflowline([pEdge]))
            
    return aFlowline_out, aEdge_out

