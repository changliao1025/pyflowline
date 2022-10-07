from pyflowline.classes.flowline import pyflowline
def break_flowline_by_length(aFlowline_in, dDistance):

    aFlowline_out=list()
    nFlowline = len(aFlowline_in)
    for i in range(nFlowline):
        pFlowline = aFlowline_in[i]
        pFlowline_out = pFlowline.break_by_length(dDistance)
            
        aFlowline_out.append(pFlowline_out)    
          

    return aFlowline_out
