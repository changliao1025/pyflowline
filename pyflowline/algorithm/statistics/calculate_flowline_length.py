
from pyflowline.classes.flowline import pyflowline
def calculate_flowline_length(aFlowline_in):

    dLength = 0.0

    nflowline = len(aFlowline_in)

    for i in range(nflowline):

        pFlowline= aFlowline_in[i]

        pFlowline.calculate_length()

        dLength = dLength + pFlowline.dLength
    

    return dLength