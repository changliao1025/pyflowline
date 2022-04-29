

from pyflowline.algorithms.auxiliary.find_index_in_list import find_flowline_in_list
def remove_duplicate_flowline(aFlowline_in):
    aFlowline_out = list()
    nFlowline = len(aFlowline_in)
    for i in range(nFlowline):
        pFlowline = aFlowline_in[i]
        iFlag_exist, lIndex = find_flowline_in_list( aFlowline_out,  pFlowline)
        if iFlag_exist ==1:
            pass
        else:
            aFlowline_out.append(pFlowline)

        pass
    return aFlowline_out