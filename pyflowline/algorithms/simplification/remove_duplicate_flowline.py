
#from pyflowline.algorithms.auxiliary.find_index_in_list import find_flowline_in_list

def remove_duplicate_flowline(aFlowline_in):

    aFlowline_out = set(aFlowline_in)
    return list(aFlowline_out)
    #aFlowline_out = []
    #for pFlowline in aFlowline_in:
    #    iFlag_exist, _ = find_flowline_in_list(aFlowline_out, pFlowline)
    #    if iFlag_exist == 0:
    #        aFlowline_out.append(pFlowline)
    #return aFlowline_out