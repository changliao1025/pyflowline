

from pyflowline.algorithms.auxiliary.find_index_in_list import find_edge_in_list
def remove_duplicate_edge(aEdge_in):
    aEdge_out = list()
    nEdge = len(aEdge_in)
    for i in range(nEdge):
        pEdge = aEdge_in[i]
        iFlag_exist, lIndex = find_edge_in_list( aEdge_out,  pEdge)
        if iFlag_exist ==1:
            pass
        else:
            aEdge_out.append(pEdge)

        pass
    return aEdge_out