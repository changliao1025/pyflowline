
import numpy as np
from rtree.index import Index as RTreeindex
from pyflowline.algorithms.auxiliary.check_head_water import check_head_water #this function should not be used since stream order of headwater is available


def update_head_water_stream_order(aFlowline_in):
    start_vertices = {flowline.pVertex_start for flowline in aFlowline_in}
    end_vertices = {flowline.pVertex_end for flowline in aFlowline_in}
    for flowline in aFlowline_in:
        pVertex_start = flowline.pVertex_start
        is_headwater = pVertex_start in start_vertices and pVertex_start not in end_vertices
        flowline.iStream_order = 1 if is_headwater else -1

    return aFlowline_in

def define_stream_order(aFlowline_in, aConfluence_in, iFlag_so_method_in=1):
    """define the stream order, but do we need to keep the confluence information?

    Args:
        aFlowline_in (_type_): _description_
        aConfluence_in (_type_): _description_

    Returns:
        _type_: _description_
    """
    nFlowline = len(aFlowline_in)
    nSegment = nFlowline
    aFlowline_out = list()


    if nFlowline == 0 :
        print ('data incomplete')
    else:
        aStream_order = np.full(nFlowline, 0, dtype=int)

        nConfleunce = len(aConfluence_in)
        aFlag_confluence_treated = np.full(nConfleunce, 0, dtype=int)
        #build rtree for confluence
        index_confluence = RTreeindex()
        for i, confluence in enumerate(aConfluence_in):
            pVertex_confluence = confluence.pVertex_confluence
            x, y = pVertex_confluence.dLongitude_degree, pVertex_confluence.dLatitude_degree
            pBound = (x - 1E-5, y - 1E-5, x + 1E-5, y + 1E-5)
            index_confluence.insert(i, pBound)
        while aFlowline_in[0].iStream_order < 0:
            for i in range(nConfleunce):
                if aFlag_confluence_treated[i] == 1:
                    continue
                pConfluence = aConfluence_in[i]
                aFlowline_upstream = pConfluence.aFlowline_upstream
                pFlowline_downstream = pConfluence.pFlowline_downstream
                if pFlowline_downstream is None:
                    continue
                iStream_segment = pFlowline_downstream.iStream_segment
                aStrord = [upstream.iStream_order for upstream in aFlowline_upstream if upstream.iStream_order >= 1]
                #if iFlag_upstream_done == 1:
                if len(aStrord) == len(aFlowline_upstream):
                    aFlag_confluence_treated[i] = 1
                    #get unique value
                    if iFlag_so_method_in == 1: # default method, strahler stream order
                        iStream_order = max(aStrord) if len(set(aStrord)) > 1 else aStrord[0] + 1
                    else:  #shreve stream order
                        iStream_order = sum(aStrord) if len(aStrord) > 0 else 1
                    #update
                    pFlowline_downstream.iStream_order = iStream_order
                    aFlowline_in[nSegment-iStream_segment].iStream_order = iStream_order
                    #update confluence
                    x, y = pFlowline_downstream.pVertex_end.dLongitude_degree, pFlowline_downstream.pVertex_end.dLatitude_degree
                    #aIntersect = list(index_confluence.search_surrounding([x, y]))
                    delta = 1E-5  # or your desired search tolerance
                    bbox = (x - delta, y - delta, x + delta, y + delta)
                    aIntersect = list(index_confluence.intersection(bbox))
                    for k in aIntersect:
                        pConfluence2 = aConfluence_in[k]
                        if pConfluence2.pVertex_confluence == pFlowline_downstream.pVertex_end:
                            for pFlowline_upstream2 in pConfluence2.aFlowline_upstream:
                                if pFlowline_upstream2.iStream_segment == iStream_segment:
                                    pFlowline_upstream2.iStream_order = iStream_order
                                    break
                        pass
        for i, flowline in enumerate(aFlowline_in):
            aFlowline_out.append(flowline)
            aStream_order[i] = flowline.iStream_order




    return aFlowline_out, aStream_order