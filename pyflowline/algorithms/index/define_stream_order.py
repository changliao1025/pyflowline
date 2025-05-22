
import numpy as np
import importlib.util
from pyflowline.algorithms.auxiliary.check_head_water import check_head_water #this function should not be used since stream order of headwater is available
iFlag_cython = importlib.util.find_spec("cython")
if iFlag_cython is not None:
    from tinyr import RTree
    iFlag_use_rtree = 1
else:
    iFlag_use_rtree =0
    pass

def update_head_water_stream_order(aFlowline_in):
    start_vertices = {flowline.pVertex_start for flowline in aFlowline_in}
    end_vertices = {flowline.pVertex_end for flowline in aFlowline_in}

    for flowline in aFlowline_in:
        pVertex_start = flowline.pVertex_start
        is_headwater = pVertex_start in start_vertices and pVertex_start not in end_vertices
        flowline.iStream_order = 1 if is_headwater else -1

    return aFlowline_in

    #nFlowline = len(aFlowline_in)
    #aFlowline_out = list()
    #for i in range(nFlowline):
    #    pFlowline = aFlowline_in[i]
    #    pVertex_start = pFlowline.pVertex_start
    #    if check_head_water(aFlowline_in, pVertex_start)==1:
    #        pFlowline.iStream_order = 1
    #    else:
    #        pFlowline.iStream_order = -1
    #
    #    aFlowline_out.append(pFlowline)
    #return aFlowline_out


def define_stream_order(aFlowline_in, aConfluence_in):
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
    if iFlag_use_rtree == 1:
        iMethod = 1
    else:
        iMethod = 2
        pass

    if nFlowline == 0 :
        print ('data incomplete')
    else:
        aStream_order = np.full(nFlowline, 0, dtype=int)
        if iMethod == 1: #the new method
            nConfleunce = len(aConfluence_in)
            aFlag_confluence_treated = np.full(nConfleunce, 0, dtype=int)
            #build rtree for confluence
            index_confluence = RTree( max_cap=5, min_cap=2)
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
                    iStream_segment = pFlowline_downstream.iStream_segment
                    #iFlag_upstream_done = 1
                    #nUpstream = len(aFlowline_upstream)
                    aStrord = [upstream.iStream_order for upstream in aFlowline_upstream if upstream.iStream_order >= 1]

                    #aStrord = list()
                    #for j in range(nUpstream):
                    #    pFlowline_upstream = aFlowline_upstream[j]
                    #    iStream_order_upstream = pFlowline_upstream.iStream_order
                    #    if iStream_order_upstream < 1:
                    #        iFlag_upstream_done = 0
                    #        break
                    #    else:
                    #        aStrord.append( iStream_order_upstream  )

                    #if iFlag_upstream_done == 1:
                    if len(aStrord) == len(aFlowline_upstream):
                        aFlag_confluence_treated[i] = 1
                        #now we can process the downstream
                        #get unique value
                        iStream_order = max(aStrord) if len(set(aStrord)) > 1 else aStrord[0] + 1

                        #update
                        pFlowline_downstream.iStream_order = iStream_order
                        aFlowline_in[nSegment-iStream_segment].iStream_order = iStream_order
                        #update confluence
                        #x = pFlowline_downstream.pVertex_end.dLongitude_degree
                        #y = pFlowline_downstream.pVertex_end.dLatitude_degree
                        #left =   x - 1E-5
                        #right =  x + 1E-5
                        #bottom = y - 1E-5
                        #top =    y + 1E-5
                        #pBound= (left, bottom, right, top)
                        #aIntersect = list(index_confluence.search(pBound))

                        #update confluence
                        x, y = pFlowline_downstream.pVertex_end.dLongitude_degree, pFlowline_downstream.pVertex_end.dLatitude_degree
                        #pBound = (x - 1E-5, y - 1E-5, x + 1E-5, y + 1E-5)
                        #aIntersect = list(index_confluence.search(pBound))
                        aIntersect = list(index_confluence.search_surrounding([x, y]))

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


        else: #the old method, not computationally efficient enough
            for i in range(nFlowline):
                pFlowline = aFlowline_in[i]
                pVertex_start=pFlowline.pVertex_start
                if check_head_water(aFlowline_in, pVertex_start)==1:
                    aStream_order[i] = 1
                    pass

            while aStream_order[0] == 0:
                for  i in range(nFlowline):
                    if aStream_order[i] !=0:
                        continue

                    pFlowline = aFlowline_in[i]
                    iFlag_upstream_done = 1
                    pVertex_start=pFlowline.pVertex_start
                    pVertex_end=pFlowline.pVertex_end
                    aStrord=list()
                    for  j in range(nFlowline):
                        pFlowline2 = aFlowline_in[j]
                        pVertex_start2=pFlowline2.pVertex_start
                        pVertex_end2=pFlowline2.pVertex_end
                        if pVertex_start == pVertex_end2:
                            if aStream_order[j] == 0:
                                iFlag_upstream_done=0
                                break
                            else:
                                aStrord.append( aStream_order[j]  )

                    if iFlag_upstream_done==1:
                        #get unique value
                        dummy = np.array(aStrord)
                        dummy1 = np.unique(dummy)

                        if len(dummy1) == 1: #all upstreams have the same order
                            aStream_order[i] = aStrord[0] + 1
                        else:
                            aStream_order[i] = np.max(dummy)

            for i in range(nFlowline):
                pFlowline = aFlowline_in[i]
                pFlowline.iStream_order =    aStream_order[i]
                aFlowline_out.append(pFlowline)

    return aFlowline_out, aStream_order