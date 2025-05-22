import numpy as np



def split_flowline_by_length(aFlowline_in, dDistance):

    aFlowline_out=list()
    nFlowline = len(aFlowline_in)
    for i in range(nFlowline):
        pFlowline = aFlowline_in[i]
        pFlowline_out = pFlowline.split_edge_by_length(dDistance)
        aFlowline_out.append(pFlowline_out)

    return aFlowline_out

def split_edge_by_length(pEdge_in, dLength_in):
    from pyflowline.classes.edge import pyedge
    aEdge_out=list()

    dLength_total = pEdge_in.dLength

    if dLength_total<= dLength_in:
        aEdge_out.append(pEdge_in)
    else:
        pVertex_start = pEdge_in.pVertex_start
        pVertex_end = pEdge_in.pVertex_end
        #nBreak = 1
        #for i in range(1,10):
        #    if np.power(2, i) * dLength_in > dLength_total:
        #        nBreak = i
        #        break
        #    else:
        #        pass
        ##now we have the break count
        #nPoint = np.power(2, nBreak) + 1
        #nEdge = nPoint-1

        n1 = pVertex_start.toNvector()
        n2 = pVertex_end.toNvector()

        mid = n1.plus(n2)

        pVertex_mid = mid.toLatLon()

        pEdge1 = pyedge(pVertex_start, pVertex_mid)

        pEdge2 = pyedge(pVertex_mid, pVertex_end)
        if pEdge1.dLength<=dLength_in:
            aEdge_out.append(pEdge1)
        else:
            aEdge_out_dummy1 = split_edge_by_length(pEdge1, dLength_in)
            for pe in aEdge_out_dummy1:
                aEdge_out.append(pe)



        if pEdge2.dLength<=dLength_in:
            aEdge_out.append(pEdge2)
        else:
            aEdge_out_dummy2 = split_edge_by_length(pEdge2, dLength_in)
            for pe in aEdge_out_dummy2:
                aEdge_out.append(pe)


    return aEdge_out

