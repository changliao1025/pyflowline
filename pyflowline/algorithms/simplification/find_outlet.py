def find_outlet(aFlowline_in, pVertex_outlet_in):
    nFlowline = len(aFlowline_in)
    dDistance_min = float("inf")
    for i in range(nFlowline):
        pFlowline = aFlowline_in[i]
        pVertex_end = pFlowline.pVertex_end
        dDistance = pVertex_end.calculate_distance(pVertex_outlet_in)
        if dDistance < dDistance_min:
            dDistance_min = dDistance
            lIndex_outlet = i

    pVertex_outlet = aFlowline_in[lIndex_outlet].pVertex_end
    return pVertex_outlet
