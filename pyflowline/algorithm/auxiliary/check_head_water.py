def check_head_water(aFlowline_in, pVertex_start_in):
    nFlowline = len(aFlowline_in)
    iFlag_head_water = -1
    iCount = 0
    for i in range(nFlowline):
        pFlowline = aFlowline_in[i]
        pVerter_start = pFlowline.pVertex_start
        pVerter_end = pFlowline.pVertex_end
        if pVerter_start == pVertex_start_in:
            iCount = iCount + 1
            pass
        if  pVerter_end == pVertex_start_in:
            iCount = iCount + 1
            pass
        pass
    if iCount == 1:
        iFlag_head_water=1
        
    return iFlag_head_water