def check_head_water(aFlowline_in, pVertex_start_in):
    """[Check whether a vertex assoacited with a flowline is a headwater or not]
    if the stream order info is avaialble, use the stream order info to check is easier and faster

    Args:
        aFlowline_in ([pyflowline]): [all the flowline]
        pVertex_start_in ([pyvertex]): [the vertex of interest]

    Returns:
        [int]: [0: not headwater; 1: is headwater]
    """    
    start_vertices = {flowline.pVertex_start for flowline in aFlowline_in}
    end_vertices = {flowline.pVertex_end for flowline in aFlowline_in}
    # Check if the vertex is a headwater
    is_headwater = pVertex_start_in in start_vertices and pVertex_start_in not in end_vertices
    return int(is_headwater)

    #nFlowline = len(aFlowline_in)
    #iFlag_head_water = -1
    #iCount = 0
    #for i in range(nFlowline):
    #    pFlowline = aFlowline_in[i]
    #    pVerter_start = pFlowline.pVertex_start
    #    pVerter_end = pFlowline.pVertex_end
    #    if pVerter_start == pVertex_start_in:
    #        iCount = iCount + 1
    #        pass
    #    if  pVerter_end == pVertex_start_in:
    #        iCount = iCount + 1
    #        pass
    #    pass
    #if iCount == 1:
    #    iFlag_head_water=1
    #    
    #return iFlag_head_water
    #Create sets of all start and end vertices

