def find_hexagon_through_edge(aHexagon, pEdge_in):

    nHexagon = len(aHexagon)
    aHexagon_out = list()
    for i in range(nHexagon):
        pHexagon = aHexagon[i]
        if pHexagon.has_this_edge(pEdge_in) ==1:
            aHexagon_out.append(pHexagon)
            pass
        else:
            pass

    return aHexagon_out