import sys
import numpy as np
from pyflowline.algorithms.auxiliary.check_head_water import check_head_water
sys.setrecursionlimit(100000)
lFlowlineIndex=0

def correct_flowline_direction(aFlowline_in, pVertex_outlet_in):
    nFlowline = len(aFlowline_in)
    dDiatance_min = float('inf')
    unfinished_flowlines = set()
    vertex_to_flowlines = {}  # New dictionary to map each vertex to its flowlines

    for i in range(nFlowline):
        pFlowline = aFlowline_in[i]
        pVertex_end = pFlowline.pVertex_end
        dDiatance = pVertex_end.calculate_distance(pVertex_outlet_in)
        if dDiatance < dDiatance_min:
            dDiatance_min = dDiatance
            lIndex_outlet = i

    for i in range(nFlowline):
        if i != lIndex_outlet:
            pFlowline = aFlowline_in[i]
            unfinished_flowlines.add(aFlowline_in[i])
            # Update the dictionary
            vertex_to_flowlines.setdefault(pFlowline.pVertex_start, []).append(pFlowline)
            vertex_to_flowlines.setdefault(pFlowline.pVertex_end, []).append(pFlowline)
            #print(pFlowline.pVertex_start.wkt)
            #print(pFlowline.pVertex_end.wkt)

    aVertex_downslope_table = [aFlowline_in[lIndex_outlet].pVertex_start]
    aFlowline_out= [aFlowline_in[lIndex_outlet]]
    while unfinished_flowlines:
        aVertex_downslope_current= []
        iCount = 0
        for pVertex_dummy in aVertex_downslope_table:
            to_remove = set()
            #for pFlowline in unfinished_flowlines:
            aFlowline_dummy = vertex_to_flowlines.get(pVertex_dummy, [])
            for pFlowline in aFlowline_dummy:  # Use the dictionary here
                if pFlowline in unfinished_flowlines:
                    if pFlowline.pVertex_end == pVertex_dummy :
                        aVertex_downslope_current.append(pFlowline.pVertex_start)
                        to_remove.add(pFlowline)
                        aFlowline_out.append(pFlowline)
                        iCount = iCount + 1
                    else:
                        if pFlowline.pVertex_start == pVertex_dummy :
                            pFlowline.reverse()
                            aVertex_downslope_current.append(pFlowline.pVertex_start)
                            to_remove.add(pFlowline)
                            aFlowline_out.append(pFlowline)
                            iCount = iCount + 1

            unfinished_flowlines -= to_remove

        if iCount == 0:
           break

        if len(unfinished_flowlines)==0:
           break
        aVertex_downslope_table = aVertex_downslope_current

    return  aFlowline_out


