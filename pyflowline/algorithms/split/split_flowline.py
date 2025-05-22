
import numpy as np
from pyflowline.classes.edge import pyedge
from pyflowline.classes.flowline import pyflowline

import importlib.util
iFlag_cython = importlib.util.find_spec("cython")
if iFlag_cython is not None:
    from pyflowline.algorithms.cython.kernel import find_vertex_on_edge
else:
    from pyflowline.algorithms.auxiliary.find_index_in_list import find_vertex_on_edge

def split_flowline(aFlowline_in, aVertex_in, iFlag_intersect = None, iFlag_use_id=None):
    """
    Split flowline based on the intersection with a list of vertex
    Input:
        aFlowline_in: list of flowline
        aVertex_in: list of vertex
        iFlag_intersect: 1: a vertex maybe on a line, but it is not a vertex of the line
        iFlag_use_id: 1:
    Output:
        aFlowline_out: list of flowline
    """
    aFlowline_out = list()
    nFlowline = len(aFlowline_in)
    aVertex_in_set = set(aVertex_in)

    for i in range(nFlowline):
        pFlowline = aFlowline_in[i]
        iStream_order = pFlowline.iStream_order
        #iStream_segment = pFlowline.iStream_segment
        iFlag_dam = pFlowline.iFlag_dam
        #nVertex = pFlowline.nVertex
        nEdge= pFlowline.nEdge
        iPart = 0
        aVertex  = list() #the actual vertex of ROI
        aVertex_all = list() #include vertex that is not ROI, but we need them to subset
        for j in range(nEdge):
            pEdge=pFlowline.aEdge[j]
            pVertex = pEdge.pVertex_start
            aVertex_all.append(pVertex)

            if pVertex in aVertex_in_set:
                iPart += 1
                aVertex.append(pVertex)

            if iFlag_intersect is not None:
                if iFlag_use_id is not None:
                    aDistance = list()
                    iFlag_exist=0
                    aVertex_dummy=list()
                    aIndex=list()
                    npoint = 0
                    for k in range(len(aVertex_in)):
                        pVertex0 = aVertex_in[k]
                        if pVertex0.lFlowlineID == pFlowline.lFlowlineID:
                            iFlag_exist =1
                            distance  = pEdge.pVertex_start.calculate_distance(pVertex0)
                            if distance == 0:
                                continue
                            iPart = iPart + 1
                            aDistance.append(distance)
                            aVertex_dummy.append(pVertex0)
                            aIndex.append(k)
                            npoint= npoint+ 1
                        else:
                            pass

                    if aDistance:
                        aIndex_order = np.argsort(aDistance)
                        for k in aIndex_order:
                            pVertex_dummy = aVertex_in[aIndex[k]]
                            aVertex.append(pVertex_dummy)
                            aVertex_all.append(pVertex_dummy)

                else:
                    iFlag_exist, npoint, aIndex = find_vertex_on_edge(aVertex_in, pEdge)
                    #they must be ordered
                    if iFlag_exist==1:
                        for m in range(npoint):
                            pVertex_dummy = aVertex_in[aIndex[m]]
                            iPart = iPart + 1
                            aVertex.append(pVertex_dummy)
                            aVertex_all.append(pVertex_dummy)

        #the last ending vertex
        pVertex = pFlowline.pVertex_end
        aVertex_all.append(pVertex)

        if pVertex in aVertex_in_set:
            iPart = iPart + 1
            aVertex.append(pVertex)
        if iPart == 0 :
            print('Something is wrong')
        else:
            if iPart ==1:
                #print('This flowline does not form any loop')
                if iFlag_use_id is not None:
                    pass
                pass
            else:
                if iPart >=2:
                    nLine = iPart-1
                    aVertex_index = [aVertex_all.index(pVertex) for pVertex in aVertex if pVertex in aVertex_all]

                    #find duplicate
                    for k in range(nLine):
                        t = aVertex_index[k]
                        s = aVertex_index[k+1]
                        if s!=t:
                            aEdge = [pyedge.create(aVertex_all[l], aVertex_all[l+1]) for l in range(t, s)]
                            # Remove None from the aEdge list
                            aEdge = [edge for edge in aEdge if edge is not None]
                            pFlowline1 = pyflowline(aEdge)
                            pFlowline1.iStream_order = iStream_order
                            pFlowline1.iFlag_dam = iFlag_dam
                            aFlowline_out.append(pFlowline1)
                        pass
                    pass
                pass

    return aFlowline_out

