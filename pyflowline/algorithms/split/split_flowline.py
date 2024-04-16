
import numpy as np
from pyflowline.classes.edge import pyedge
from pyflowline.classes.flowline import pyflowline

import importlib.util
iFlag_cython = importlib.util.find_spec("cython") 
if iFlag_cython is not None:
    from pyflowline.algorithms.cython.kernel import find_vertex_on_edge
    from pyflowline.algorithms.cython.kernel import find_vertex_in_list
  
else:
    from pyflowline.algorithms.auxiliary.find_index_in_list import find_vertex_on_edge    
    from pyflowline.algorithms.auxiliary.find_vertex_in_list import find_vertex_in_list


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
            #get the start first
            #iFlag_exist, lIndex = find_vertex_in_list( aVertex_in,  pVertex)
            #if iFlag_exist == 1:                
            #    iPart = iPart + 1
            #    aVertex.append(pVertex)   
            #    pass

            if pVertex in aVertex_in_set:                
                iPart += 1
                aVertex.append(pVertex)

            if iFlag_intersect is not None:
                if iFlag_use_id is not None:
                    aDistance = list()
                    iFlag_exist=0
                    aVertex_dummy=list()
                    aIndex=list()
                    npoint =0
                    for k in range(len(aVertex_in)):
                        pVertex = aVertex_in[k]
                        if pVertex.lFlowlineID == pFlowline.lFlowlineID:
                            iFlag_exist =1
                            iPart = iPart + 1                            
                            distance  = pEdge.pVertex_start.calculate_distance(pVertex)
                            aDistance.append(distance)
                            aVertex_dummy.append(pVertex)
                            aIndex.append(k) 
                            npoint= npoint+ 1
                            
                    #sort needed
                    #if iFlag_exist == 1 :
                    #    x = np.array(aDistance)
                    #    b = np.argsort(x)
                    #    c = np.array(aIndex)
                    #    d= c[b]
                    #    aIndex_order = list(d)
                    #    #then push back
                    #    for k in range(npoint):
                    #        pVertex_dummy = aVertex_in[ aIndex_order[k]  ] 
                    #        aVertex.append(pVertex_dummy)
                    #        aVertex_all.append(pVertex_dummy)
                    #        pass
                    if aDistance:
                        aIndex_order = np.argsort(aDistance)
                        for k in aIndex_order:
                            pVertex_dummy = aVertex_in[k] 
                            aVertex.append(pVertex_dummy)
                            aVertex_all.append(pVertex_dummy)

                else:
                    iFlag_exist, npoint, aIndex = find_vertex_on_edge( aVertex_in, pEdge)
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
        #iFlag_exist, lIndex = find_vertex_in_list( aVertex_in,  pVertex)
        #if iFlag_exist == 1:
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
                    #rebuild index
                    #aVertex_index =list()
                    #for m in range(iPart):
                    #    pVertex= aVertex[m]
                    #    iFlag_exist, lIndex = find_vertex_in_list( aVertex_all,  pVertex)
                    #    if iFlag_exist ==1:
                    #        aVertex_index.append(lIndex)
                    #        pass
                    
                    aVertex_index = [aVertex_all.index(pVertex) for pVertex in aVertex if pVertex in aVertex_all]
          

                    #find duplicate
                    for k in range(nLine):
                        t = aVertex_index[k]
                        s = aVertex_index[k+1]
                        if s!=t:
                            #aEdge=list()
                            #for l in range(t,s):
                            #    pVertex0 = aVertex_all[l]  
                            #    pVertex1 = aVertex_all[l+1]  
                            #    pEdge = pyedge(pVertex0, pVertex1)
                            #    aEdge.append(pEdge)
                            #    pass
                            aEdge=[pyedge(aVertex_all[l], aVertex_all[l+1]) for l in range(t,s)]

                            pFlowline1 = pyflowline(aEdge)
                            pFlowline1.iStream_order = iStream_order
                            pFlowline1.iFlag_dam = iFlag_dam
                            aFlowline_out.append(pFlowline1)
                        pass                    
                    pass
                pass   

    return aFlowline_out

