
import numpy as np
from pyflowline.classes.edge import pyedge
from pyflowline.classes.flowline import pyflowline
from pyflowline.algorithms.auxiliary.find_index_in_list import find_vertex_in_list, find_vertex_on_edge
def split_flowline(aFlowline_in, aVertex_in, iFlag_intersect = None, iFlag_use_id=None):
    aFlowline_out = list()
    nFlowline = len(aFlowline_in)
    for i in range(nFlowline):
        pFlowline = aFlowline_in[i]
        iStream_order = pFlowline.iStream_order
        iStream_segment = pFlowline.iStream_segment
        lFlowlineID = pFlowline.lFlowlineID        
        iFlag_dam = pFlowline.iFlag_dam
        nVertex = pFlowline.nVertex
        nEdge= pFlowline.nEdge
        iPart = 0        
        aVertex  = list() #the actual vertex of ROI
        aVertex_all = list() #include vertex that is not ROI, but we need them to subset 
        for j in range(nEdge):
            pEdge=pFlowline.aEdge[j]
            pVertex = pEdge.pVertex_start
            aVertex_all.append(pVertex)
            #get the start first
            iFlag_exist, lIndex = find_vertex_in_list( aVertex_in,  pVertex)
            if iFlag_exist == 1:                
                iPart = iPart + 1
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
                        if pVertex.lFlowlineID == lFlowlineID:
                            iFlag_exist =1
                            iPart = iPart + 1                            
                            distance  = pEdge.pVertex_start.calculate_distance(pVertex)
                            aDistance.append(distance)
                            aVertex_dummy.append(pVertex)
                            aIndex.append(k) 
                            npoint= npoint+ 1
                            
                    #sort needed
                    if iFlag_exist == 1 :
                        x = np.array(aDistance)
                        b = np.argsort(x)
                        c = np.array(aIndex)
                        d= c[b]
                        aIndex_order = list(d)
                        #then push back
                        for k in range(npoint):
                            pVertex_dummy = aVertex_in[ aIndex_order[k]  ] 
                            aVertex.append(pVertex_dummy)
                            aVertex_all.append(pVertex_dummy)

                else:
                    iFlag_exist, npoint, aIndex = find_vertex_on_edge( aVertex_in,  pEdge)
                    #they must be ordered
                    if iFlag_exist==1:
                        for m in range(npoint):
                            pVertex_dummy = aVertex_in[aIndex[m]]
                            iPart = iPart + 1
                            aVertex.append(pVertex_dummy)
                            aVertex_all.append(pVertex_dummy)
                        pass

        #the last ending vertex
        pVertex = pFlowline.pVertex_end
        aVertex_all.append(pVertex)
        iFlag_exist, lIndex = find_vertex_in_list( aVertex_in,  pVertex)
        if iFlag_exist == 1:
            iPart = iPart + 1
            aVertex.append(pVertex)
            pass

        if iPart == 0 :
            print('Something is wrong')
            pass
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
                    aVertex_index =list()
                    for m in range(iPart):
                        pVertex= aVertex[m]
                        iFlag_exist, lIndex = find_vertex_in_list( aVertex_all,  pVertex)
                        if iFlag_exist ==1:
                            aVertex_index.append(lIndex)
                    
                    #find duplicate
                    for k in range(nLine):
                        t = aVertex_index[k]
                        s = aVertex_index[k+1]
                        if s!=t:
                            aEdge=list()
                            for l in range(t,s):
                                pVertex0 = aVertex_all[l]  
                                pVertex1 = aVertex_all[l+1]  
                                pEdge = pyedge(pVertex0, pVertex1)
                                aEdge.append(pEdge)
                                pass

                            pFlowline1 = pyflowline(aEdge)
                            pFlowline1.iStream_order = iStream_order
                            pFlowline1.iFlag_dam = iFlag_dam
                            aFlowline_out.append(pFlowline1)
                        pass                    
                    pass
                pass   

    return aFlowline_out

