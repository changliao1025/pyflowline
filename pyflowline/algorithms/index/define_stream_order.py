
import numpy as np 
import importlib
from pyflowline.algorithms.auxiliary.check_head_water import check_head_water #this function should not be used since stream order of headwater is available 
iFlag_cython = importlib.util.find_spec("cython") 
if iFlag_cython is not None:
    from pyflowline.external.tinyr.tinyr.tinyr import RTree
    iFlag_use_rtree = 1
else:
    iFlag_use_rtree =0
    pass

def update_head_water_stream_order(aFlowline_in):
    nFlowline = len(aFlowline_in)
    aFlowline_out = list()
    for i in range(nFlowline):
        pFlowline = aFlowline_in[i]          
        pVertex_start = pFlowline.pVertex_start  
        if check_head_water(aFlowline_in, pVertex_start)==1:
            pFlowline.iStream_order = 1
            pass
        else:
            pFlowline.iStream_order = -1
        
        aFlowline_out.append(pFlowline)

    return aFlowline_out


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
    iMethod = 1       
    if nFlowline == 0 :
        print ('data incomplete')
    else:       
        aStream_order = np.full(nFlowline, 0, dtype=int)  
        if iMethod == 1: #the new method
            nConfleunce = len(aConfluence_in)
            aFlag_confluence_treated = np.full(nConfleunce, 0, dtype=int)
            #build rtree for confluence
            index_confluence = RTree( max_cap=5, min_cap=2)
            for i in range(nConfleunce):
                lID = i 
                pVertex_confluence = aConfluence_in[i].pVertex_confluence 
                x = pVertex_confluence.dLongitude_degree
                y = pVertex_confluence.dLatitude_degree   
                left =   x - 1E-5
                right =  x + 1E-5
                bottom = y - 1E-5
                top =    y + 1E-5
                pBound= (left, bottom, right, top)                
                index_confluence.insert(lID, pBound)  #

          
            while aFlowline_in[0].iStream_order < 0:
                for i in range(nConfleunce):
                    if aFlag_confluence_treated[i] == 1:
                        continue

                    pConfluence = aConfluence_in[i]
                    aFlowline_upstream = pConfluence.aFlowline_upstream 
                    pFlowline_downstream = pConfluence.pFlowline_downstream
                    iStream_segment = pFlowline_downstream.iStream_segment
                    iFlag_upstream_done = 1
                    nUpstream = len(aFlowline_upstream)
                    aStrord = list()
                    for j in range(nUpstream):
                        pFlowline_upstream = aFlowline_upstream[j]
                        iStream_order_upstream = pFlowline_upstream.iStream_order
                        if iStream_order_upstream < 1:
                            iFlag_upstream_done = 0
                            break
                        else:
                            aStrord.append( iStream_order_upstream  )
                        pass

                    if iFlag_upstream_done == 1:
                        aFlag_confluence_treated[i] = 1
                        #now we can process the downstream
                        #get unique value
                        dummy = np.array(aStrord)
                        dummy1 = np.unique(dummy)
                        if len(dummy1) == 1: #all upstreams have the same order
                            iStream_order = aStrord[0] + 1
                        else:
                            iStream_order = np.max(dummy)  
                        pass
                        
                        #update
                        pFlowline_downstream.iStream_order = iStream_order
                        aFlowline_in[nSegment-iStream_segment].iStream_order = iStream_order
                        #update confluence
                        x = pFlowline_downstream.pVertex_end.dLongitude_degree
                        y = pFlowline_downstream.pVertex_end.dLatitude_degree   
                        left =   x - 1E-5
                        right =  x + 1E-5
                        bottom = y - 1E-5
                        top =    y + 1E-5
                        pBound= (left, bottom, right, top)  
                        aIntersect = list(index_confluence.search(pBound))
                        for k in aIntersect:
                            pConfluence2 = aConfluence_in[k]
                            if pConfluence2.pVertex_confluence == pFlowline_downstream.pVertex_end:
                                for pFlowline_upstream2 in pConfluence2.aFlowline_upstream:
                                    if pFlowline_upstream2.iStream_segment == iStream_segment:
                                        pFlowline_upstream2.iStream_order = iStream_order
                                        break
                            pass

                                             
            for i in range(nFlowline):
                pFlowline = aFlowline_in[i]      
                pFlowline.iStream_order =    aFlowline_in[i].iStream_order   
                aFlowline_out.append(pFlowline)
                aStream_order[i] = pFlowline.iStream_order
                
            pass
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