import numpy as np
def remove_flowline_loop(aFlowline_in):
    aFlowline_out=list()
    nFlowline = len(aFlowline_in)
    def find_paralle_stream( pVertex_start_in):        
        ndownstream=0
        aDownstream=list()
        aStream_order_out=list()
        for j in range(nFlowline):
            pFlowline = aFlowline_in[j]
            pVertex_start = pFlowline.pVertex_start
            pVertex_end = pFlowline.pVertex_end
            if pVertex_start == pVertex_start_in: 
                ndownstream= ndownstream+1
                aDownstream.append(j)
                aStream_order_out.append(  pFlowline.iStream_order  )
                pass                
        return ndownstream, aDownstream, aStream_order_out

    lID=0
    aFlag = np.full(nFlowline, 0, dtype=int)
    for i in range(nFlowline):
        pFlowline = aFlowline_in[i]      
        pVertex_start = pFlowline.pVertex_start
        pVertex_end = pFlowline.pVertex_end
        iStream_order = pFlowline.iStream_order        
        ndownstream , aDownstream, aStream_order = find_paralle_stream( pVertex_start)
        if ndownstream == 1:
            if aFlag[i] !=1:
                pFlowline.lIndex = lID
                aFlowline_out.append(pFlowline)
                lID = lID + 1
                aFlag[i]=1
            pass
        else:                     
            #more than one, so we only take the current one or high order one            
            if(ndownstream>1):    
                sort_index = np.argsort(aStream_order)            
                sort_index = sort_index[::-1]
                lIndex = aDownstream[ sort_index[0] ]
                pFlowline_down =  aFlowline_in[lIndex]          
                if aFlag[ lIndex ]  !=1:
                    pFlowline_down.lIndex = lID
                    aFlowline_out.append(pFlowline_down)
                    lID = lID + 1
                    aFlag[lIndex]=1
                    pass
                
                for k in range(ndownstream):
                    aFlag[ aDownstream[k] ] = 1
            pass        
        pass 

    return aFlowline_out