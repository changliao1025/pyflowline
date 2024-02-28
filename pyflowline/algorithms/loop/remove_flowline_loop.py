import numpy as np
def remove_flowline_loop(aFlowline_in):
    """_summary_

    Args:
        aFlowline_in (_type_): _description_

    Returns:
        List: List of flowline, ordered from outlet to headwater
    """
    aFlowline_out=list()
    nFlowline = len(aFlowline_in)

    # Create a dictionary that maps pVertex_start to a list of flowlines
    flowline_dict = {}
    for i in range(nFlowline):
        pFlowline = aFlowline_in[i]
        pVertex_start = pFlowline.pVertex_start
        if pVertex_start not in flowline_dict:
            flowline_dict[pVertex_start] = []
        flowline_dict[pVertex_start].append((i, pFlowline))

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
        #ndownstream , aDownstream, aStream_order = find_paralle_stream( pVertex_start)

        # Get all parallel streams
        parallel_streams = flowline_dict.get(pVertex_start, [])
        ndownstream = len(parallel_streams)
        aDownstream = [j for j, _ in parallel_streams]
        aStream_order = [flowline.iStream_order for _, flowline in parallel_streams]

        if ndownstream == 1:
            if aFlag[i] !=1:
                pFlowline.lFlowlineIndex = lID
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
                    pFlowline_down.lFlowlineIndex = lID
                    aFlowline_out.append(pFlowline_down)
                    lID = lID + 1
                    aFlag[lIndex]=1
                    pass
                
                for k in range(ndownstream):
                    aFlag[ aDownstream[k] ] = 1
           

    return aFlowline_out