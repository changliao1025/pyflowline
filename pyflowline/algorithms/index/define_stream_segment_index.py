from pyflowline.algorithms.auxiliary.check_head_water import check_head_water
def define_stream_segment_index(aFlowline_in):
    """build stream segment index, because they are ordered from outlet to headwater, so the index is from 1 to n

    Args:
        aFlowline_in (_type_): _description_

    Returns:
        _type_: _description_
    """
    nFlowline = len(aFlowline_in)
    aFlowline_out = list()
    aStream_segment = list()
    if nFlowline == 0 :
        print ('data incomplete')
    else:  
        for i in range(nFlowline):
            pFlowline = aFlowline_in[i]
            pFlowline.iStream_segment = nFlowline - i              
            aStream_segment.append(pFlowline.iStream_segment)
            aFlowline_out.append(pFlowline)

        pass
    return aFlowline_out, aStream_segment