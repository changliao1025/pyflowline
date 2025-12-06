def define_stream_topology(aFlowline_in, aConfluence_in):
    """build the upstream -downstream relationship

    Args:
        aFlowline_in (_type_): It has segment information
        aConfluence_in (_type_): It has no segment information

    Returns:
        _type_: _description_
    """
    nFlowline = len(aFlowline_in)
    

    if nFlowline == 0 :
        print ('data incomplete')
        return []

    for pConfluence in aConfluence_in:
        aFlowline_upstream = pConfluence.aFlowline_upstream
        pFlowline_downstream = pConfluence.pFlowline_downstream
        if pFlowline_downstream is None:
            continue
        else:
            iStream_segment_downstream = pFlowline_downstream.iStream_segment
            lFlowline_index_downstream = nFlowline - iStream_segment_downstream
            aFlowline_in[lFlowline_index_downstream].aFlowline_upstream = []

            #process upstream
            for pFlowline_upstream in aFlowline_upstream:
                iStream_segment_upstream = pFlowline_upstream.iStream_segment
                lFlowline_index_upstream = nFlowline - iStream_segment_upstream
                aFlowline_in[lFlowline_index_upstream].lFlowlineIndex_downstream = lFlowline_index_downstream
                aFlowline_in[lFlowline_index_downstream].aFlowline_upstream.append(iStream_segment_upstream)

    return aFlowline_in

    #return aFlowline_out