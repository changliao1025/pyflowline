def remove_small_river(aFlowline_in, dThreshold_in):
    """Remove small river that meet the threshold and headwater requirement, also dam flowline are reserved

    Args:
        aFlowline_in (_type_): _description_
        dThreshold_in (_type_): _description_

    Returns:
        List: Theortically, the flowline should be ordered from outlet to headwater
    """  
         
    if len(aFlowline_in) == 1:
        return [aFlowline_in[0]]
    else:
        aFlowline_out=list()   
        lID = 0        
        for pFlowline in aFlowline_in:            
            iFlag_dam = pFlowline.iFlag_dam           
            dLength = pFlowline.calculate_length()
            if iFlag_dam == 1 or pFlowline.iStream_order != 1 or dLength > dThreshold_in:
                pFlowline.lFlowlineIndex = lID
                aFlowline_out.append(pFlowline)
                lID += 1
            #if iFlag_dam ==1:
            #    pFlowline.lFlowlineIndex = lID
            #    aFlowline_out.append(pFlowline)
            #    lID = lID +1       
            #else:
            #    if pFlowline.iStream_order == 1:
            #        if dLength > dThreshold_in :
            #            pFlowline.lFlowlineIndex = lID
            #            aFlowline_out.append(pFlowline)
            #            lID = lID + 1                        
            #    else: #this one might be not used     
            #        pFlowline.lFlowlineIndex = lID
            #        aFlowline_out.append(pFlowline)
            #        lID = lID +1       
            #        pass        
                    #if check_head_water(aFlowline_in, pVertex_start)==1:
                    #    if dLength > dThreshold_in :
                    #        pFlowline.lIndex = lID
                    #        aFlowline_out.append(pFlowline)
                    #        lID = lID + 1 
                    #        pass
                    #    else:
                    #        pass
                    #else:        
                    #    pFlowline.lIndex = lID
                    #    aFlowline_out.append(pFlowline)
                    #    lID = lID +1       
                    #    pass            
            
           

    return aFlowline_out