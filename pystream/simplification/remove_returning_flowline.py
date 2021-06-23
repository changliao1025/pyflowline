from pystream.shared.hexagon import pyhexagon
def remove_returning_flowline(aHexagon_intersect_in, pVertex_outlet_in):
    aHexagon_out = list()
    #checking input data 

    #input hexagon should have at least one flowline inside

    #input flowline should have order information and segment index?

    
    nHexagon=  len(aHexagon_intersect_in)

    def checkIfDuplicates(listOfElems):
        ''' Check if given list contains any duplicates '''    
        for elem in listOfElems:
            if listOfElems.count(elem) > 1:
                return True
        return False


    def simplify_list(aHexagon_flowline_in):

        nHexagon = len(aHexagon_flowline_in)
        aHexagon_flowline_out = list()
        #check unique
        iFlag_unique = checkIfDuplicates(aHexagon_flowline_in)
        if iFlag_unique == True:
            return aHexagon_flowline_in
        else:
            for i in range(nHexagon):
                elem=aHexagon_flowline_in[i]
                if aHexagon_flowline_in.count(elem) > 1:
                    dummy = aHexagon_flowline_in.index(elem)
                    start = dummy[0]
                    end = dummy[-1]
                    del aHexagon_flowline_in[start, end-1] 
                    aHexagon_flowline_in = aHexagon_flowline_in
                    break
                    pass
                else:
                    pass

                pass 

            #check again
            iFlag_unique = checkIfDuplicates(aHexagon_flowline_in)
            if iFlag_unique == True:
                return aHexagon_flowline_in
            else:
                simplify_list(aHexagon_flowline_in)

            

    def retrieve_flowline_intersect_index(iSegment_in, lID_in, pVertex_end_in):

        iFlag_found = 1
        pVertex_end_current = pVertex_end_in
        aHexagon_flowline=list()
        #aHexagon_flowline.append(lID_in)
        aSegment_upstream = list()
        aVertex_end_upstream = list()
        while iFlag_found == 1:
            iFlag_found = 0 
            for j in range(nHexagon):
                pHexagon = aHexagon_intersect_in[j]
                lID = pHexagon.lID
                aFlowline= pHexagon.aFlowline
                nFlowline = len(aFlowline)
                for i in range(nFlowline):
                    pFlowline = aFlowline[i]
                    pVertex_start = pFlowline.pVertex_start
                    pVertex_end = pFlowline.pVertex_end
                    iSegment = pFlowline.iSegment
                    if pVertex_end == pVertex_end_current:
                        if iSegment == iSegment_in:
                            #found it
                            pVertex_end_current = pVertex_start

                            aHexagon_flowline.append(lID)
                            iFlag_found = 1
                            break
                            pass
                        else:
                            #this is a new upstream segment location
                            iFlag_found = 0
                            #save the upstream information for new step
                            aSegment_upstream.append(iSegment)
                            aVertex_end_upstream.append(pVertex_start)
                            break
                            pass

                    pass
        
        #should have finished search
        if iFlag_found == 0:

            #reverse 
            aHexagon_flowline = aHexagon_flowline[::-1]
            #simplify list
            aHexagon_simple = simplify_list(aHexagon_flowline)

            nUpstream = len(aSegment_upstream)
            for i in range(nUpstream):
                retrieve_flowline_intersect_index(aSegment_upstream[i], aVertex_end_upstream[i])
                pass
            pass

        return 

    #starting from the outlet
    
    for j in range(nHexagon):
        pHexagon = aHexagon_intersect_in[j]
        aFlowline= pHexagon.aFlowline
        nFlowline = len(aFlowline)
        for i in range(nFlowline):
            pFlowline = aFlowline[i]
            pVertex_start = pFlowline.pVertex_start
            pVertex_end = pFlowline.pVertex_end
            dDiatance = pVertex_end.calculate_distance( pVertex_outlet_in)
            if  dDiatance < 100.0:
                #found it
                nsegment = pFlowline.iSegment
                lID_outlet = pHexagon.lID
                pVertex_outlet = pVertex_end
                break
                pass    
            else:
                #print(dDiatance)
                pass

            pass     

    #loop through       
    iSegment = nsegment
    
    retrieve_flowline_intersect_index(iSegment, lID_outlet, pVertex_outlet)


    return aHexagon_out
