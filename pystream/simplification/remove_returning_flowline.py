import copy
import numpy as np
from pystream.shared.hexagon import pyhexagon
from pystream.shared.flowline import pyflowline
from pystream.format.convert_coordinates_to_flowline import convert_coordinates_to_flowline
from pystream.add_unique_vertex import add_unique_vertex

from pystream.find_vertex_in_list import find_vertex_in_list
def remove_returning_flowline(aHexagon_in, aHexagon_intersect_in, pVertex_outlet_in):
    aHexagon_out = list()
    aFlowline_out=list()
    #checking input data 

    #input hexagon should have at least one flowline inside

    #input flowline should have order information and segment index?

    
    nHexagon=  len(aHexagon_intersect_in)

    def checkIfDuplicates(listOfElems):
        ''' Check if given list contains any duplicates '''    
        iFlag_unique = 1
        for elem in listOfElems:
            if listOfElems.count(elem) > 1:
                iFlag_unique = 0
                break
            else:
                pass
        
        return iFlag_unique


    def simplify_list(aHexagon_flowline_in):
        aHexagon_flowline_out = copy.deepcopy(aHexagon_flowline_in)

        nHexagon = len(aHexagon_flowline_out)        
        #check unique
        iFlag_unique = checkIfDuplicates(aHexagon_flowline_out)
        if iFlag_unique == 1:
            return aHexagon_flowline_out
        else:
            for i in range(nHexagon):
                elem=aHexagon_flowline_out[i]
                if aHexagon_flowline_out.count(elem) > 1:
                    dummy = [i for i, x in enumerate(aHexagon_flowline_out) if x == elem]
                    
                    start = dummy[0]
                    end = dummy[-1]
                    del aHexagon_flowline_out[start: end]                     
                    break
                    pass
                else:
                    pass

                pass             
            
            aHexagon_flowline_out = simplify_list(aHexagon_flowline_out)

        return   aHexagon_flowline_out

    def retrieve_flowline_intersect_index(iSegment_in, iStream_order_in, pVertex_end_in):

        iFlag_found = 1
        pVertex_end_current = pVertex_end_in
        aHexagon_flowline=list()
        #aHexagon_flowline.append(lID_in)
        aSegment_upstream = list()
        aVertex_end_upstream = list()
        aStream_order = list()
        while iFlag_found == 1:
            iFlag_found = 0 
            for j in range(nHexagon):
                pHexagon = aHexagon_intersect_in[j]
                lID = pHexagon.lIndex
                aFlowline= pHexagon.aFlowline
                nFlowline = len(aFlowline)
                for i in range(nFlowline):
                    pFlowline = aFlowline[i]
                    pVertex_start = pFlowline.pVertex_start
                    pVertex_end = pFlowline.pVertex_end
                    iSegment = pFlowline.iSegment
                    iStream_order = pFlowline.iStream_order
                    if pVertex_end == pVertex_end_current:
                        if iSegment == iSegment_in:
                            #found it
                            # #check length as well
                            dLength = pFlowline.dLength

                            iFlag_found2, dummy = pHexagon.which_edge_cross_this_vertex(pVertex_start)
                            iFlag_found3, dummy = pHexagon.which_edge_cross_this_vertex(pVertex_end)

                            if iFlag_found2 == 1 and iFlag_found3 == 1: #on the edge
                                if dLength < pHexagon.dLength :
                                    pVertex_end_current = pVertex_start  
                                    aHexagon_flowline.append(lID)
                                    iFlag_found = 1
                                    pass
                                else:
                                    pVertex_end_current = pVertex_start    
                                    aHexagon_flowline.append(lID)
                                    iFlag_found = 1
                                    pass

                                pass

                            else:
                                pVertex_end_current = pVertex_start    
                                aHexagon_flowline.append(lID)
                                iFlag_found = 1
                                pass                            
                                
                            break
                            pass
                        else:
                            iFlag_found2, dummy = pHexagon.which_edge_cross_this_vertex(pVertex_start)
                            iFlag_found3, dummy = pHexagon.which_edge_cross_this_vertex(pVertex_end)
                            if pVertex_end == pVertex_end_in: 
                                iFlag_found = 1
                                pass
                            else:
                                #this is a new upstream segment location
                                iFlag_found = 0
                                #save the upstream information for new step
                                aSegment_upstream.append(iSegment)
                                aVertex_end_upstream.append(pVertex_end)
                                aStream_order.append(iStream_order)
                               
                            
                            pass

                    pass

                if iFlag_found == 1:
                    break
        
        #should have finished search
        if iFlag_found == 0:

            #reverse 
            aHexagon_flowline = aHexagon_flowline[::-1]
            #simplify list
            #print(aHexagon_flowline)
            aHexagon_simple = simplify_list(aHexagon_flowline)
            #save the output
            nHexagon2 = len(aHexagon_simple)
            aCoordinates = list()
            if nHexagon2 >1:
                for i in range(nHexagon2):
                    lIndex = aHexagon_simple[i]
                    x = aHexagon_in[lIndex].dX_center
                    y = aHexagon_in[lIndex].dY_center
                    aCoordinates.append([x,y])
                    pass

                pFlowline = convert_coordinates_to_flowline(aCoordinates)
                pFlowline.iSegment = iSegment_in
                pFlowline.iStream_order=iStream_order_in
                aFlowline_out.append(pFlowline)

            print(iSegment_in, ': ',aHexagon_simple)

            sort_index = np.argsort(aStream_order)
            
            sort_index = sort_index[::-1]
            print(sort_index)
            nUpstream = len(aSegment_upstream)
            for i in sort_index:
                retrieve_flowline_intersect_index(aSegment_upstream[i], aStream_order[i], aVertex_end_upstream[i])
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
                iStream_order = pFlowline.iStream_order
                lID_outlet = pHexagon.lIndex
                pVertex_outlet = pVertex_end
                break
                pass    
            else:
                #print(dDiatance)
                pass

            pass     

    #loop through       
    iSegment = nsegment
        
    retrieve_flowline_intersect_index(iSegment, iStream_order,  pVertex_outlet)

    #remove parallel

    nsegment2 = len(aFlowline_out)
    aFlowline_out_no_parallel = list()
    aVertex_all = list()

    #for i in range(iStream_order,0,-1):
    #    for j in range(nsegment2):
    #        pFlowline = aFlowline_out[j]
    #        iord = pFlowline.iStream_order
    #        pVertex_start = pFlowline.pVertex_start
    #        pVertex_end = pFlowline.pVertex_end
    #        nVertex= pFlowline.nVertex
    #        aVertex = pFlowline.aVertex
    #        if iord == iStream_order:
    #            for k in aVertex:
    #                add_unique_vertex(aVertex_all, k)
    #                pass
#
    #            aFlowline_out_no_parallel.append(pFlowline)
    #            pass
    #        else:
    #            if iord==i:
    #                if nVertex > 2:
    #                    iFlag_exist = 0 
    #                    for m in range(1, nVertex-1,1):
    #                        pVertex = aVertex[m]
    #                        iFlag_exist, dummy =  find_vertex_in_list(aVertex_all, pVertex)
    #                        if iFlag_exist ==1:
    #                            #there is parallel river in this grid
    #                            #so we will abandon the whole flowline and its upstream?
    #                            break
    #                            pass
#
    #                        pass
#
    #                    if iFlag_exist ==1:                            
    #                        pass
    #                    else:
    #                        aFlowline_out_no_parallel.append(pFlowline)
    #                        for m in range( nVertex):
    #                            add_unique_vertex(aVertex_all, pVertex)
#
    #                else:
    #                    add_unique_vertex(aVertex_all, pVertex_start)
    #                    add_unique_vertex(aVertex_all, pVertex_end)
    #                    aFlowline_out_no_parallel.append(pFlowline)
    #                    pass
#
    #                pass
#
    #        pass
    #    
    #    pass
    



    return aHexagon_out, aFlowline_out, aFlowline_out_no_parallel
