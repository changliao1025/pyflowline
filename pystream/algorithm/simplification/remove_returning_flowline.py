import copy
import numpy as np
from pystream.shared.hexagon import pyhexagon
from pystream.shared.flowline import pyflowline
from pystream.format.convert_coordinates_to_flowline import convert_coordinates_to_flowline

from pystream.algorithm.auxiliary.find_vertex_in_list import find_vertex_in_list
def remove_returning_flowline(iMesh_type, aCell_in, aCell_intersect_in, pVertex_outlet_in):
    aCell_out = list()
    aFlowline_out=list()
    #checking input data 

    #input hexagon should have at least one flowline inside

    #input flowline should have order information and segment index?

    
    nCell=  len(aCell_intersect_in)

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


    def simplify_list(aCell_flowline_in):
        aCell_flowline_out = copy.deepcopy(aCell_flowline_in)

        nCell2 = len(aCell_flowline_out)        
        #check unique
        iFlag_unique = checkIfDuplicates(aCell_flowline_out)
        if iFlag_unique == 1:
            return aCell_flowline_out
        else:
            for i in range(nCell2):
                elem=aCell_flowline_out[i]
                if aCell_flowline_out.count(elem) > 1:
                    dummy = [i for i, x in enumerate(aCell_flowline_out) if x == elem]
                    
                    start = dummy[0]
                    end = dummy[-1]
                    del aCell_flowline_out[start: end]                     
                    break
                    pass
                else:
                    pass

                pass             
            
            aCell_flowline_out = simplify_list(aCell_flowline_out)

        return   aCell_flowline_out

    def retrieve_flowline_intersect_index(iSegment_in, iStream_order_in, pVertex_end_in):

        iFlag_found = 1
        pVertex_end_current = pVertex_end_in
        aCell_flowline=list()
  
        aSegment_upstream = list()
        aVertex_end_upstream = list()
        aStream_order = list()
        iFlag_skip = 0
        iFlag_previous_overlap=0
        while iFlag_found == 1:
            iFlag_found = 0 
            
            for j in range(nCell):
                pCell = aCell_intersect_in[j]
                lID = pCell.lIndex
                aFlowline= pCell.aFlowline
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

                            iFlag_found2, dummy2 = pCell.which_edge_cross_this_vertex(pVertex_start)
                            iFlag_found3, dummy3 = pCell.which_edge_cross_this_vertex(pVertex_end)

                            if iFlag_found2 == 1 and iFlag_found3 == 1: #on the edge
                                #same edge
                                if dummy2.is_overlap(dummy3) ==1:
                                    
                                    pVertex_end_current = pVertex_start    
                                    aCell_flowline.append(lID)
                                    iFlag_found = 1                                                                       
                                    iFlag_previous_overlap =1
                                    
                                    pass
                                else:
                                    if dLength < 0.5*pCell.dLength : #because it taks a short cut                                        
                                        pVertex_end_current = pVertex_start   
                                        iFlag_found = 1  

                                        if iMesh_type ==1 or iMesh_type ==4 or iMesh_type==5:
                                            pass
                                        else:
                                            if iFlag_previous_overlap ==1: 
                                                aCell_flowline.append(lID)              
                                        
                                        iFlag_previous_overlap=0
                                        
                                        pass
                                    else:
                                        pVertex_end_current = pVertex_start    
                                        aCell_flowline.append(lID)
                                        iFlag_found = 1
                                        iFlag_previous_overlap=0
                                        pass

                                    pass
                                

                                pass

                            else:
                                pVertex_end_current = pVertex_start    
                                aCell_flowline.append(lID)                                
                                iFlag_found = 1
                                iFlag_previous_overlap=0
                                pass                            
                                
                            break
                            pass
                        else:
                            iFlag_found2, dummy2 = pCell.which_edge_cross_this_vertex(pVertex_start)
                            iFlag_found3, dummy3 = pCell.which_edge_cross_this_vertex(pVertex_end)
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
            aCell_flowline = aCell_flowline[::-1]
            #simplify list
         
            aCell_simple = simplify_list(aCell_flowline)
            #save the output
            nCell3 = len(aCell_simple)
            aCoordinates = list()
            if nCell3 >1:
                for i in range(nCell3):
                    lIndex = aCell_simple[i]
                    x = aCell_in[lIndex].dX_center
                    y = aCell_in[lIndex].dY_center
                    aCoordinates.append([x,y])
                    pass

                pFlowline = convert_coordinates_to_flowline(aCoordinates)
                pFlowline.iSegment = iSegment_in
                pFlowline.iStream_order=iStream_order_in
                aFlowline_out.append(pFlowline)

            print(iSegment_in, ': ',aCell_simple)

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
    dDiatance_min=0.0
    iFlag_first=1
    for j in range(nCell):
        pCell = aCell_intersect_in[j]
        aFlowline= pCell.aFlowline
        nFlowline = len(aFlowline)
        for i in range(nFlowline):
            pFlowline = aFlowline[i]
            pVertex_start = pFlowline.pVertex_start
            pVertex_end = pFlowline.pVertex_end
            dDiatance = pVertex_end.calculate_distance( pVertex_outlet_in)
            if iFlag_first ==1:
                dDiatance_min = dDiatance               
                lID_outlet = pCell.lIndex               
                lID_outlet2 = j
                lID_outlet3 = i
                iFlag_first=0
            else:
                if  dDiatance < dDiatance_min:
                    dDiatance_min = dDiatance                
                    lID_outlet = pCell.lIndex                
                    lID_outlet2 = j
                    lID_outlet3 = i
                    pass    
                else:
                    #print(dDiatance)
                    pass

            pass     

    #loop through       
    pFlowline = aCell_intersect_in[lID_outlet2].aFlowline[lID_outlet3]
    iSegment =pFlowline.iSegment
    iStream_order=pFlowline.iStream_order
    pVertex_outlet=pFlowline.pVertex_end
        
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
    



    return aCell_out, aFlowline_out, aFlowline_out_no_parallel
