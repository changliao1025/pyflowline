import copy
import numpy as np
from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_flowline
from pyflowline.algorithms.auxiliary.find_index_in_list import check_if_duplicates

def remove_returning_flowline(iMesh_type_in, aCell_intersect_in, pVertex_outlet_in):    
    aFlowline_out=list()    
    nCell =  len(aCell_intersect_in)
    def simplify_list(aCell_flowline_in):
        aCell_flowline_out = copy.deepcopy(aCell_flowline_in)
        nCell2 = len(aCell_flowline_out)          
        iFlag_unique = check_if_duplicates(aCell_flowline_out)
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
                else:
                    pass

                pass             
            
            aCell_flowline_out = simplify_list(aCell_flowline_out)
        return   aCell_flowline_out
    
    def simplify_list2(aCell_flowline_in):
        aCell_flowline_out = copy.deepcopy(aCell_flowline_in)
        nCell = len(aCell_flowline_out)

        # Check if list contains duplicates
        iFlag_unique = check_if_duplicates(aCell_flowline_out)
        if iFlag_unique == 1:
            return aCell_flowline_out

        # Create a dictionary to store counts of elements
        elem_counts = {}
        for elem in aCell_flowline_out:
            elem_counts[elem] = elem_counts.get(elem, 0) + 1

        i = 0
        while i < nCell - 1:
            elem = aCell_flowline_out[i]
            if elem_counts[elem] > 1:
                indices = [j for j, x in enumerate(aCell_flowline_out) if x == elem]
                start = indices[0]
                end = indices[-1]
                del aCell_flowline_out[start + 1 : end]
                nCell = len(aCell_flowline_out)
            else:
                i += 1

        return aCell_flowline_out

    def retrieve_flowline_intersect_index(iSegment_in, iStream_order_in, pVertex_end_in):
        lCellID =-1
        if iSegment_in == 608:
            print('here')
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
            #if lCellID == 110191:
            #    print('here')
            #    pass         
            for j in range(nCell):
                pCell = aCell_intersect_in[j]
                lCellID = pCell.lCellID
                #if lCellID == 110921:
                #    print('here')
                aFlowline= pCell.aFlowline
                nFlowline = len(aFlowline)
                for i in range(nFlowline):
                    pFlowline = aFlowline[i]
                    pVertex_start = pFlowline.pVertex_start
                    pVertex_end = pFlowline.pVertex_end
                    iStream_segment = pFlowline.iStream_segment
                    iStream_order = pFlowline.iStream_order
                    if pVertex_end == pVertex_end_current: #be careful because this function is strict
                        if iStream_segment == iSegment_in:                            
                            dLength = pFlowline.dLength
                            iFlag_found2, dummy2 = pCell.which_edge_cross_this_vertex(pVertex_start)
                            iFlag_found3, dummy3 = pCell.which_edge_cross_this_vertex(pVertex_end)

                            if iFlag_found2 == 1 and iFlag_found3 == 1: #on the edge
                                #same edge
                                if dummy2.is_overlap(dummy3) ==1:                                    
                                    pVertex_end_current = pVertex_start    
                                    aCell_flowline.append(lCellID)
                                    #print(lCellID)
                                    iFlag_found = 1                
                                    iFlag_previous_overlap =1
                                    
                                    pass
                                else:
                                    if dLength < (0.1*pCell.dLength) : #because it taks a short cut                                        
                                        pVertex_end_current = pVertex_start   
                                        iFlag_found = 1                                          
                                        if iFlag_previous_overlap ==1: 
                                            aCell_flowline.append(lCellID)              
                                            #print(lCellID)
                                            iFlag_previous_overlap=0
                                            pass
                                        else:
                                            iFlag_previous_overlap=1
                                            pass                                             
                                        pass
                                    else:
                                        pVertex_end_current = pVertex_start    
                                        aCell_flowline.append(lCellID)
                                        #print(lCellID)
                                        iFlag_found = 1                                      
                                        pass

                                    pass                                

                                pass

                            else:
                                pVertex_end_current = pVertex_start    
                                aCell_flowline.append(lCellID)                                
                                iFlag_found = 1
                                iFlag_previous_overlap=0
                                pass                            
                                
                            break
                            
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
                                aSegment_upstream.append(iStream_segment)
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
            #if 112378 in aCell_flowline:
            #    print(aCell_flowline)
            #    print('here')
            #    pass
            #simplify list         
            aCell_simple = simplify_list(aCell_flowline)
            #save the output
            nCell3 = len(aCell_simple)
            aCoordinates = list()            
            if nCell3 >1:
                for i in range(nCell3):
                    lCellID = aCell_simple[i]
                    for j in range( nCell):
                        if aCell_intersect_in[j].lCellID == lCellID:    
                            x = aCell_intersect_in[j].dLongitude_center_degree
                            y = aCell_intersect_in[j].dLatitude_center_degree                
                            aCoordinates.append([x,y])
                    pass
                pFlowline = convert_gcs_coordinates_to_flowline(aCoordinates)
                pFlowline.iStream_segment = iSegment_in
                pFlowline.iStream_order=iStream_order_in
                aFlowline_out.append(pFlowline)

            sort_index = np.argsort(aStream_order)            
            #sort_index = sort_index[::-1] #can we process low order first?
            #nUpstream = len(aSegment_upstream)
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
                lCellID_outlet = pCell.lCellID               
                lCellID_outlet2 = j
                lCellID_outlet3 = i
                iFlag_first=0
            else:
                if  dDiatance < dDiatance_min:
                    dDiatance_min = dDiatance                
                    lCellID_outlet = pCell.lCellID                
                    lCellID_outlet2 = j
                    lCellID_outlet3 = i
                    pass    
                else:
                    pass

            pass     

    #loop through       
    pFlowline = aCell_intersect_in[lCellID_outlet2].aFlowline[lCellID_outlet3]
    iStream_segment =pFlowline.iStream_segment
    iStream_order=pFlowline.iStream_order
    pVertex_outlet=pFlowline.pVertex_end
        
    retrieve_flowline_intersect_index(iStream_segment, iStream_order, pVertex_outlet)

    return aFlowline_out, lCellID_outlet, pVertex_outlet
