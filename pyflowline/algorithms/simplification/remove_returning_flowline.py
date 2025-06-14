import copy
import numpy as np
from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_flowline

def remove_returning_flowline(iMesh_type_in, aCell_intersect_in, pVertex_outlet_in):
    aFlowline_out=list()
    nCell =  len(aCell_intersect_in)
    #lCellID_to_cell = {cell.lCellID: cell for cell in aCell_intersect_in}
    def simplify_list(aCell_flowline_in):
        aCell_flowline_out = copy.deepcopy(aCell_flowline_in)
        nCell2 = len(aCell_flowline_out)
        #iFlag_unique = check_if_duplicates(aCell_flowline_out)
        if int(len(aCell_flowline_out) == len(set(aCell_flowline_out))) :
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

            aCell_flowline_out = simplify_list(aCell_flowline_out)
        return   aCell_flowline_out

    def retrieve_flowline_intersect_index(iSegment_in, iStream_order_in, pVertex_end_in):
        lCellID =-1
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
                lCellID = pCell.lCellID
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
                                    iFlag_found = 1
                                    iFlag_previous_overlap =1

                                    pass
                                else:
                                    if dLength < (0.1*pCell.dLength) : #because it taks a short cut
                                        pVertex_end_current = pVertex_start
                                        iFlag_found = 1
                                        if iFlag_previous_overlap ==1:
                                            aCell_flowline.append(lCellID)
                                            iFlag_previous_overlap=0
                                            pass
                                        else:
                                            iFlag_previous_overlap=1
                                            pass
                                        pass
                                    else:
                                        pVertex_end_current = pVertex_start
                                        aCell_flowline.append(lCellID)
                                        iFlag_found = 1


                            else:
                                pVertex_end_current = pVertex_start
                                aCell_flowline.append(lCellID)
                                iFlag_found = 1
                                iFlag_previous_overlap=0

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

    dDistance_min = float('inf')
    lCellID_outlet = lCellIndex_outlet = lFlowlineIndex_outlet = None
    for j, pCell in enumerate(aCell_intersect_in):
        for i, pFlowline in enumerate(pCell.aFlowline):
            dDiatance = pFlowline.pVertex_end.calculate_distance(pVertex_outlet_in)
            if dDiatance < dDistance_min:
                dDistance_min = dDiatance
                lCellID_outlet = pCell.lCellID
                lCellIndex_outlet = j
                lFlowlineIndex_outlet = i

    #loop through
    pFlowline = aCell_intersect_in[lCellIndex_outlet].aFlowline[lFlowlineIndex_outlet]
    iStream_segment =pFlowline.iStream_segment
    iStream_order=pFlowline.iStream_order
    pVertex_outlet=pFlowline.pVertex_end

    retrieve_flowline_intersect_index(iStream_segment, iStream_order, pVertex_outlet)

    #update outlet location since the conceptual flowline use center
    pVertex_outlet =  aFlowline_out[0].pVertex_end

    return aFlowline_out, lCellID_outlet, pVertex_outlet
