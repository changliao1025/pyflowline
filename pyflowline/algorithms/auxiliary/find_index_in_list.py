import copy
import numpy as np

def find_vertex_in_list(aVertex_in, pVertex_in):
    """[find the index of a vertex in a list]

    Args:
        aVertex_in ([type]): [description]
        pVertex_in ([type]): [description]

    Returns:
        [type]: [description]
    """

    iFlag_exist = 0
    lIndex= -1
    nVertex= len(aVertex_in)

    if nVertex > 0 :
        for i in np.arange( nVertex):
            pVertex = aVertex_in[i]
            if pVertex == pVertex_in:
                iFlag_exist = 1      
                lIndex = i 
                break                
            else:
                pass

        pass        
        
    else:
        pass
    
    return iFlag_exist, lIndex



def find_vertex_on_edge(aVertex_in, pEdge_in):
    iFlag_exist = 0
    aIndex= list()
    aIndex_order=list()
    aDistance=list()
    nVertex= len(aVertex_in)
    npoint = 0    
    if nVertex > 0 :
        for i in np.arange( nVertex):
            pVertex = aVertex_in[i]
            iFlag_overlap, dDistance, diff = pEdge_in.check_vertex_on_edge(pVertex)
            if iFlag_overlap == 1:                
                iFlag_exist = 1                      
                aDistance.append(dDistance)
                aIndex.append(i) 
                npoint = npoint + 1          
            else:                
                if  diff < 1.0:
                    iFlag_overlap = pEdge_in.check_vertex_on_edge(pVertex)
                   
                pass

        #re-order 
        if iFlag_exist == 1 :
            x = np.array(aDistance)
            b = np.argsort(x)
            c = np.array(aIndex)
            d= c[b]
            aIndex_order = list(d)        
    else:
        pass
    
    return iFlag_exist, npoint , aIndex_order


def find_edge_in_list(aEdge_in, pEdge_in):
    """[find the index of an edge in a list]

    Args:
        aEdge_in ([type]): [description]
        pEdge_in ([type]): [description]

    Returns:
        [type]: [description]
    """
    iFlag_exist = 0
    lIndex= -1
    nEdge= len(aEdge_in)

    if nEdge > 0 :
        for i in np.arange( nEdge):
            pEdge = aEdge_in[i]
            if pEdge == pEdge_in:
                iFlag_exist = 1      
                lIndex = i 
                break                
            else:
                pass

        pass        
        
    else:
        pass
    
    return iFlag_exist, lIndex


def find_flowline_in_list(aFlowline_in, pFlowline_in):
    """[find the index of a flowline in a list]

    Args:
        aFlowline_in ([type]): [description]
        pFlowline_in ([type]): [description]

    Returns:
        [type]: [description]
    """
    iFlag_exist = 0
    lIndex= -1
    nFlowline= len(aFlowline_in)

    if nFlowline > 0 :
        for i in np.arange( nFlowline):
            pFlowline = aFlowline_in[i]
            if pFlowline == pFlowline_in:
                iFlag_exist = 1      
                lIndex = i 
                break                
            else:
                pass

        pass        
        
    else:
        pass
    
    return iFlag_exist, lIndex

def find_hexagon_through_edge(aHexagon_in, pEdge_in):
    """find the hexagons which contain an edge

    Args:
        aHexagon_in ([type]): [description]
        pEdge_in ([type]): [description]

    Returns:
        [type]: [description]
    """

    nHexagon = len(aHexagon_in)
    aHexagon_out = list()
    for i in range(nHexagon):
        pHexagon = aHexagon_in[i]
        if pHexagon.has_this_edge(pEdge_in) ==1:
            aHexagon_out.append(pHexagon)
            pass
        else:
            pass

    return aHexagon_out

def check_if_duplicates(aList_in):
    """[Check if given list contains any duplicates]

    Returns:
        [type]: [description]
    """
    iFlag_unique = 1
    for elem in aList_in:
        if aList_in.count(elem) > 1:
            iFlag_unique = 0
            break
        else:
            pass
    
    return iFlag_unique


def add_unique_vertex(aVertex_in, pVertex_in):
    """[add a vertex to a list if it is not already included]

    Args:
        aVertex_in ([type]): [description]
        pVertex_in ([type]): [description]

    Returns:
        [type]: [description]
    """
    iFlag_exist = 0
    nVertex = len(aVertex_in)     

    iFlag_exist, dummy =  find_vertex_in_list(aVertex_in, pVertex_in)

    if iFlag_exist == 1:
        pass
    else:
        aVertex_in.append(pVertex_in)
        pass

    return aVertex_in, iFlag_exist


def find_list_in_list(aList_in, pList_in):
    c = copy.deepcopy(pList_in)
    c.sort()
    nList = len(aList_in)
    iFlag = 0
    for i in range(nList):
        a = aList_in[i]
        b = copy.deepcopy(a)
        b.sort()
        if (b == c ):
            iFlag =1 
            break
        else:
            iFlag = 0
   
    return iFlag
