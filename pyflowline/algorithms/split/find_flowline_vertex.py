

import importlib
iFlag_cython = importlib.util.find_spec("cython") 
if iFlag_cython is not None:
    from pyflowline.algorithms.cython.kernel import add_unique_vertex
else:
    from pyflowline.algorithms.auxiliary.find_index_in_list import add_unique_vertex

def find_flowline_vertex(aFlowline_in, dThreshold_in=1.0E-6): 
    """use the set to add unique vertex into a list

    Args:
        aFlowline_in (_type_): _description_
        dThreshold_in (_type_, optional): _description_. Defaults to 1.0E-6.

    Returns:
        _type_: _description_
    """
    aVertex = set()
    nFlowline = len(aFlowline_in) 
    #set up index and id first
    lVertexID = 1
    
    for i in range(0, nFlowline):
        pFlowline = aFlowline_in[i]        
        #update the start and end vertex
        pFlowline.pVertex_start.lVertexID = lVertexID
        lVertexID = lVertexID + 1
        pFlowline.pVertex_end.lVertexID = lVertexID
        lVertexID = lVertexID + 1
        aVertex.add(pFlowline.pVertex_start) 
        aVertex.add(pFlowline.pVertex_end)         
        pass

    #conver the set back to list
    aVertex=list(aVertex)
    #for i in range(0, nFlowline):      
    #    pFlowline = aFlowline_in[i]
    #    pVertex_start = pFlowline.pVertex_start
    #    pVertex_end = pFlowline.pVertex_end        
    #    aVertex, dummy = add_unique_vertex(aVertex, pVertex_start,dThreshold_in)
    #    aVertex, dummy = add_unique_vertex(aVertex, pVertex_end, dThreshold_in)
    #    pass

    return aVertex

def find_flowline_vertex_old(aFlowline_in, dThreshold_in=1.0E-6): 
    aVertex = list()
    nFlowline = len(aFlowline_in)    
    aVertex=list(aVertex)
    for i in range(0, nFlowline):      
        pFlowline = aFlowline_in[i]
        pVertex_start = pFlowline.pVertex_start
        pVertex_end = pFlowline.pVertex_end        
        aVertex, dummy = add_unique_vertex(aVertex, pVertex_start,dThreshold_in)
        aVertex, dummy = add_unique_vertex(aVertex, pVertex_end, dThreshold_in)
        pass

    return aVertex
