import importlib.util
from typing import List, Optional
iFlag_cython = importlib.util.find_spec("cython")
if iFlag_cython is not None:
    from pyflowline.algorithms.cython.kernel import add_unique_vertex
else:
    from pyflowline.algorithms.auxiliary.find_index_in_list import add_unique_vertex

# Import for the new graph-based approach
from pyflowline.classes.rivergraph import pyrivergraph
from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.flowline import pyflowline

def find_flowline_vertex(aFlowline_in, dThreshold_in=1.0E-6):
    """use the set to add unique vertex into a list

    Args:
        aFlowline_in (_type_): _description_
        dThreshold_in (_type_, optional): _description_. Defaults to 1.0E-6.

    Returns:
        _type_: _description_
    """
    aVertex = set()
    lVertexID = 1

    for pFlowline in aFlowline_in:
        pFlowline.pVertex_start.lVertexID = lVertexID
        pFlowline.pVertex_start.lFlowlineID = pFlowline.lFlowlineID
        lVertexID += 1
        pFlowline.pVertex_end.lVertexID = lVertexID
        pFlowline.pVertex_end.lFlowlineID = pFlowline.lFlowlineID
        lVertexID += 1
        aVertex.update([pFlowline.pVertex_start, pFlowline.pVertex_end])

    #conver the set back to list
    aVertex=list(aVertex)

    return aVertex

def find_flowline_vertex_graph(aFlowline_in: List[pyflowline],
                              dThreshold_in: float = 1.0E-6) -> List[pyvertex]:
    """
    Graph-based replacement for find_flowline_vertex using pyrivergraph.

    This function provides the same functionality as find_flowline_vertex but uses
    the pyrivergraph class for vertex extraction. It offers better performance for
    large networks and enables additional graph-based analysis capabilities.

    Args:
        aFlowline_in (List[pyflowline]): List of flowline objects
        dThreshold_in (float, optional): Distance threshold (maintained for compatibility,
                                       not used in graph-based approach). Defaults to 1.0E-6.

    Returns:
        List[pyvertex]: List of unique vertices with assigned vertex IDs

    Example:
        >>> flowlines = [flowline1, flowline2, flowline3]
        >>> vertices = find_flowline_vertex_graph(flowlines)
        >>> print(f"Found {len(vertices)} unique vertices")
    """
    if not aFlowline_in:
        return []

    # Create river graph and extract vertices
    river_graph = pyrivergraph(aFlowline_in)
    return river_graph.get_vertices()

def find_flowline_vertex_enhanced(aFlowline_in: List[pyflowline],
                                 dThreshold_in: float = 1.0E-6,
                                 use_graph: bool = True) -> List[pyvertex]:
    """
    Enhanced vertex extraction with option to use graph-based or set-based approach.

    This function allows choosing between the original set-based approach and the
    new graph-based approach. The graph-based approach is recommended for better
    performance and additional capabilities.

    Args:
        aFlowline_in (List[pyflowline]): List of flowline objects
        dThreshold_in (float, optional): Distance threshold for set-based approach.
                                       Defaults to 1.0E-6.
        use_graph (bool, optional): If True, use graph-based approach; if False,
                                  use original set-based approach. Defaults to True.

    Returns:
        List[pyvertex]: List of unique vertices with assigned vertex IDs

    Example:
        >>> # Use graph-based approach (recommended)
        >>> vertices = find_flowline_vertex_enhanced(flowlines, use_graph=True)
        >>>
        >>> # Use original set-based approach
        >>> vertices = find_flowline_vertex_enhanced(flowlines, use_graph=False)
    """
    if use_graph:
        return find_flowline_vertex_graph(aFlowline_in, dThreshold_in)
    else:
        return find_flowline_vertex(aFlowline_in, dThreshold_in)

# Compatibility alias - can be used to gradually migrate existing code
find_flowline_vertex_new = find_flowline_vertex_graph

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
