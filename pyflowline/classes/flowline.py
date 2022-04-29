from abc import ABCMeta, abstractmethod

import copy

import json
from json import JSONEncoder
import numpy as np
from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.edge import pyedge


class FlowlineClassEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.float32):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, list):
            pass  
        if isinstance(obj, pyvertex):
            return json.loads(obj.tojson()) #lVertexID
        if isinstance(obj, pyedge):
            return obj.lEdgeID         
            
        return JSONEncoder.default(self, obj)

class pyflowline(object):
    __metaclass__ = ABCMeta 

    lFlowlineID=-1
    lIndex=-1
    lIndex_upstream=-1
    lIndex_downstream=-1

    iFlag_dam = 0
    lNHDPlusID=-1

    pVertex_start=None
    pVertex_end=None
    aEdge=None
    aVertex=None

    dLength=0.0
    dSinuosity=0.0

    iStream_segment=-1
    iStream_order =-1

    nEdge=0
    nVertex=0
    iFlag_right = 0
    iFlag_left = 0
    aFlowlineID_start_start = None
    aFlowlineID_start_end = None
    aFlowlineID_end_start = None
    aFlowlineID_end_end = None
    
    def __init__(self, aEdge):    
        self.aEdge = aEdge
        nEdge  = len(aEdge)
        self.nEdge = nEdge
        self.pVertex_start = aEdge[0].pVertex_start
        self.pVertex_end =  aEdge[ nEdge-1  ].pVertex_end
        nVertex = nEdge +1
        self.aVertex=list()
        for i in range(nEdge):
            self.aVertex.append( aEdge[i].pVertex_start )
            pass

        self.aVertex.append( aEdge[nEdge-1].pVertex_end )
        self.nVertex = nVertex
        self.dLength= self.calculate_length()
        self.aFlowlineID_start_start = list()
        self.aFlowlineID_start_end = list()
        self.aFlowlineID_end_start = list()
        self.aFlowlineID_end_end = list()
     
        return

    def calculate_length(self):
        dLength =0.0
        #loop though
        for edge in self.aEdge:
            edge.calculate_length()
            dLength = dLength + edge.dLength

        #assing
        self.dLength= dLength

        return dLength

    
    
    def check_upstream(self, other):
        iFlag_upstream =-1
        v0 = self.pVertex_start
        v1 = self.pVertex_end

        v2 = other.pVertex_start
        v3 = other.pVertex_end

        if v0 == v3:
            iFlag_upstream =1
        else:
            iFlag_upstream=0

        return iFlag_upstream

    def check_downstream(self, other):
        iFlag_downstream =-1
        v0 = self.pVertex_start
        v1 = self.pVertex_end
        v2 = other.pVertex_start
        v3 = other.pVertex_end
        if v1 == v2:
            iFlag_downstream =1
        else:
            iFlag_downstream=0

        return iFlag_downstream
    
    def reverse(self):
        '''
        reverse the direction of a flowline
        '''
        aVertex = self.aVertex 
        nVertex = self.nVertex
        aVertex_new = list()
        for i in range(nVertex-1,-1,-1) :
            aVertex_new.append( aVertex[i] )

        self.aVertex = aVertex_new
        nVertex  = len(aVertex)
        aEdge = list()
        for i in range(nVertex-1):
            pEdge = pyedge( self.aVertex[i], self.aVertex[i+1] )
            aEdge.append(pEdge)
            pass
        
        self.aEdge = aEdge
        nEdge = len(aEdge)
        self.pVertex_start = aEdge[0].pVertex_start
        self.pVertex_end =  aEdge[ nEdge-1  ].pVertex_end

    def merge_upstream(self, other):
        pFlowline_out = copy.deepcopy(other)    
        
        iFlag_dam1 = other.iFlag_dam
        pVertex_start1 = other.pVertex_start
        pVertex_end1 = other.pVertex_end
        nVertex1 = other.nVertex
        nEdge1 = other.nEdge

        pVertex_start2 = self.pVertex_start
        pVertex_end2 = self.pVertex_end
        nVertex2 = self.nVertex
        nEdge2 = self.nEdge
        iFlag_dam2 = self.iFlag_dam


        if pVertex_end1 == pVertex_start2:
            #this is the supposed operation because they should connect

            nVertex = nVertex1 + nVertex2 - 1
            nEdge = nVertex -1 
            aEdge = copy.deepcopy(other.aEdge )
            for i in range(nEdge2):
                aEdge.append( self.aEdge[i] )
                pass

            aVertex = copy.deepcopy(other.aVertex)
            for i in range(1, nVertex2):
                aVertex.append( self.aVertex[i] )
                pass
            
            pFlowline_out.iFlag_dam = max(iFlag_dam1, iFlag_dam2)
            pFlowline_out.aEdge = aEdge
            pFlowline_out.aVertex = aVertex
            pFlowline_out.nEdge = nEdge
            pFlowline_out.nVertex = nVertex
            pFlowline_out.dLength = self.dLength + other.dLength
            pFlowline_out.pVertex_start = pVertex_start1
            pFlowline_out.pVertex_end = pVertex_end2
            pass
        else:
            pass

        return pFlowline_out

    def calculate_flowline_sinuosity(self):
        pVertex_start = self.pVertex_start
        pVertex_end = self.pVertex_end
        dDistance = pVertex_start.calculate_distance(pVertex_end)
        self.dSinuosity = self.dLength / dDistance
        return

    def __eq__(self, other):                       
        iFlag_overlap = 0 
        nEdge1 = self.nEdge
        nEdge2 = other.nEdge
        if nEdge1 == nEdge2:
            for i in np.arange( nEdge1):
                pEdge1 = self.aEdge[i]
                pEdge2 = other.aEdge[i]
                if pEdge1 == pEdge2:
                    iFlag_overlap =1 
                else:
                    iFlag_overlap =0 
                    break                
            
        else:
            iFlag_overlap = 0

        return iFlag_overlap

    def __ne__(self, other):
        return not self.__eq__(other)
    
    def tojson(self):
        aSkip = ['aEdge', \
                'aVertex','aFlowlineID_start_start','aFlowlineID_start_end',
                'aFlowlineID_end_start','aFlowlineID_end_end']

        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)
        sJson = json.dumps(obj,  \
            sort_keys=True, \
                indent = 4, \
                    ensure_ascii=True, \
                        cls=FlowlineClassEncoder)
        return sJson
        