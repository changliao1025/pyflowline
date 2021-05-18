from abc import ABCMeta, abstractmethod
import json

from osgeo import gdal, osr, ogr
from pystream.shared.vertex import pyvertex
from pystream.shared.edge import pyedge

class pyflowline(object):
    __metaclass__ = ABCMeta 
    aEdge=None
    aVertex=None
    dLength=0.0

    lIndex=-1
    lIndex_upstream=-1
    lIndex_downstream=-1

    #def __init__(self, aEdge):    
    #    self.aEdge = aEdge
    #    nEdge  = len(aEdge)
    #    self.nEdge = nEdge
    #    self.pVertex_start = aEdge[0].pVertex_start
    #    self.pVertex_end =  aEdge[ nEdge-1  ].pVertex_end
    #    return

    def __init__(self, aCoordinate):
        

        self.aVertex = aCoordinate
        nVertex  = len(aCoordinate)
        self.nVertex = nVertex
        self.nEdge = nVertex-1
        self.aEdge=list()
        for i in range(nVertex-1):
            pEdge = pyedge( aCoordinate[i], aCoordinate[i+1] )
            self.aEdge.append(pEdge)
            pass

        self.pVertex_start = self.aEdge[0].pVertex_start
        self.pVertex_end =  self.aEdge[nVertex-2  ].pVertex_end

        return

    def calculate_length(self):
        dLength =0.0
        #loop though
        for edge in self.aEdge:
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
