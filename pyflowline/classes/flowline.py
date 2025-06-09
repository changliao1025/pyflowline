
import os
import copy
import json
from json import JSONEncoder
import numpy as np
from osgeo import ogr, gdal, osr
from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.edge import pyedge
from pyflowline.formats.export_vertex import export_vertex_as_polygon

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
            return json.loads(obj.tojson())
        if isinstance(obj, pyedge):
            return obj.lEdgeID

        return JSONEncoder.default(self, obj)

class pyflowline(object):
    """The pyflowline class

    Args:
        object (object): None

    Returns:
        pyflowline: The flowline object
    """

    lFlowlineID=-1
    lFlowlineID_downstream=-1 #if braided, then we need a list
    lFlowlineIndex=-1
    lIndex_upstream=-1
    lIndex_downstream=-1

    iFlag_keep = 1 #used for simplification algorithm

    iFlag_dam = 0
    lNHDPlusID=-1

    pVertex_start=None
    pVertex_end=None
    aEdge=None
    aVertex=None

    dLength=0.0
    dSinuosity=0.0
    dDrainage_area=0.0

    iStream_segment=-1
    iStream_order = -1

    nEdge=0
    nVertex=0
    iFlag_right = 0
    iFlag_left = 0
    aFlowlineID_start_start = None
    aFlowlineID_start_end = None
    aFlowlineID_end_start = None
    aFlowlineID_end_end = None

    #for stream topology
    lFlowlineIndex_downstream = None #only store the index, not the actual objects
    aFlowline_upstream = None

    pBound=None

    def __init__(self, aEdge):
        """
        Initilize a flowline object

        Args:
            aEdge (list [pyedge]): A list of edge objects
        """
        self.aEdge = aEdge
        nEdge = len(aEdge)
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
        self.iFlag_keep = 1

        #also build wkt string
        pGeometry = ogr.Geometry(ogr.wkbLineString)
        for i in range(nVertex):
            pGeometry.AddPoint(self.aVertex[i].dLongitude_degree, self.aVertex[i].dLatitude_degree)

        self.wkt = pGeometry.ExportToWkt()
        pGeometry = None

        self.calculate_flowline_bound()

        return

    def __hash__(self):
        return hash((self.pVertex_start, self.pVertex_end))

    def calculate_length(self):
        """
        Calcualte the length

        Returns:
            float: The length of the flowline
        """
        #dLength =0.0
        #loop though
        #for i in range(self.nEdge):
        #for edge in self.aEdge:
        #    self.aEdge[i].calculate_length()
        #    dLength = dLength + self.aEdge[i].dLength

        #assing
        #self.dLength= dLength


        self.dLength = sum(edge.dLength for edge in self.aEdge)
        return self.dLength

    def calculate_flowline_bound(self):
        dLat_min = 90
        dLat_max = -90
        dLon_min = 180
        dLon_max = -180
        for i in range(self.nVertex):
            dLon_max = np.max( [dLon_max, self.aVertex[i].dLongitude_degree] )
            dLon_min = np.min( [dLon_min, self.aVertex[i].dLongitude_degree] )
            dLat_max = np.max( [dLat_max, self.aVertex[i].dLatitude_degree] )
            dLat_min = np.min( [dLat_min, self.aVertex[i].dLatitude_degree] )

        self.pBound = (float(dLon_min), float(dLat_min), float(dLon_max), float(dLat_max))
        return self.pBound

    def check_upstream(self, other):
        """
        Check whether another flowline is upstream or not

        Args:
            other (pyflowline): The other flowline

        Returns:
            int: 1 if it is, 0 if not
        """
        iFlag_upstream = -1
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
        """
        Check whether another flowline is downstream or not

        Args:
            other (pyflowline): The other flowline

        Returns:
            int: 1 if it is, 0 if not
        """
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
        """
        Reverse a flowline
        """
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
        """
        Merge two flowlines as one

        Args:
            other (pyflowline): The other flowline

        Returns:
            pyflowline: The merged flowline
        """
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

        #is this necessary
        #pFlowline_out.iStream_segment = self.iStream_segment

        return pFlowline_out

    def split_edge_by_length(self, dDistance):
        """
        Split a flowline using the length threshold

        Args:
            dDistance (float): The length threshold for each edge

        Returns:
            pyflowline: The updated flowline
        """
        aEdge=list()
        pFlowline_out=None
        for edge in self.aEdge:
            #edge.calculate_length()
            if edge.dLength > dDistance:
                #break it
                aEdge0=edge.split_by_length(dDistance)
                for edge0 in aEdge0:
                    aEdge.append(edge0)
                pass
            else:
                aEdge.append(edge)
                pass
        pFlowline_out=pyflowline(aEdge)
        #copy the attributes
        pFlowline_out.copy_attributes(self)

        return pFlowline_out

    def split_by_length(self, dDistance):
        aFlowline = list()
        if self.dLength<=dDistance:
            aFlowline.append(self)
        else:
            #split using the half of the edge
            nEdge = len(self.aEdge)
            first_leg = [0, int(nEdge/2)]
            second_leg = [int(nEdge/2), nEdge]
            aEdge0 = self.aEdge[first_leg[0]:first_leg[1]]
            aEdge1 = self.aEdge[second_leg[0]:second_leg[1]]
            pFlowline0 = pyflowline(aEdge0)
            pFlowline1 = pyflowline(aEdge1)
            pFlowline0.copy_attributes(self)
            pFlowline1.copy_attributes(self)
            if pFlowline0.dLength > dDistance:
                aFlowline.extend(pFlowline0.split_by_length(dDistance))
            else:
                aFlowline.append(pFlowline0)
            if pFlowline1.dLength > dDistance:
                aFlowline.extend(pFlowline1.split_by_length(dDistance))
            else:
                aFlowline.append(pFlowline1)
            pass

        return aFlowline

    def copy_attributes(self, other):
        """
        Copy the attributes from another flowline

        Args:
            other (pyflowline): The other flowline
        """
        self.lFlowlineID = other.lFlowlineID
        self.lFlowlineID_downstream = other.lFlowlineID_downstream
        self.lFlowlineIndex = other.lFlowlineIndex
        self.iFlag_dam = other.iFlag_dam
        self.lNHDPlusID = other.lNHDPlusID
        self.iStream_segment = other.iStream_segment
        self.iStream_order = other.iStream_order
        self.dDrainage_area = other.dDrainage_area
        self.iFlag_right = other.iFlag_right
        self.iFlag_left = other.iFlag_left
        self.iFlag_keep = other.iFlag_keep
        self.lFlowlineIndex_downstream = other.lFlowlineIndex_downstream
        return

    def calculate_flowline_sinuosity(self):
        """
        Calculate the sinuosoty of a flowline
        """
        pVertex_start = self.pVertex_start
        pVertex_end = self.pVertex_end
        dDistance = pVertex_start.calculate_distance(pVertex_end)
        self.dSinuosity = self.dLength / dDistance
        return

    def calculate_distance_to_vertex(self, pVertex):
        dDistance_min_vertex = 1.0E10
        pVertex_min_vertex = None
        for i in range(self.nVertex):
            dDistance = pVertex.calculate_distance(self.aVertex[i])
            if dDistance < dDistance_min_vertex:
                dDistance_min_vertex = dDistance
                pVertex_min_vertex = self.aVertex[i]
            pass

        dDistance_min_edge = 1.0E10
        pVertex_min_edge = None
        for i in range(self.nEdge):
            dDistance, pVertex_out_edge = self.aEdge[i].calculate_distance_to_vertex(pVertex)
            if dDistance < dDistance_min_edge:
                dDistance_min_edge = dDistance
                pVertex_min_edge = pVertex_out_edge
            pass

        if dDistance_min_vertex < dDistance_min_edge:
            dDistance_min = dDistance_min_vertex
            pVertex_out = pVertex_min_vertex
        else:
            dDistance_min = dDistance_min_edge
            pVertex_out = pVertex_min_edge
            pass

        return dDistance_min, pVertex_out

    def calculate_buffer_zone_polygon(self, dRadius,sFilename_out = None, sFolder_out=None):
        """
        Calculate the buffer zone polygon

        Args:
            dRadius (float): The buffer zone distance

        Returns:
            list: A list of buffer zone points
        """
        pMultiPolygon = ogr.Geometry(ogr.wkbMultiPolygon)
        aVertex_out = list()
        aCircle_out = list()
        for i in range(self.nEdge):
            edge = self.aEdge[i]
            aVertex, aVertex_center, aVertex_circle, aCircle = edge.calculate_buffer_zone_polygon(dRadius)

            aCircle_out.append(aCircle)
            #create a polygon feature for each edge
            ring = ogr.Geometry(ogr.wkbLinearRing)
            for pVertex in aVertex:
                ring.AddPoint(pVertex.dLongitude_degree, pVertex.dLatitude_degree)

            ring.AddPoint(aVertex[0].dLongitude_degree, aVertex[0].dLatitude_degree) #close the ring
            # Create a polygon and add the ring to it
            pPolygon = ogr.Geometry(ogr.wkbPolygon)
            pPolygon.AddGeometry(ring)
            # Add the polygon to the MultiPolygon
            pMultiPolygon.AddGeometry(pPolygon)
            #save out for debug
            if sFolder_out is not None:
                sFilename_dummy= os.path.join(sFolder_out, 'buffer_zone_edge_%d.geojson' % i)
                export_vertex_as_polygon(aVertex, sFilename_dummy)

        pUnionPolygon = pMultiPolygon.UnionCascaded()
        for i in range(pUnionPolygon.GetGeometryRef(0).GetPointCount()):
            lon, lat, _ = pUnionPolygon.GetGeometryRef(0).GetPoint(i)
            point2= dict()
            point2['dLongitude_degree'] = lon
            point2['dLatitude_degree'] =  lat
            pVertex2 = pyvertex(point2)
            aVertex_out.append(pVertex2)

        if sFilename_out is not None:
            export_vertex_as_polygon(aVertex_out, sFilename_out)

        return aVertex_out, aCircle_out

    def __eq__(self, other):
        """
        Check whether two flowline are equivalent

        Args:
            other (pyflowline): The other flowline

        Returns:
            int: 1 if equivalent, 0 if not
        """
        if len(self.aEdge) != len(other.aEdge):
            return 0

        return int(all(edge1 == edge2 for edge1, edge2 in zip(self.aEdge, other.aEdge)))

    def __ne__(self, other):
        """
        Check whether two flowline are equivalent

        Args:
            other (pyflowline): The other flowline

        Returns:
            int: 0 if equivalent, 1 if not
        """
        return not self.__eq__(other)

    def tojson(self):
        """
        Convert a pyflowline object to a json string

        Returns:
            json str: A json string
        """
        aSkip = ['aEdge',
                'aVertex','aFlowlineID_start_start','aFlowlineID_start_end',
                'aFlowlineID_end_start','aFlowlineID_end_end']

        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)
            pass


        sJson = json.dumps(obj,
            sort_keys=True,
                indent = 4,
                    ensure_ascii=True,
                        cls=FlowlineClassEncoder)
        return sJson
