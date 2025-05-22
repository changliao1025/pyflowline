import os
import json
from json import JSONEncoder
import importlib.util
import numpy as np
from pyearth.gis.geometry.calculate_intersect_on_great_circle import calculate_intersect_on_great_circle

from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.circle import pycircle
from pyflowline.algorithms.split.split_by_length import split_edge_by_length
from pyflowline.formats.export_vertex import export_vertex_as_polygon
from pyflowline.algorithms.find_minimal_enclosing_polygon import find_minimal_enclosing_polygon

iFlag_cython = importlib.util.find_spec("cython")
if iFlag_cython is not None:
    from pyflowline.algorithms.cython.kernel import calculate_distance_based_on_longitude_latitude
    from pyflowline.algorithms.cython.kernel import calculate_angle_betwen_vertex
    from pyflowline.algorithms.cython.kernel import calculate_distance_to_plane
else:
    from pyearth.gis.geometry.calculate_distance_based_on_longitude_latitude import calculate_distance_based_on_longitude_latitude
    from pyearth.gis.geometry.calculate_angle_betwen_vertex import calculate_angle_betwen_vertex
    from pyearth.gis.geometry.calculate_distance_to_plane import calculate_distance_to_plane

class EdgeClassEncoder(JSONEncoder):
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


        return JSONEncoder.default(self, obj)

class pyedge(object):
    """The pyedge class

    Args:
        object (object): None

    Returns:
        pyedge: A edge object
    """


    lEdgeID=-1
    lEdgeIndex=-1
    pVertex_start = None
    pVertex_end = None
    dLength=0.0
    lIndex_upstream=-1
    lIndex_downstream=-1
    pBound=None

    def __init__(self, pVertex_start_in, pVertex_end_in):
        """
        Initilize a pyedge object

        Args:
            pVertex_start_in (pyvertex): The starting vertex
            pVertex_end_in (pyvertex): The ending vertex
        """
        if pVertex_start_in == pVertex_end_in:
            print('The two vertices are the same')
            return None
        else:
            self.pVertex_start = pVertex_start_in
            self.pVertex_end = pVertex_end_in
            self.dLength = self.calculate_length()
            self.calculate_edge_bound()

        return

    @classmethod
    def create(cls, pVertex_start_in, pVertex_end_in):
        """
        Factory method to create a pyedge object

        Args:
            pVertex_start_in (pyvertex): The starting vertex
            pVertex_end_in (pyvertex): The ending vertex

        Returns:
            pyedge or None: A pyedge object if valid, otherwise None
        """
        if pVertex_start_in == pVertex_end_in:
            print("The two vertices are the same. Returning None.")
            return None
        return cls(pVertex_start_in, pVertex_end_in)

    def calculate_edge_bound(self):

        dLon_max = np.max( [self.pVertex_start.dLongitude_degree, self.pVertex_end.dLongitude_degree] )
        dLon_min = np.min( [self.pVertex_start.dLongitude_degree, self.pVertex_end.dLongitude_degree] )
        dLat_max = np.max( [self.pVertex_start.dLatitude_degree, self.pVertex_end.dLatitude_degree] )
        dLat_min = np.min( [self.pVertex_start.dLatitude_degree, self.pVertex_end.dLatitude_degree] )

        self.pBound = (dLon_min, dLat_min, dLon_max, dLat_max)
        return self.pBound

    def calculate_length(self):
        """
        Calcualate the length of the edge

        Returns:
            float: The length of the edge
        """
        dLength =0.0
        dLength = self.pVertex_start.calculate_distance( self.pVertex_end)
        self.dLength= dLength
        return dLength

    def check_shared_vertex(self, other):
        """
        Check whether two edges are sharing the same vertex

        Args:
            other (pyedge): The other edge object to be checked

        Returns:
            int: Flag, 1: shared; 0: non-sharing
        """
        iFlag_shared = -1
        v0 = self.pVertex_start
        v1 = self.pVertex_end
        v2 = other.pVertex_start
        v3 = other.pVertex_end
        if v0 == v2 or v0 == v3 or v1==v2 or v1==v3:
            iFlag_shared = 1
        else:
            iFlag_shared = 0

        return iFlag_shared

    def check_upstream(self, other):
        """
        Check whether another edge is the upstream of current edge

        Args:
            other (pyedge): The other edge object to be checked

        Returns:
            int: Flag, 1: upstream; 0: non-upstream
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
        Check whether another edge is the downstream of current edge


        Args:
            other (pyedge): The other edge object to be checked

        Returns:
            int: Flag, 1: downstream; 0: non-downstream
        """
        iFlag_downstream = -1
        v0 = self.pVertex_start
        v1 = self.pVertex_end
        v2 = other.pVertex_start
        v3 = other.pVertex_end
        if v1 == v2:
            iFlag_downstream =1
        else:
            iFlag_downstream=0

        return iFlag_downstream

    def split_by_length(self,dLength_in):
        """
        Split an edge using the threshold

        Args:
            dLength_in (float): The length threshold

        Returns:
            list [pyedge]: A list of edge objects, length of 1 if it meets the requirement
        """
        aEdge_out=list()
        if self.dLength <= dLength_in:
            aEdge_out.append(self)
        else:
            #find location from up to down
            aEdge_out = split_edge_by_length(self, dLength_in)
            pass
        return aEdge_out

    def reverse(self):
        """
        Reverse an edge
        """
        v0 = self.pVertex_start
        v1 = self.pVertex_end
        self.pVertex_start = v1
        self.pVertex_end = v0
        return

    def is_overlap(self, pEdge_in):
        """
        Check if two edges overlap each other

        Args:
            pEdge_in (pyedge): The other edge to be checked

        Returns:
            int: 1 if overlap; 0 if not
        """
        iFlag_overlap = 0
        pVertex_start1 = self.pVertex_start
        pVertex_end1 = self.pVertex_end
        pVertex_start2 = pEdge_in.pVertex_start
        pVertex_end2 = pEdge_in.pVertex_end

        if pVertex_start1 == pVertex_start2 and pVertex_end1 == pVertex_end2:
            iFlag_overlap = 1
        else:
            if  pVertex_start1 == pVertex_end2 and pVertex_end1 == pVertex_start2:
                iFlag_overlap = 1
            else:
                iFlag_overlap = 0

        return iFlag_overlap

    def check_vertex_on_edge(self, pVertex_in):
        """
        Check if a vertex on an edge

        Args:
            pVertex_in (pyvertex): The vertex to be checked

        Returns:
            tuple[int, float, float]: 1 if it is on; 0 if not. Length and distance are calculated if on.
        """
        iFlag =0
        dDistance = -1
        dDistance_plane = 9999
        pVertex_start = self.pVertex_start
        pVertex_end = self.pVertex_end
        self.dLength = pVertex_start.calculate_distance(pVertex_end)
        if pVertex_in != pVertex_start and pVertex_in!=pVertex_end:
            d1 = pVertex_start.calculate_distance(pVertex_in)
            d2 = pVertex_end.calculate_distance(pVertex_in)
            d3 = d1 + d2 - self.dLength
            angle3deg = calculate_angle_betwen_vertex(\
                 pVertex_start.dLongitude_degree, pVertex_start.dLatitude_degree,\
                 pVertex_in.dLongitude_degree, pVertex_in.dLatitude_degree,\
                 pVertex_end.dLongitude_degree,pVertex_end.dLatitude_degree)

            dDistance_plane = calculate_distance_to_plane(\
                 pVertex_start.dLongitude_degree, pVertex_start.dLatitude_degree,\
                 pVertex_in.dLongitude_degree, pVertex_in.dLatitude_degree,\
                 pVertex_end.dLongitude_degree,pVertex_end.dLatitude_degree)
            #lons = [pVertex_start.dLongitude_degree,pVertex_in.dLongitude_degree,pVertex_end.dLongitude_degree]
            #lats = [pVertex_start.dLatitude_degree, pVertex_in.dLatitude_degree, pVertex_end.dLatitude_degree]
            #dArea = calculate_polygon_area(lons, lats)

            if  angle3deg > 178 and d3 < 1.0: #care
                iFlag = 1
                dDistance = d1
            else:
                iFlag = 0

        else:
            iFlag = 0

        return iFlag, dDistance, dDistance_plane

    def calculate_distance_to_vertex(self, pVertex_in):
        pVertex_start = self.pVertex_start
        pVertex_end = self.pVertex_end
        pVertex_out = None
        dDistance_min = 1.0E9

        d1 = pVertex_start.calculate_distance(pVertex_in)
        d2 = pVertex_end.calculate_distance(pVertex_in)
        d3 = d1 + d2 - self.dLength
        angle3deg = calculate_angle_betwen_vertex(\
                 pVertex_start.dLongitude_degree, pVertex_start.dLatitude_degree,\
                 pVertex_in.dLongitude_degree, pVertex_in.dLatitude_degree,\
                 pVertex_end.dLongitude_degree,pVertex_end.dLatitude_degree)

        if angle3deg > 90: #care
            dDistance_plane = calculate_distance_to_plane(\
                 pVertex_start.dLongitude_degree, pVertex_start.dLatitude_degree,\
                 pVertex_in.dLongitude_degree, pVertex_in.dLatitude_degree,\
                 pVertex_end.dLongitude_degree,pVertex_end.dLatitude_degree)

            if dDistance_plane < d1 and dDistance_plane < d2:
                dLongitude_intersect, dLatitude_intersect = calculate_intersect_on_great_circle(\
                    pVertex_start.dLongitude_degree, pVertex_start.dLatitude_degree,\
                 pVertex_in.dLongitude_degree, pVertex_in.dLatitude_degree,\
                 pVertex_end.dLongitude_degree,pVertex_end.dLatitude_degree)

                point= dict()
                point['dLongitude_degree'] = dLongitude_intersect
                point['dLatitude_degree'] = dLatitude_intersect
                pVertex_out=pyvertex(point)

                dDistance_min = pVertex_out.calculate_distance(pVertex_in)
                if dDistance_min < d1 and dDistance_min < d2:
                    pass
                else:
                    if d1 < d2:
                        pVertex_out = pVertex_start
                        dDistance_min = d1
                    else:
                        pVertex_out = pVertex_end
                        dDistance_min = d2
                    pass
            else:
                print('error')
                pass
        else:
            if d1 < d2:
                pVertex_out = pVertex_start
                dDistance_min = d1
            else:
                pVertex_out = pVertex_end
                dDistance_min = d2
            pass


        return dDistance_min, pVertex_out

    def calculate_buffer_zone_polygon(self, dRadius, nPoint = 36, sFilename_out=None, sFolder_out=None):
        # Create a geodesic object

        if self.dLength < dRadius * 2.0:
            aEdge = [self]
        else:
            aEdge = self.split_by_length(dRadius)

        nEdge = len(aEdge)
        aVertex_out = list()
        aVertex_center = list()
        aVertex_circle = list()
        aLongitude_degree = []
        aLatitude_degree = []
        aCircle = list()
        for i in range(nEdge):
            pEdge = aEdge[i]
            pVertex_start = pEdge.pVertex_start
            pVertex_end = pEdge.pVertex_end

            aVertex_center.append(pVertex_start)
            aVertex_center.append(pVertex_end) #might, but should be ok have duplicate
            # Calculate the geodesic buffer

            aVertex_start = pVertex_start.calculate_buffer_zone_circle(dRadius, nPoint)
            #create a circle object
            pEdge.pCircle_start = pycircle(pVertex_start, aVertex_start)
            aVertex_circle.extend(aVertex_start)
            for pVertex_buffer0 in aVertex_start:
                aLongitude_degree.append(pVertex_buffer0.dLongitude_degree)
                aLatitude_degree.append(pVertex_buffer0.dLatitude_degree)

            aVertex_end = pVertex_end.calculate_buffer_zone_circle(dRadius, nPoint)
            pEdge.pCircle_end = pycircle(pVertex_end, aVertex_end)
            aVertex_circle.extend(aVertex_end)
            for pVertex_buffer1 in aVertex_end:
                aLongitude_degree.append(pVertex_buffer1.dLongitude_degree)
                aLatitude_degree.append(pVertex_buffer1.dLatitude_degree)

            if sFolder_out is not None:
                sFilename_dummy = os.path.join(sFolder_out, 'buffer_zone_start_%d.geojson' % i)
                export_vertex_as_polygon(aVertex_start, sFilename_dummy)
                sFilename_dummy = os.path.join(sFolder_out, 'buffer_zone_end_%d.geojson'% i)
                export_vertex_as_polygon(aVertex_end, sFilename_dummy)

            aCircle.append(pEdge.pCircle_start)
            aCircle.append(pEdge.pCircle_end)

        #combine them using convex function
        pPolygon_out = find_minimal_enclosing_polygon(aLongitude_degree, aLatitude_degree)
        #return as a list of vectex is more friendly

        for p in pPolygon_out:
            point0= dict()
            point0['dLongitude_degree'] = p[0]
            point0['dLatitude_degree'] = p[1]
            pVertex_out = pyvertex(point0)
            aVertex_out.append(pVertex_out)
        if sFilename_out is not None:
            export_vertex_as_polygon(aVertex_out, sFilename_out)

        return aVertex_out, aVertex_center, aVertex_circle, aCircle

    def __eq__(self, other):
        """
        Check if two edges are equivalent

        Args:
            other (pyedge): The other edge
            how about direction?

        Returns:
            int: 1 if equivalent; 0 if not
        """
        return int( self.pVertex_start == other.pVertex_start and self.pVertex_end == other.pVertex_end )

    def __ne__(self, other):
        """
        Check if two edges are equivalent

        Args:
            other (pyedge): The other edge

        Returns:
            int: 0 if equivalent; 1 if not
        """
        return not self.__eq__(other)

    def tojson(self):
        """
        Convert an edge object to a json string

        Returns:
            json str: A json string
        """

        obj = self.__dict__.copy()

        sJson = json.dumps(obj, \
            sort_keys=True, \
                indent = 4, \
                    ensure_ascii=True, \
                        cls=EdgeClassEncoder)
        return sJson

