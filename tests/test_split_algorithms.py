"""Unit tests for split algorithms in pyflowline and pyearthriver."""

import unittest

from pyearthriver.algorithms.topology.segmentation import (
    split_edge_by_length,
    split_flowline_by_length,
    split_flowline_to_edge,
)
from pyearthriver.classes.edge import pyedge
from pyearthriver.classes.flowline import pyflowline
from pyearthriver.classes.rivergraph import pyrivergraph
from pyearthriver.classes.vertex import pyvertex


class TestSplitAlgorithms(unittest.TestCase):
    """Test suite for flowline and edge split algorithms."""

    def setUp(self):
        self.v1 = pyvertex({"dLongitude_degree": -77.0, "dLatitude_degree": 40.0})
        self.v2 = pyvertex({"dLongitude_degree": -76.5, "dLatitude_degree": 40.0})
        self.v3 = pyvertex({"dLongitude_degree": -76.0, "dLatitude_degree": 40.0})
        self.edge1 = pyedge(self.v1, self.v2)
        self.edge1.set_edge_id(1)
        self.edge2 = pyedge(self.v2, self.v3)
        self.edge2.set_edge_id(2)

    def test_split_edge_short(self):
        """Edge shorter than threshold should not be split."""
        segments = split_edge_by_length(self.edge1, 100000.0)
        self.assertEqual(len(segments), 1)
        self.assertIs(segments[0], self.edge1)

    def test_split_edge_long(self):
        """Edge longer than threshold should be split using pyearthriver spherical slerp."""
        dLength_total = self.edge1.dLength  # ~42.6 km
        threshold = 15000.0  # 15 km
        segments = split_edge_by_length(self.edge1, threshold)

        self.assertGreater(len(segments), 1)
        total_len = sum(seg.dLength for seg in segments)
        self.assertAlmostEqual(total_len, dLength_total, delta=1.0)

        for seg in segments:
            self.assertLessEqual(seg.dLength, threshold * 1.01)
            self.assertIsInstance(seg.pVertex_start, pyvertex)
            self.assertIsInstance(seg.pVertex_end, pyvertex)
            self.assertEqual(seg.lEdgeID, 1)

    def test_split_edge_invalid_distance(self):
        """Negative or zero threshold should raise ValueError."""
        with self.assertRaises(ValueError):
            split_edge_by_length(self.edge1, 0.0)
        with self.assertRaises(ValueError):
            split_edge_by_length(self.edge1, -500.0)

    def test_split_flowline_short(self):
        """Flowline shorter than threshold should not be modified."""
        fl = pyflowline([self.edge1])
        fl.set_flowline_id(10)
        result = split_flowline_by_length([fl], 100000.0)
        self.assertEqual(len(result), 1)
        self.assertIs(result[0], fl)

    def test_split_flowline_multi_edge(self):
        """Multi-edge flowline longer than threshold should be split properly."""
        fl = pyflowline([self.edge1, self.edge2])
        fl.set_flowline_id(42)
        fl.set_stream_order(3)

        threshold = 20000.0
        result = split_flowline_by_length([fl], threshold)

        self.assertGreater(len(result), 2)
        total_len = sum(seg.dLength for seg in result)
        self.assertAlmostEqual(total_len, fl.dLength, delta=1.0)

        for seg in result:
            self.assertLessEqual(seg.dLength, threshold * 1.01)
            self.assertEqual(seg.lFlowlineID, 42)
            self.assertEqual(seg.iStream_order, 3)

    def test_split_flowline_to_edge(self):
        """Flowline should be split into individual 1-edge flowlines with attributes preserved."""
        fl = pyflowline([self.edge1, self.edge2])
        fl.set_flowline_id(99)
        fl.set_stream_order(4)

        flowlines_out, edges_out = split_flowline_to_edge([fl])
        self.assertEqual(len(flowlines_out), 2)
        self.assertEqual(len(edges_out), 2)

        for fl_seg in flowlines_out:
            self.assertEqual(fl_seg.nEdge, 1)
            self.assertEqual(fl_seg.lFlowlineID, 99)
            self.assertEqual(fl_seg.iStream_order, 4)

    def test_rivergraph_split_flowline_by_length(self):
        """pRivergraph.split_flowline_by_length should update both output and graph state."""
        fl = pyflowline([self.edge1, self.edge2])
        fl.set_flowline_id(55)
        fl.set_stream_order(2)

        rg = pyrivergraph([fl])
        self.assertEqual(len(rg.aFlowline), 1)

        threshold = 20000.0
        result = rg.split_flowline_by_length(threshold)

        self.assertGreater(len(result), 2)
        self.assertEqual(len(rg.aFlowline), len(result))
        for fl_seg in rg.aFlowline:
            self.assertLessEqual(fl_seg.dLength, threshold * 1.01)
            self.assertEqual(fl_seg.lFlowlineID, 55)
            self.assertEqual(fl_seg.iStream_order, 2)

    def test_rivergraph_split_flowline_to_edge(self):
        """pRivergraph.split_flowline_to_edge should split into 1-edge flowlines and update graph."""
        fl = pyflowline([self.edge1, self.edge2])
        fl.set_flowline_id(77)
        fl.set_stream_order(3)

        rg = pyrivergraph([fl])
        fl_out, edges_out = rg.split_flowline_to_edge()

        self.assertEqual(len(fl_out), 2)
        self.assertEqual(len(edges_out), 2)
        self.assertEqual(len(rg.aFlowline), 2)
        for fl_seg in rg.aFlowline:
            self.assertEqual(fl_seg.nEdge, 1)
            self.assertEqual(fl_seg.lFlowlineID, 77)
            self.assertEqual(fl_seg.iStream_order, 3)


if __name__ == "__main__":
    unittest.main()
