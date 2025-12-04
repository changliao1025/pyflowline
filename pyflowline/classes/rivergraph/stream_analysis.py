"""
Stream analysis module for pyrivergraph.

This module provides stream ordering algorithms, hydrological calculations,
and stream network analysis capabilities including Strahler and Shreve
ordering methods, flow accumulation, and watershed delineation.
"""

import time
import logging
import math
from typing import List, Dict, Any, Optional, Set, Tuple, Union
from collections import defaultdict, deque
from enum import Enum

from .state_management import GraphStateManager
from .utils import calculate_distance_2d, topological_sort, calculate_flow_accumulation

logger = logging.getLogger(__name__)


class StreamOrderMethod(Enum):
    """Available stream ordering methods."""
    STRAHLER = "strahler"
    SHREVE = "shreve"
    HORTON = "horton"


class StreamClassification(Enum):
    """Stream classification categories."""
    HEADWATER = "headwater"
    TRIBUTARY = "tributary"
    MAIN_STEM = "main_stem"
    OUTLET = "outlet"


class StreamAnalysisResult:
    """Container for stream analysis results."""

    def __init__(self):
        self.stream_orders = {}
        self.classifications = {}
        self.flow_accumulation = {}
        self.watershed_properties = {}
        self.execution_time = 0.0
        self.method_used = None

    def add_stream_order(self, flowline_id: int, order: int):
        """Add stream order for a flowline."""
        self.stream_orders[flowline_id] = order

    def add_classification(self, flowline_id: int, classification: StreamClassification):
        """Add classification for a flowline."""
        self.classifications[flowline_id] = classification

    def add_flow_accumulation(self, flowline_id: int, accumulation: float):
        """Add flow accumulation for a flowline."""
        self.flow_accumulation[flowline_id] = accumulation

    def get_summary(self) -> Dict[str, Any]:
        """Get a summary of stream analysis results."""
        order_distribution = {}
        if self.stream_orders:
            for order in self.stream_orders.values():
                order_distribution[order] = order_distribution.get(order, 0) + 1

        return {
            'total_streams': len(self.stream_orders),
            'max_order': max(self.stream_orders.values()) if self.stream_orders else 0,
            'order_distribution': order_distribution,
            'method_used': self.method_used.value if self.method_used else None,
            'execution_time': self.execution_time
        }


class StreamAnalyzer:
    """
    Comprehensive stream analyzer for river graphs with stream ordering,
    classification, and hydrological analysis capabilities.
    """

    def __init__(self, graph_instance: Any, state_manager: Optional[GraphStateManager] = None):
        """
        Initialize the stream analyzer.

        Args:
            graph_instance: The pyrivergraph instance to analyze
            state_manager: Optional state manager for tracking analysis history
        """
        self.graph = graph_instance
        self.state_manager = state_manager
        self.tolerance = 1e-6  # Distance tolerance for geometric operations

        # Cached analysis results
        self._stream_order_cache = {}
        self._flow_accumulation_cache = {}
        self._topology_cache = None

        logger.debug("Initialized StreamAnalyzer")

    def calculate_stream_order(self,
                             method: StreamOrderMethod = StreamOrderMethod.STRAHLER,
                             force_recalculate: bool = False) -> StreamAnalysisResult:
        """
        Calculate stream order using the specified method.

        Args:
            method: Stream ordering method to use
            force_recalculate: Whether to force recalculation even if cached

        Returns:
            StreamAnalysisResult with stream orders and analysis
        """
        start_time = time.time()
        result = StreamAnalysisResult()
        result.method_used = method

        try:
            # Check cache first
            cache_key = method.value
            if not force_recalculate and cache_key in self._stream_order_cache:
                cached_result = self._stream_order_cache[cache_key]
                result.stream_orders = cached_result.copy()
                result.execution_time = time.time() - start_time
                logger.debug(f"Using cached stream orders for method {method.value}")
                return result

            logger.info(f"Calculating stream order using {method.value} method")

            if method == StreamOrderMethod.STRAHLER:
                orders = self._calculate_strahler_order()
            elif method == StreamOrderMethod.SHREVE:
                orders = self._calculate_shreve_order()
            elif method == StreamOrderMethod.HORTON:
                orders = self._calculate_horton_order()
            else:
                raise ValueError(f"Unknown stream ordering method: {method}")

            # Store results
            for flowline_id, order in orders.items():
                result.add_stream_order(flowline_id, order)

            # Cache results
            self._stream_order_cache[cache_key] = orders.copy()

            # Update flowline objects with stream orders
            self._update_flowline_stream_orders(orders)

            result.execution_time = time.time() - start_time
            logger.info(f"Stream order calculation completed in {result.execution_time:.3f}s")

            return result

        except Exception as e:
            logger.error(f"Stream order calculation failed: {e}")
            result.execution_time = time.time() - start_time
            return result

    def classify_streams(self) -> Dict[int, StreamClassification]:
        """
        Classify streams based on their position and characteristics.

        Returns:
            Dictionary mapping flowline IDs to classifications
        """
        try:
            classifications = {}

            # Build topology if not cached
            if self._topology_cache is None:
                self._build_topology_cache()

            for i, flowline in enumerate(self.graph.aFlowline):
                if flowline is None:
                    continue

                classification = self._classify_single_stream(i)
                classifications[i] = classification

            return classifications

        except Exception as e:
            logger.error(f"Stream classification failed: {e}")
            return {}


        """
        Analyze watershed properties including drainage density, bifurcation ratios, etc.

        Returns:
            Dictionary with watershed analysis results
        """
        try:
            properties = {}

            # Calculate basic metrics
            total_length = self._calculate_total_stream_length()
            properties['total_stream_length'] = total_length

            # Calculate drainage density
            drainage_area = self._calculate_total_drainage_area()
            if drainage_area > 0:
                properties['drainage_density'] = total_length / drainage_area
            else:
                properties['drainage_density'] = 0.0

            # Calculate bifurcation ratios
            bifurcation_ratios = self._calculate_bifurcation_ratios()
            properties['bifurcation_ratios'] = bifurcation_ratios

            # Calculate length ratios
            length_ratios = self._calculate_length_ratios()
            properties['length_ratios'] = length_ratios

            # Calculate stream frequency
            if drainage_area > 0:
                stream_count = len([f for f in self.graph.aFlowline if f is not None])
                properties['stream_frequency'] = stream_count / drainage_area
            else:
                properties['stream_frequency'] = 0.0

            # Calculate network properties
            properties['network_properties'] = self._analyze_network_properties()

            return properties

        except Exception as e:
            logger.error(f"Watershed property analysis failed: {e}")
            return {}