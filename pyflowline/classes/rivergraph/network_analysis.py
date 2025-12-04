"""
Network analysis module for pyrivergraph.

This module provides topology analysis and pattern detection capabilities
for river network graphs, including connectivity analysis, network metrics,
and topological pattern identification.
"""

import time
import logging
import math
from typing import List, Dict, Any, Optional, Set, Tuple, Union
from collections import defaultdict, deque
from enum import Enum

from .state_management import GraphStateManager
from .utils import calculate_distance_2d, topological_sort

logger = logging.getLogger(__name__)


class NetworkMetric(Enum):
    """Available network metrics for analysis."""
    CONNECTIVITY = "connectivity"
    CENTRALITY = "centrality"
    CLUSTERING = "clustering"
    DIAMETER = "diameter"
    DENSITY = "density"
    EFFICIENCY = "efficiency"


class TopologyPattern(Enum):
    """Topological patterns in river networks."""
    DENDRITIC = "dendritic"
    PARALLEL = "parallel"
    TRELLIS = "trellis"
    RECTANGULAR = "rectangular"
    RADIAL = "radial"
    ANNULAR = "annular"


class NetworkAnalysisResult:
    """Container for network analysis results."""

    def __init__(self):
        self.metrics = {}
        self.patterns = {}
        self.statistics = {}
        self.execution_time = 0.0
        self.analysis_performed = []

    def add_metric(self, metric_name: str, value: Any, description: str = ""):
        """Add a metric to the results."""
        self.metrics[metric_name] = {
            'value': value,
            'description': description
        }

    def add_pattern(self, pattern_name: str, confidence: float, details: Dict[str, Any] = None):
        """Add a detected pattern to the results."""
        self.patterns[pattern_name] = {
            'confidence': confidence,
            'details': details or {}
        }

    def add_statistic(self, stat_name: str, value: Any):
        """Add a statistic to the results."""
        self.statistics[stat_name] = value

    def get_summary(self) -> Dict[str, Any]:
        """Get a summary of analysis results."""
        return {
            'metrics_count': len(self.metrics),
            'patterns_detected': len(self.patterns),
            'statistics_count': len(self.statistics),
            'execution_time': self.execution_time,
            'analysis_performed': self.analysis_performed
        }


class NetworkAnalyzer:
    """
    Comprehensive network analyzer for river graphs with topology analysis,
    pattern detection, and network metrics calculation.
    """

    def __init__(self, graph_instance: Any, state_manager: Optional[GraphStateManager] = None):
        """
        Initialize the network analyzer.

        Args:
            graph_instance: The pyrivergraph instance to analyze
            state_manager: Optional state manager for tracking analysis history
        """
        self.graph = graph_instance
        self.state_manager = state_manager
        self.tolerance = 1e-6  # Distance tolerance for geometric operations

        # Cached analysis results
        self._adjacency_cache = None
        self._connectivity_cache = None
        self._metrics_cache = {}

        logger.debug("Initialized NetworkAnalyzer")

    def analyze_network(self,
                       metrics: List[NetworkMetric] = None,
                       detect_patterns: bool = True,
                       calculate_statistics: bool = True) -> NetworkAnalysisResult:
        """
        Perform comprehensive network analysis.

        Args:
            metrics: List of metrics to calculate (None for all)
            detect_patterns: Whether to detect topological patterns
            calculate_statistics: Whether to calculate basic statistics

        Returns:
            NetworkAnalysisResult with analysis findings
        """
        start_time = time.time()
        result = NetworkAnalysisResult()

        try:
            logger.info("Starting comprehensive network analysis")

            # Calculate basic statistics
            if calculate_statistics:
                self._calculate_basic_statistics(result)
                result.analysis_performed.append("basic_statistics")

            # Calculate network metrics
            if metrics is None:
                metrics = list(NetworkMetric)

            for metric in metrics:
                self._calculate_metric(metric, result)
                result.analysis_performed.append(f"metric_{metric.value}")

            # Detect topological patterns
            if detect_patterns:
                self._detect_patterns(result)
                result.analysis_performed.append("pattern_detection")

            result.execution_time = time.time() - start_time

            # Log summary
            summary = result.get_summary()
            logger.info(f"Network analysis completed: {summary}")

            return result

        except Exception as e:
            logger.error(f"Network analysis failed: {e}")
            result.execution_time = time.time() - start_time
            return result

    def analyze_connectivity(self) -> Dict[str, Any]:
        """
        Analyze network connectivity properties.

        Returns:
            Dictionary with connectivity analysis results
        """
        try:
            if self._connectivity_cache is not None:
                return self._connectivity_cache

            connectivity = {}

            # Find connected components
            components = self._find_connected_components()
            connectivity['component_count'] = len(components)
            connectivity['largest_component_size'] = max(len(comp) for comp in components) if components else 0
            connectivity['component_sizes'] = [len(comp) for comp in components]

            # Calculate connectivity ratio
            total_flowlines = len([f for f in self.graph.aFlowline if f is not None])
            connectivity['connectivity_ratio'] = connectivity['largest_component_size'] / total_flowlines if total_flowlines > 0 else 0

            # Find articulation points (critical flowlines)
            articulation_points = self._find_articulation_points()
            connectivity['articulation_points'] = articulation_points
            connectivity['articulation_count'] = len(articulation_points)

            # Find bridges (critical connections)
            bridges = self._find_bridges()
            connectivity['bridges'] = bridges
            connectivity['bridge_count'] = len(bridges)

            self._connectivity_cache = connectivity
            return connectivity

        except Exception as e:
            logger.error(f"Connectivity analysis failed: {e}")
            return {}

    def analyze_centrality(self) -> Dict[str, Any]:
        """
        Analyze centrality measures for network elements.

        Returns:
            Dictionary with centrality analysis results
        """
        try:
            centrality = {}

            # Build adjacency representation
            adjacency = self._get_adjacency_list()

            # Calculate degree centrality
            degree_centrality = self._calculate_degree_centrality(adjacency)
            centrality['degree_centrality'] = degree_centrality

            # Calculate betweenness centrality
            betweenness_centrality = self._calculate_betweenness_centrality(adjacency)
            centrality['betweenness_centrality'] = betweenness_centrality

            # Calculate closeness centrality
            closeness_centrality = self._calculate_closeness_centrality(adjacency)
            centrality['closeness_centrality'] = closeness_centrality

            # Find most central nodes
            centrality['most_central_degree'] = max(degree_centrality.items(), key=lambda x: x[1]) if degree_centrality else None
            centrality['most_central_betweenness'] = max(betweenness_centrality.items(), key=lambda x: x[1]) if betweenness_centrality else None
            centrality['most_central_closeness'] = max(closeness_centrality.items(), key=lambda x: x[1]) if closeness_centrality else None

            return centrality

        except Exception as e:
            logger.error(f"Centrality analysis failed: {e}")
            return {}

    def analyze_flow_patterns(self) -> Dict[str, Any]:
        """
        Analyze flow patterns in the network.

        Returns:
            Dictionary with flow pattern analysis results
        """
        try:
            flow_patterns = {}

            # Analyze flow accumulation
            flow_accumulation = self._calculate_flow_accumulation()
            flow_patterns['flow_accumulation'] = flow_accumulation

            # Find main stems (highest flow accumulation paths)
            main_stems = self._find_main_stems(flow_accumulation)
            flow_patterns['main_stems'] = main_stems

            # Analyze tributary patterns
            tributary_patterns = self._analyze_tributary_patterns()
            flow_patterns['tributary_patterns'] = tributary_patterns

            # Calculate drainage density
            drainage_density = self._calculate_drainage_density()
            flow_patterns['drainage_density'] = drainage_density

            # Analyze bifurcation ratios
            bifurcation_ratios = self._calculate_bifurcation_ratios()
            flow_patterns['bifurcation_ratios'] = bifurcation_ratios

            return flow_patterns

        except Exception as e:
            logger.error(f"Flow pattern analysis failed: {e}")
            return {}

    def detect_network_anomalies(self) -> List[Dict[str, Any]]:
        """
        Detect anomalies in the network structure.

        Returns:
            List of detected anomalies with details
        """
        anomalies = []

        try:
            # Detect isolated flowlines
            isolated = self._detect_isolated_flowlines()
            if isolated:
                anomalies.append({
                    'type': 'isolated_flowlines',
                    'count': len(isolated),
                    'flowlines': isolated,
                    'severity': 'medium'
                })

            # Detect dangling flowlines
            dangling = self._detect_dangling_flowlines()
            if dangling:
                anomalies.append({
                    'type': 'dangling_flowlines',
                    'count': len(dangling),
                    'flowlines': dangling,
                    'severity': 'low'
                })

            # Detect unusual confluence patterns
            unusual_confluences = self._detect_unusual_confluences()
            if unusual_confluences:
                anomalies.append({
                    'type': 'unusual_confluences',
                    'count': len(unusual_confluences),
                    'confluences': unusual_confluences,
                    'severity': 'medium'
                })

            # Detect flow direction inconsistencies
            flow_inconsistencies = self._detect_flow_inconsistencies()
            if flow_inconsistencies:
                anomalies.append({
                    'type': 'flow_inconsistencies',
                    'count': len(flow_inconsistencies),
                    'flowlines': flow_inconsistencies,
                    'severity': 'high'
                })

            return anomalies

        except Exception as e:
            logger.error(f"Anomaly detection failed: {e}")
            return []

    # Private analysis methods

    def _calculate_basic_statistics(self, result: NetworkAnalysisResult):
        """Calculate basic network statistics."""
        try:
            # Count elements
            flowline_count = len([f for f in self.graph.aFlowline if f is not None])
            vertex_count = len([v for v in self.graph.aVertex if v is not None])
            confluence_count = len([c for c in self.graph.aConfluence if c is not None])

            result.add_statistic('flowline_count', flowline_count)
            result.add_statistic('vertex_count', vertex_count)
            result.add_statistic('confluence_count', confluence_count)

            # Calculate total length
            total_length = 0
            for flowline in self.graph.aFlowline:
                if flowline is not None and hasattr(flowline, 'dLength'):
                    total_length += flowline.dLength

            result.add_statistic('total_length', total_length)
            result.add_statistic('average_flowline_length', total_length / flowline_count if flowline_count > 0 else 0)

            # Calculate stream order statistics
            stream_orders = []
            for flowline in self.graph.aFlowline:
                if flowline is not None and hasattr(flowline, 'iStreamOrder'):
                    stream_orders.append(flowline.iStreamOrder)

            if stream_orders:
                result.add_statistic('max_stream_order', max(stream_orders))
                result.add_statistic('average_stream_order', sum(stream_orders) / len(stream_orders))
                result.add_statistic('stream_order_distribution', {order: stream_orders.count(order) for order in set(stream_orders)})

        except Exception as e:
            logger.error(f"Basic statistics calculation failed: {e}")

    def _calculate_metric(self, metric: NetworkMetric, result: NetworkAnalysisResult):
        """Calculate a specific network metric."""
        try:
            if metric == NetworkMetric.CONNECTIVITY:
                connectivity = self.analyze_connectivity()
                result.add_metric('connectivity_ratio', connectivity.get('connectivity_ratio', 0),
                                "Ratio of largest connected component to total network")
                result.add_metric('component_count', connectivity.get('component_count', 0),
                                "Number of disconnected components")

            elif metric == NetworkMetric.CENTRALITY:
                centrality = self.analyze_centrality()
                if centrality.get('most_central_degree'):
                    result.add_metric('most_central_node', centrality['most_central_degree'][0],
                                    "Node with highest degree centrality")

            elif metric == NetworkMetric.DENSITY:
                density = self._calculate_network_density()
                result.add_metric('network_density', density,
                                "Ratio of actual connections to possible connections")

            elif metric == NetworkMetric.DIAMETER:
                diameter = self._calculate_network_diameter()
                result.add_metric('network_diameter', diameter,
                                "Longest shortest path in the network")

            elif metric == NetworkMetric.EFFICIENCY:
                efficiency = self._calculate_network_efficiency()
                result.add_metric('network_efficiency', efficiency,
                                "Average inverse shortest path length")

        except Exception as e:
            logger.error(f"Metric calculation failed for {metric}: {e}")

    def _detect_patterns(self, result: NetworkAnalysisResult):
        """Detect topological patterns in the network."""
        try:
            # Detect dendritic pattern
            dendritic_confidence = self._detect_dendritic_pattern()
            if dendritic_confidence > 0.5:
                result.add_pattern('dendritic', dendritic_confidence,
                                 {'description': 'Tree-like branching pattern'})

            # Detect parallel pattern
            parallel_confidence = self._detect_parallel_pattern()
            if parallel_confidence > 0.3:
                result.add_pattern('parallel', parallel_confidence,
                                 {'description': 'Parallel drainage channels'})

            # Detect trellis pattern
            trellis_confidence = self._detect_trellis_pattern()
            if trellis_confidence > 0.3:
                result.add_pattern('trellis', trellis_confidence,
                                 {'description': 'Rectangular grid-like pattern'})

        except Exception as e:
            logger.error(f"Pattern detection failed: {e}")

    def _get_adjacency_list(self) -> Dict[int, List[int]]:
        """Build adjacency list representation of the network."""
        if self._adjacency_cache is not None:
            return self._adjacency_cache

        adjacency = defaultdict(list)

        try:
            for i, flowline in enumerate(self.graph.aFlowline):
                if flowline is None:
                    continue

                if hasattr(flowline, 'aVertex') and flowline.aVertex and len(flowline.aVertex) >= 2:
                    start_vertex = flowline.aVertex[0]
                    end_vertex = flowline.aVertex[-1]

                    if (hasattr(start_vertex, 'lVertexIndex') and
                        hasattr(end_vertex, 'lVertexIndex')):
                        start_id = start_vertex.lVertexIndex
                        end_id = end_vertex.lVertexIndex

                        adjacency[start_id].append(end_id)
                        adjacency[end_id].append(start_id)  # Undirected for connectivity analysis

            self._adjacency_cache = dict(adjacency)
            return self._adjacency_cache

        except Exception as e:
            logger.error(f"Adjacency list construction failed: {e}")
            return {}

    def _find_connected_components(self) -> List[Set[int]]:
        """Find connected components using DFS."""
        adjacency = self._get_adjacency_list()
        visited = set()
        components = []

        def dfs(node: int, component: Set[int]):
            if node in visited:
                return
            visited.add(node)
            component.add(node)

            for neighbor in adjacency.get(node, []):
                if neighbor not in visited:
                    dfs(neighbor, component)

        for node in adjacency:
            if node not in visited:
                component = set()
                dfs(node, component)
                if component:
                    components.append(component)

        return components

    def _find_articulation_points(self) -> List[int]:
        """Find articulation points (cut vertices) in the network."""
        adjacency = self._get_adjacency_list()
        if not adjacency:
            return []

        visited = set()
        discovery = {}
        low = {}
        parent = {}
        articulation_points = set()
        time = [0]  # Use list to allow modification in nested function

        def dfs(u: int):
            children = 0
            visited.add(u)
            discovery[u] = low[u] = time[0]
            time[0] += 1

            for v in adjacency.get(u, []):
                if v not in visited:
                    children += 1
                    parent[v] = u
                    dfs(v)

                    low[u] = min(low[u], low[v])

                    # Root is articulation point if it has more than one child
                    if parent.get(u) is None and children > 1:
                        articulation_points.add(u)

                    # Non-root is articulation point if removing it disconnects the graph
                    if parent.get(u) is not None and low[v] >= discovery[u]:
                        articulation_points.add(u)

                elif v != parent.get(u):
                    low[u] = min(low[u], discovery[v])

        for node in adjacency:
            if node not in visited:
                dfs(node)

        return list(articulation_points)

    def _find_bridges(self) -> List[Tuple[int, int]]:
        """Find bridges (cut edges) in the network."""
        adjacency = self._get_adjacency_list()
        if not adjacency:
            return []

        visited = set()
        discovery = {}
        low = {}
        parent = {}
        bridges = []
        time = [0]

        def dfs(u: int):
            visited.add(u)
            discovery[u] = low[u] = time[0]
            time[0] += 1

            for v in adjacency.get(u, []):
                if v not in visited:
                    parent[v] = u
                    dfs(v)

                    low[u] = min(low[u], low[v])

                    # Edge (u, v) is a bridge if low[v] > discovery[u]
                    if low[v] > discovery[u]:
                        bridges.append((u, v))

                elif v != parent.get(u):
                    low[u] = min(low[u], discovery[v])

        for node in adjacency:
            if node not in visited:
                dfs(node)

        return bridges

    def _calculate_degree_centrality(self, adjacency: Dict[int, List[int]]) -> Dict[int, float]:
        """Calculate degree centrality for all nodes."""
        degree_centrality = {}
        max_degree = len(adjacency) - 1 if len(adjacency) > 1 else 1

        for node in adjacency:
            degree = len(adjacency[node])
            degree_centrality[node] = degree / max_degree

        return degree_centrality

    def _calculate_betweenness_centrality(self, adjacency: Dict[int, List[int]]) -> Dict[int, float]:
        """Calculate betweenness centrality using Brandes' algorithm."""
        betweenness = {node: 0.0 for node in adjacency}

        for s in adjacency:
            # Single-source shortest paths
            stack = []
            paths = {node: [] for node in adjacency}
            sigma = {node: 0 for node in adjacency}
            sigma[s] = 1
            distance = {node: -1 for node in adjacency}
            distance[s] = 0
            queue = deque([s])

            while queue:
                v = queue.popleft()
                stack.append(v)

                for w in adjacency.get(v, []):
                    if distance[w] < 0:
                        queue.append(w)
                        distance[w] = distance[v] + 1

                    if distance[w] == distance[v] + 1:
                        sigma[w] += sigma[v]
                        paths[w].append(v)

            # Accumulation
            delta = {node: 0 for node in adjacency}
            while stack:
                w = stack.pop()
                for v in paths[w]:
                    delta[v] += (sigma[v] / sigma[w]) * (1 + delta[w])
                if w != s:
                    betweenness[w] += delta[w]

        # Normalize
        n = len(adjacency)
        if n > 2:
            normalization = 2.0 / ((n - 1) * (n - 2))
            for node in betweenness:
                betweenness[node] *= normalization

        return betweenness

    def _calculate_closeness_centrality(self, adjacency: Dict[int, List[int]]) -> Dict[int, float]:
        """Calculate closeness centrality for all nodes."""
        closeness = {}

        for node in adjacency:
            # Calculate shortest paths from this node to all others
            distances = self._shortest_paths_from_node(node, adjacency)

            if distances:
                total_distance = sum(distances.values())
                if total_distance > 0:
                    closeness[node] = (len(distances) - 1) / total_distance
                else:
                    closeness[node] = 0.0
            else:
                closeness[node] = 0.0

        return closeness

    def _shortest_paths_from_node(self, start: int, adjacency: Dict[int, List[int]]) -> Dict[int, int]:
        """Calculate shortest paths from a node to all other nodes using BFS."""
        distances = {start: 0}
        queue = deque([start])

        while queue:
            current = queue.popleft()

            for neighbor in adjacency.get(current, []):
                if neighbor not in distances:
                    distances[neighbor] = distances[current] + 1
                    queue.append(neighbor)

        return distances

    def _calculate_network_density(self) -> float:
        """Calculate network density."""
        adjacency = self._get_adjacency_list()
        if len(adjacency) < 2:
            return 0.0

        actual_edges = sum(len(neighbors) for neighbors in adjacency.values()) // 2
        possible_edges = len(adjacency) * (len(adjacency) - 1) // 2

        return actual_edges / possible_edges if possible_edges > 0 else 0.0

    def _calculate_network_diameter(self) -> int:
        """Calculate network diameter (longest shortest path)."""
        adjacency = self._get_adjacency_list()
        max_distance = 0

        for node in adjacency:
            distances = self._shortest_paths_from_node(node, adjacency)
            if distances:
                max_distance = max(max_distance, max(distances.values()))

        return max_distance

    def _calculate_network_efficiency(self) -> float:
        """Calculate network efficiency."""
        adjacency = self._get_adjacency_list()
        if len(adjacency) < 2:
            return 0.0

        total_efficiency = 0.0
        pair_count = 0

        for node in adjacency:
            distances = self._shortest_paths_from_node(node, adjacency)
            for other_node, distance in distances.items():
                if other_node != node and distance > 0:
                    total_efficiency += 1.0 / distance
                    pair_count += 1

        return total_efficiency / pair_count if pair_count > 0 else 0.0

    def _calculate_flow_accumulation(self) -> Dict[int, float]:
        """Calculate flow accumulation for each flowline."""
        flow_accumulation = {}

        try:
            # Simple flow accumulation based on upstream area
            for i, flowline in enumerate(self.graph.aFlowline):
                if flowline is None:
                    continue

                # Start with base flow (could be based on drainage area)
                base_flow = 1.0
                if hasattr(flowline, 'dDrainageArea'):
                    base_flow = flowline.dDrainageArea

                # Add upstream contributions
                upstream_flow = self._calculate_upstream_flow(i)
                flow_accumulation[i] = base_flow + upstream_flow

            return flow_accumulation

        except Exception as e:
            logger.error(f"Flow accumulation calculation failed: {e}")
            return {}

    def _calculate_upstream_flow(self, flowline_id: int) -> float:
        """Calculate upstream flow contribution for a flowline."""
        upstream_flow = 0.0

        try:
            flowline = self.graph.aFlowline[flowline_id]
            if flowline is None or not hasattr(flowline, 'aVertex') or not flowline.aVertex:
                return 0.0

            # Get start vertex
            start_vertex = flowline.aVertex[0]
            if not hasattr(start_vertex, 'lVertexIndex'):
                return 0.0

            start_vertex_id = start_vertex.lVertexIndex

            # Find upstream flowlines
            for i, other_flowline in enumerate(self.graph.aFlowline):
                if i == flowline_id or other_flowline is None:
                    continue

                if hasattr(other_flowline, 'aVertex') and other_flowline.aVertex:
                    end_vertex = other_flowline.aVertex[-1]
                    if (hasattr(end_vertex, 'lVertexIndex') and
                        end_vertex.lVertexIndex == start_vertex_id):
                        # This flowline flows into our start vertex
                        if hasattr(other_flowline, 'dDrainageArea'):
                            upstream_flow += other_flowline.dDrainageArea
                        else:
                            upstream_flow += 1.0

            return upstream_flow

        except Exception as e:
            logger.error(f"Upstream flow calculation failed: {e}")
            return 0.0

    def _find_main_stems(self, flow_accumulation: Dict[int, float]) -> List[int]:
        """Find main stem flowlines based on flow accumulation."""
        if not flow_accumulation:
            return []

        # Sort by flow accumulation
        sorted_flowlines = sorted(flow_accumulation.items(), key=lambda x: x[1], reverse=True)

        # Take top 10% as main stems
        main_stem_count = max(1, len(sorted_flowlines) // 10)
        return [flowline_id for flowline_id, _ in sorted_flowlines[:main_stem_count]]

    def _analyze_tributary_patterns(self) -> Dict[str, Any]:
        """Analyze tributary junction patterns."""
        patterns = {
            'junction_angles': [],
            'tributary_lengths': [],
            'junction_types': defaultdict(int)
        }

        try:
            for i, confluence in enumerate(self.graph.aConfluence):
                if confluence is None:
                    continue

                if hasattr(confluence, 'aFlowline') and confluence.aFlowline:
                    junction_flowlines = len(confluence.aFlowline)

                    if junction_flowlines == 2:
                        patterns['junction_types']['simple'] += 1
                    elif junction_flowlines == 3:
                        patterns['junction_types']['triple'] += 1
                    elif junction_flowlines > 3:
                        patterns['junction_types']['complex'] += 1

                    # Calculate junction angles if possible
                    if junction_flowlines >= 2:
                        angles = self._calculate_junction_angles(confluence)
                        patterns['junction_angles'].extend(angles)

            return patterns

        except Exception as e:
            logger.error(f"Tributary pattern analysis failed: {e}")
            return patterns

    def _calculate_junction_angles(self, confluence: Any) -> List[float]:
        """Calculate angles at a confluence junction."""
        angles = []

        try:
            if not hasattr(confluence, 'aFlowline') or len(confluence.aFlowline) < 2:
                return angles

            # Get directions of flowlines at confluence
            directions = []
            for flowline in confluence.aFlowline:
                if flowline and hasattr(flowline, 'aVertex') and flowline.aVertex:
                    # Calculate direction based on last few vertices
                    if len(flowline.aVertex) >= 2:
                        v1 = flowline.aVertex[-2]
                        v2 = flowline.aVertex[-1]

                        if (hasattr(v1, 'dX') and hasattr(v1, 'dY') and
                            hasattr(v2, 'dX') and hasattr(v2, 'dY')):
                            dx = v2.dX - v1.dX
                            dy = v2.dY - v1.dY
                            angle = math.atan2(dy, dx)
                            directions.append(angle)

            # Calculate angles between directions
            for i in range(len(directions)):
                for j in range(i + 1, len(directions)):
                    angle_diff = abs(directions[i] - directions[j])
                    angle_diff = min(angle_diff, 2 * math.pi - angle_diff)
                    angles.append(math.degrees(angle_diff))

            return angles

        except Exception as e:
            logger.error(f"Junction angle calculation failed: {e}")
            return angles

    def _calculate_drainage_density(self) -> float:
        """Calculate drainage density (total length / drainage area)."""
        try:
            total_length = 0.0
            total_area = 0.0

            for flowline in self.graph.aFlowline:
                if flowline is None:
                    continue

                if hasattr(flowline, 'dLength'):
                    total_length += flowline.dLength

                if hasattr(flowline, 'dDrainageArea'):
                    total_area += flowline.dDrainageArea

            return total_length / total_area if total_area > 0 else 0.0

        except Exception as e:
            logger.error(f"Drainage density calculation failed: {e}")
            return 0.0

    def _calculate_bifurcation_ratios(self) -> Dict[int, float]:
        """Calculate bifurcation ratios for different stream orders."""
        try:
            # Count streams by order
            stream_counts = defaultdict(int)

            for flowline in self.graph.aFlowline:
                if flowline is None:
                    continue

                if hasattr(flowline, 'iStreamOrder'):
                    stream_counts[flowline.iStreamOrder] += 1

            # Calculate bifurcation ratios
            bifurcation_ratios = {}
            orders = sorted(stream_counts.keys())

            for i in range(len(orders) - 1):
                current_order = orders[i]
                next_order = orders[i + 1]

                if stream_counts[next_order] > 0:
                    bifurcation_ratios[current_order] = stream_counts[current_order] / stream_counts[next_order]

            return bifurcation_ratios

        except Exception as e:
            logger.error(f"Bifurcation ratio calculation failed: {e}")
            return {}

    def _detect_dendritic_pattern(self) -> float:
        """Detect dendritic (tree-like) drainage pattern."""
        try:
            # Dendritic patterns have high bifurcation ratios and tree-like structure
            bifurcation_ratios = self._calculate_bifurcation_ratios()

            if not bifurcation_ratios:
                return 0.0

            # Check if bifurcation ratios are in typical dendritic range (3-5)
            dendritic_score = 0.0
            for ratio in bifurcation_ratios.values():
                if 3.0 <= ratio <= 5.0:
                    dendritic_score += 0.3
                elif 2.0 <= ratio <= 6.0:
                    dendritic_score += 0.1

            # Check for tree-like connectivity (no cycles)
            connectivity = self.analyze_connectivity()
            if connectivity.get('bridge_count', 0) == len([f for f in self.graph.aFlowline if f is not None]) - 1:
                dendritic_score += 0.4

            return min(dendritic_score, 1.0)

        except Exception as e:
            logger.error(f"Dendritic pattern detection failed: {e}")
            return 0.0

    def _detect_parallel_pattern(self) -> float:
        """Detect parallel drainage pattern."""
        try:
            # Parallel patterns have flowlines running in similar directions
            directions = []

            for flowline in self.graph.aFlowline:
                if flowline is None or not hasattr(flowline, 'aVertex') or not flowline.aVertex:
                    continue

                if len(flowline.aVertex) >= 2:
                    start = flowline.aVertex[0]
                    end = flowline.aVertex[-1]

                    if (hasattr(start, 'dX') and hasattr(start, 'dY') and
                        hasattr(end, 'dX') and hasattr(end, 'dY')):
                        dx = end.dX - start.dX
                        dy = end.dY - start.dY

                        if dx != 0 or dy != 0:
                            angle = math.atan2(dy, dx)
                            directions.append(angle)

            if len(directions) < 2:
                return 0.0

            # Calculate direction variance
            mean_direction = sum(directions) / len(directions)
            variance = sum((d - mean_direction) ** 2 for d in directions) / len(directions)

            # Lower variance indicates more parallel pattern
            parallel_score = max(0, 1.0 - variance / (math.pi / 4))

            return parallel_score

        except Exception as e:
            logger.error(f"Parallel pattern detection failed: {e}")
            return 0.0

    def _detect_trellis_pattern(self) -> float:
        """Detect trellis (rectangular grid-like) drainage pattern."""
        try:
            # Trellis patterns have perpendicular tributary junctions
            junction_angles = []

            for confluence in self.graph.aConfluence:
                if confluence is None:
                    continue

                angles = self._calculate_junction_angles(confluence)
                junction_angles.extend(angles)

            if not junction_angles:
                return 0.0

            # Count near-perpendicular angles (80-100 degrees)
            perpendicular_count = sum(1 for angle in junction_angles if 80 <= angle <= 100)
            trellis_score = perpendicular_count / len(junction_angles)

            return trellis_score

        except Exception as e:
            logger.error(f"Trellis pattern detection failed: {e}")
            return 0.0

    def _detect_isolated_flowlines(self) -> List[int]:
        """Detect flowlines that are not connected to the main network."""
        try:
            components = self._find_connected_components()
            if not components:
                return []

            # Find the largest component
            largest_component = max(components, key=len)

            # Find flowlines not in the largest component
            isolated_flowlines = []

            for i, flowline in enumerate(self.graph.aFlowline):
                if flowline is None:
                    continue

                # Check if flowline vertices are in largest component
                if hasattr(flowline, 'aVertex') and flowline.aVertex:
                    flowline_in_main = False

                    for vertex in flowline.aVertex:
                        if hasattr(vertex, 'lVertexIndex'):
                            if vertex.lVertexIndex in largest_component:
                                flowline_in_main = True
                                break

                    if not flowline_in_main:
                        isolated_flowlines.append(i)

            return isolated_flowlines

        except Exception as e:
            logger.error(f"Isolated flowline detection failed: {e}")
            return []

    def _detect_dangling_flowlines(self) -> List[int]:
        """Detect flowlines with only one connection (dead ends)."""
        try:
            dangling_flowlines = []

            for i, flowline in enumerate(self.graph.aFlowline):
                if flowline is None or not hasattr(flowline, 'aVertex') or not flowline.aVertex:
                    continue

                # Count connections at each end
                start_connections = 0
                end_connections = 0

                if len(flowline.aVertex) >= 2:
                    start_vertex = flowline.aVertex[0]
                    end_vertex = flowline.aVertex[-1]

                    if hasattr(start_vertex, 'lVertexIndex') and hasattr(end_vertex, 'lVertexIndex'):
                        start_id = start_vertex.lVertexIndex
                        end_id = end_vertex.lVertexIndex

                        # Count connections to other flowlines
                        for j, other_flowline in enumerate(self.graph.aFlowline):
                            if i == j or other_flowline is None:
                                continue

                            if hasattr(other_flowline, 'aVertex') and other_flowline.aVertex:
                                if len(other_flowline.aVertex) >= 2:
                                    other_start = other_flowline.aVertex[0]
                                    other_end = other_flowline.aVertex[-1]

                                    if (hasattr(other_start, 'lVertexIndex') and
                                        hasattr(other_end, 'lVertexIndex')):

                                        if other_start.lVertexIndex == start_id or other_end.lVertexIndex == start_id:
                                            start_connections += 1
                                        if other_start.lVertexIndex == end_id or other_end.lVertexIndex == end_id:
                                            end_connections += 1

                        # Flowline is dangling if it has only one connection
                        if start_connections == 0 or end_connections == 0:
                            dangling_flowlines.append(i)

            return dangling_flowlines

        except Exception as e:
            logger.error(f"Dangling flowline detection failed: {e}")
            return []

    def _detect_unusual_confluences(self) -> List[int]:
        """Detect confluences with unusual characteristics."""
        try:
            unusual_confluences = []

            for i, confluence in enumerate(self.graph.aConfluence):
                if confluence is None:
                    continue

                unusual = False

                # Check for too many connections (>4 is unusual)
                if hasattr(confluence, 'aFlowline') and confluence.aFlowline:
                    if len(confluence.aFlowline) > 4:
                        unusual = True

                    # Check for very acute angles
                    angles = self._calculate_junction_angles(confluence)
                    if any(angle < 15 for angle in angles):
                        unusual = True

                if unusual:
                    unusual_confluences.append(i)

            return unusual_confluences

        except Exception as e:
            logger.error(f"Unusual confluence detection failed: {e}")
            return []

    def _detect_flow_inconsistencies(self) -> List[int]:
        """Detect flowlines with inconsistent flow directions."""
        try:
            inconsistent_flowlines = []

            # This is a simplified check - in practice, you'd need elevation data
            # or other flow direction indicators

            for i, flowline in enumerate(self.graph.aFlowline):
                if flowline is None:
                    continue

                # Check if flowline has flow direction attribute
                if hasattr(flowline, 'iFlowDirection'):
                    # Check consistency with connected flowlines
                    # This would require more sophisticated analysis
                    pass

            return inconsistent_flowlines

        except Exception as e:
            logger.error(f"Flow inconsistency detection failed: {e}")
            return []