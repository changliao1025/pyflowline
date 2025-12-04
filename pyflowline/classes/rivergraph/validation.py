"""
Validation module for pyrivergraph.

This module provides comprehensive validation and consistency checking
mechanisms for river network graphs, including topology validation,
geometric consistency checks, and data integrity verification.
"""

import time
import logging
from typing import List, Dict, Any, Optional, Set, Tuple, Union
from collections import defaultdict
from enum import Enum

from .state_management import GraphStateManager
from .utils import calculate_distance_2d, find_cycles_dfs

logger = logging.getLogger(__name__)


class ValidationLevel(Enum):
    """Validation levels for different checking intensities."""
    BASIC = "basic"
    STANDARD = "standard"
    COMPREHENSIVE = "comprehensive"
    STRICT = "strict"


class ValidationResult:
    """Container for validation results."""

    def __init__(self):
        self.is_valid = True
        self.errors = []
        self.warnings = []
        self.info = []
        self.execution_time = 0.0
        self.checks_performed = []

    def add_error(self, message: str, element_id: Optional[int] = None,
                  element_type: Optional[str] = None):
        """Add an error to the validation result."""
        self.is_valid = False
        error = {
            'message': message,
            'element_id': element_id,
            'element_type': element_type,
            'severity': 'error'
        }
        self.errors.append(error)

    def add_warning(self, message: str, element_id: Optional[int] = None,
                   element_type: Optional[str] = None):
        """Add a warning to the validation result."""
        warning = {
            'message': message,
            'element_id': element_id,
            'element_type': element_type,
            'severity': 'warning'
        }
        self.warnings.append(warning)

    def add_info(self, message: str, element_id: Optional[int] = None,
                element_type: Optional[str] = None):
        """Add an info message to the validation result."""
        info = {
            'message': message,
            'element_id': element_id,
            'element_type': element_type,
            'severity': 'info'
        }
        self.info.append(info)

    def get_summary(self) -> Dict[str, Any]:
        """Get a summary of validation results."""
        return {
            'is_valid': self.is_valid,
            'error_count': len(self.errors),
            'warning_count': len(self.warnings),
            'info_count': len(self.info),
            'execution_time': self.execution_time,
            'checks_performed': self.checks_performed
        }


class GraphValidator:
    """
    Comprehensive validator for river network graphs with multiple
    validation levels and detailed error reporting.
    """

    def __init__(self, graph_instance: Any, state_manager: Optional[GraphStateManager] = None):
        """
        Initialize the graph validator.

        Args:
            graph_instance: The pyrivergraph instance to validate
            state_manager: Optional state manager for tracking validation history
        """
        self.graph = graph_instance
        self.state_manager = state_manager
        self.tolerance = 1e-6  # Distance tolerance for geometric checks

        logger.debug("Initialized GraphValidator")

    def validate_graph(self, level: ValidationLevel = ValidationLevel.STANDARD,
                      check_topology: bool = True,
                      check_geometry: bool = True,
                      check_attributes: bool = True,
                      check_consistency: bool = True) -> ValidationResult:
        """
        Perform comprehensive graph validation.

        Args:
            level: Validation level (basic, standard, comprehensive, strict)
            check_topology: Whether to check topological consistency
            check_geometry: Whether to check geometric validity
            check_attributes: Whether to check attribute consistency
            check_consistency: Whether to check cross-element consistency

        Returns:
            ValidationResult with detailed findings
        """
        start_time = time.time()
        result = ValidationResult()

        try:
            logger.info(f"Starting graph validation at level: {level.value}")

            # Basic checks (always performed)
            if level in [ValidationLevel.BASIC, ValidationLevel.STANDARD,
                        ValidationLevel.COMPREHENSIVE, ValidationLevel.STRICT]:
                self._check_basic_structure(result)
                result.checks_performed.append("basic_structure")

            # Standard checks
            if level in [ValidationLevel.STANDARD, ValidationLevel.COMPREHENSIVE,
                        ValidationLevel.STRICT]:
                if check_topology:
                    self._check_topology(result)
                    result.checks_performed.append("topology")

                if check_geometry:
                    self._check_geometry(result)
                    result.checks_performed.append("geometry")

                if check_attributes:
                    self._check_attributes(result)
                    result.checks_performed.append("attributes")

            # Comprehensive checks
            if level in [ValidationLevel.COMPREHENSIVE, ValidationLevel.STRICT]:
                if check_consistency:
                    self._check_consistency(result)
                    result.checks_performed.append("consistency")

                self._check_network_connectivity(result)
                result.checks_performed.append("connectivity")

                self._check_flow_direction(result)
                result.checks_performed.append("flow_direction")

            # Strict checks
            if level == ValidationLevel.STRICT:
                self._check_strict_requirements(result)
                result.checks_performed.append("strict_requirements")

                self._check_performance_indicators(result)
                result.checks_performed.append("performance_indicators")

            result.execution_time = time.time() - start_time

            # Log summary
            summary = result.get_summary()
            logger.info(f"Validation completed: {summary}")

            return result

        except Exception as e:
            result.add_error(f"Validation failed with exception: {e}")
            result.execution_time = time.time() - start_time
            logger.error(f"Graph validation failed: {e}")
            return result

    def validate_flowline(self, flowline_id: int) -> ValidationResult:
        """
        Validate a specific flowline.

        Args:
            flowline_id: ID of the flowline to validate

        Returns:
            ValidationResult for the flowline
        """
        result = ValidationResult()

        try:
            if flowline_id < 0 or flowline_id >= len(self.graph.aFlowline):
                result.add_error(f"Invalid flowline ID: {flowline_id}", flowline_id, "flowline")
                return result

            flowline = self.graph.aFlowline[flowline_id]
            if flowline is None:
                result.add_error(f"Flowline {flowline_id} is None", flowline_id, "flowline")
                return result

            # Check flowline structure
            self._validate_flowline_structure(flowline, flowline_id, result)

            # Check flowline geometry
            self._validate_flowline_geometry(flowline, flowline_id, result)

            # Check flowline attributes
            self._validate_flowline_attributes(flowline, flowline_id, result)

            return result

        except Exception as e:
            result.add_error(f"Flowline validation failed: {e}", flowline_id, "flowline")
            return result

    def validate_vertex(self, vertex_id: int) -> ValidationResult:
        """
        Validate a specific vertex.

        Args:
            vertex_id: ID of the vertex to validate

        Returns:
            ValidationResult for the vertex
        """
        result = ValidationResult()

        try:
            if vertex_id < 0 or vertex_id >= len(self.graph.aVertex):
                result.add_error(f"Invalid vertex ID: {vertex_id}", vertex_id, "vertex")
                return result

            vertex = self.graph.aVertex[vertex_id]
            if vertex is None:
                result.add_error(f"Vertex {vertex_id} is None", vertex_id, "vertex")
                return result

            # Check vertex structure
            self._validate_vertex_structure(vertex, vertex_id, result)

            # Check vertex coordinates
            self._validate_vertex_coordinates(vertex, vertex_id, result)

            return result

        except Exception as e:
            result.add_error(f"Vertex validation failed: {e}", vertex_id, "vertex")
            return result

    def validate_confluence(self, confluence_id: int) -> ValidationResult:
        """
        Validate a specific confluence.

        Args:
            confluence_id: ID of the confluence to validate

        Returns:
            ValidationResult for the confluence
        """
        result = ValidationResult()

        try:
            if confluence_id < 0 or confluence_id >= len(self.graph.aConfluence):
                result.add_error(f"Invalid confluence ID: {confluence_id}", confluence_id, "confluence")
                return result

            confluence = self.graph.aConfluence[confluence_id]
            if confluence is None:
                result.add_error(f"Confluence {confluence_id} is None", confluence_id, "confluence")
                return result

            # Check confluence structure
            self._validate_confluence_structure(confluence, confluence_id, result)

            return result

        except Exception as e:
            result.add_error(f"Confluence validation failed: {e}", confluence_id, "confluence")
            return result

    def check_cycles(self) -> ValidationResult:
        """
        Check for cycles in the river network.

        Returns:
            ValidationResult with cycle information
        """
        result = ValidationResult()

        try:
            # Build adjacency list for cycle detection
            adjacency = defaultdict(list)

            for i, flowline in enumerate(self.graph.aFlowline):
                if flowline is None:
                    continue

                if hasattr(flowline, 'aVertex') and flowline.aVertex and len(flowline.aVertex) >= 2:
                    start_vertex = flowline.aVertex[0]
                    end_vertex = flowline.aVertex[-1]

                    if (hasattr(start_vertex, 'lVertexIndex') and
                        hasattr(end_vertex, 'lVertexIndex')):
                        adjacency[start_vertex.lVertexIndex].append(end_vertex.lVertexIndex)

            # Find cycles using DFS
            cycles = find_cycles_dfs(adjacency)

            if cycles:
                for cycle in cycles:
                    result.add_warning(f"Cycle detected: {' -> '.join(map(str, cycle))}")
            else:
                result.add_info("No cycles detected in the network")

            return result

        except Exception as e:
            result.add_error(f"Cycle detection failed: {e}")
            return result

    # Private validation methods

    def _check_basic_structure(self, result: ValidationResult):
        """Check basic graph structure."""
        try:
            # Check if required lists exist
            if not hasattr(self.graph, 'aFlowline'):
                result.add_error("Graph missing aFlowline list")
                return

            if not hasattr(self.graph, 'aVertex'):
                result.add_error("Graph missing aVertex list")
                return

            if not hasattr(self.graph, 'aConfluence'):
                result.add_error("Graph missing aConfluence list")
                return

            # Check list types
            if not isinstance(self.graph.aFlowline, list):
                result.add_error("aFlowline is not a list")

            if not isinstance(self.graph.aVertex, list):
                result.add_error("aVertex is not a list")

            if not isinstance(self.graph.aConfluence, list):
                result.add_error("aConfluence is not a list")

            # Check for empty graph
            if len(self.graph.aFlowline) == 0:
                result.add_warning("Graph has no flowlines")

            if len(self.graph.aVertex) == 0:
                result.add_warning("Graph has no vertices")

        except Exception as e:
            result.add_error(f"Basic structure check failed: {e}")

    def _check_topology(self, result: ValidationResult):
        """Check topological consistency."""
        try:
            # Check flowline-vertex relationships
            for i, flowline in enumerate(self.graph.aFlowline):
                if flowline is None:
                    continue

                if hasattr(flowline, 'aVertex') and flowline.aVertex:
                    for vertex in flowline.aVertex:
                        if hasattr(vertex, 'lVertexIndex'):
                            vertex_id = vertex.lVertexIndex
                            if vertex_id >= len(self.graph.aVertex):
                                result.add_error(
                                    f"Flowline {i} references non-existent vertex {vertex_id}",
                                    i, "flowline"
                                )
                            elif self.graph.aVertex[vertex_id] is None:
                                result.add_error(
                                    f"Flowline {i} references deleted vertex {vertex_id}",
                                    i, "flowline"
                                )

            # Check confluence-vertex relationships
            for i, confluence in enumerate(self.graph.aConfluence):
                if confluence is None:
                    continue

                if hasattr(confluence, 'lVertexIndex'):
                    vertex_id = confluence.lVertexIndex
                    if vertex_id >= len(self.graph.aVertex):
                        result.add_error(
                            f"Confluence {i} references non-existent vertex {vertex_id}",
                            i, "confluence"
                        )
                    elif self.graph.aVertex[vertex_id] is None:
                        result.add_error(
                            f"Confluence {i} references deleted vertex {vertex_id}",
                            i, "confluence"
                        )

        except Exception as e:
            result.add_error(f"Topology check failed: {e}")

    def _check_geometry(self, result: ValidationResult):
        """Check geometric validity."""
        try:
            # Check flowline geometry
            for i, flowline in enumerate(self.graph.aFlowline):
                if flowline is None:
                    continue

                self._validate_flowline_geometry(flowline, i, result)

            # Check vertex coordinates
            for i, vertex in enumerate(self.graph.aVertex):
                if vertex is None:
                    continue

                self._validate_vertex_coordinates(vertex, i, result)

        except Exception as e:
            result.add_error(f"Geometry check failed: {e}")

    def _check_attributes(self, result: ValidationResult):
        """Check attribute consistency."""
        try:
            # Check flowline attributes
            for i, flowline in enumerate(self.graph.aFlowline):
                if flowline is None:
                    continue

                self._validate_flowline_attributes(flowline, i, result)

        except Exception as e:
            result.add_error(f"Attribute check failed: {e}")

    def _check_consistency(self, result: ValidationResult):
        """Check cross-element consistency."""
        try:
            # Check vertex usage consistency
            vertex_usage = defaultdict(list)

            for i, flowline in enumerate(self.graph.aFlowline):
                if flowline is None:
                    continue

                if hasattr(flowline, 'aVertex') and flowline.aVertex:
                    for vertex in flowline.aVertex:
                        if hasattr(vertex, 'lVertexIndex'):
                            vertex_usage[vertex.lVertexIndex].append(i)

            # Check for unused vertices
            for i, vertex in enumerate(self.graph.aVertex):
                if vertex is not None and i not in vertex_usage:
                    result.add_warning(f"Vertex {i} is not used by any flowline", i, "vertex")

        except Exception as e:
            result.add_error(f"Consistency check failed: {e}")

    def _check_network_connectivity(self, result: ValidationResult):
        """Check network connectivity."""
        try:
            # Build connectivity graph
            connected_components = self._find_connected_components()

            if len(connected_components) > 1:
                result.add_warning(f"Network has {len(connected_components)} disconnected components")
                for i, component in enumerate(connected_components):
                    result.add_info(f"Component {i+1} has {len(component)} flowlines")
            else:
                result.add_info("Network is fully connected")

        except Exception as e:
            result.add_error(f"Connectivity check failed: {e}")

    def _check_flow_direction(self, result: ValidationResult):
        """Check flow direction consistency."""
        try:
            # Check for consistent flow direction
            inconsistent_flows = 0

            for i, flowline in enumerate(self.graph.aFlowline):
                if flowline is None:
                    continue

                if hasattr(flowline, 'aVertex') and flowline.aVertex and len(flowline.aVertex) >= 2:
                    # Check if flow direction is consistent with elevation
                    start_vertex = flowline.aVertex[0]
                    end_vertex = flowline.aVertex[-1]

                    if (hasattr(start_vertex, 'dElevation') and
                        hasattr(end_vertex, 'dElevation')):
                        if start_vertex.dElevation < end_vertex.dElevation:
                            result.add_warning(
                                f"Flowline {i} flows uphill (start: {start_vertex.dElevation}, "
                                f"end: {end_vertex.dElevation})", i, "flowline"
                            )
                            inconsistent_flows += 1

            if inconsistent_flows == 0:
                result.add_info("All flowlines have consistent flow direction")
            else:
                result.add_warning(f"Found {inconsistent_flows} flowlines with inconsistent flow direction")

        except Exception as e:
            result.add_error(f"Flow direction check failed: {e}")

    def _check_strict_requirements(self, result: ValidationResult):
        """Check strict requirements for high-quality networks."""
        try:
            # Check for minimum vertex count per flowline
            for i, flowline in enumerate(self.graph.aFlowline):
                if flowline is None:
                    continue

                if hasattr(flowline, 'aVertex') and flowline.aVertex:
                    if len(flowline.aVertex) < 2:
                        result.add_error(
                            f"Flowline {i} has fewer than 2 vertices", i, "flowline"
                        )

            # Check for duplicate vertices in flowlines
            for i, flowline in enumerate(self.graph.aFlowline):
                if flowline is None:
                    continue

                if hasattr(flowline, 'aVertex') and flowline.aVertex:
                    vertex_ids = []
                    for vertex in flowline.aVertex:
                        if hasattr(vertex, 'lVertexIndex'):
                            vertex_ids.append(vertex.lVertexIndex)

                    if len(vertex_ids) != len(set(vertex_ids)):
                        result.add_error(
                            f"Flowline {i} has duplicate vertices", i, "flowline"
                        )

        except Exception as e:
            result.add_error(f"Strict requirements check failed: {e}")

    def _check_performance_indicators(self, result: ValidationResult):
        """Check performance indicators and potential optimization opportunities."""
        try:
            # Check for very short flowlines
            short_flowlines = 0
            for i, flowline in enumerate(self.graph.aFlowline):
                if flowline is None:
                    continue

                if hasattr(flowline, 'dLength') and flowline.dLength < 1.0:  # Less than 1 meter
                    short_flowlines += 1

            if short_flowlines > 0:
                result.add_info(f"Found {short_flowlines} very short flowlines (< 1m)")

            # Check for very long flowlines
            long_flowlines = 0
            for i, flowline in enumerate(self.graph.aFlowline):
                if flowline is None:
                    continue

                if hasattr(flowline, 'dLength') and flowline.dLength > 10000.0:  # More than 10km
                    long_flowlines += 1

            if long_flowlines > 0:
                result.add_info(f"Found {long_flowlines} very long flowlines (> 10km)")

        except Exception as e:
            result.add_error(f"Performance indicators check failed: {e}")

    def _validate_flowline_structure(self, flowline: Any, flowline_id: int, result: ValidationResult):
        """Validate flowline structure."""
        try:
            # Check required attributes
            if not hasattr(flowline, 'lFlowlineIndex'):
                result.add_error(f"Flowline {flowline_id} missing lFlowlineIndex", flowline_id, "flowline")

            if not hasattr(flowline, 'aVertex'):
                result.add_error(f"Flowline {flowline_id} missing aVertex", flowline_id, "flowline")
            elif not isinstance(flowline.aVertex, list):
                result.add_error(f"Flowline {flowline_id} aVertex is not a list", flowline_id, "flowline")
            elif len(flowline.aVertex) < 2:
                result.add_error(f"Flowline {flowline_id} has fewer than 2 vertices", flowline_id, "flowline")

        except Exception as e:
            result.add_error(f"Flowline structure validation failed: {e}", flowline_id, "flowline")

    def _validate_flowline_geometry(self, flowline: Any, flowline_id: int, result: ValidationResult):
        """Validate flowline geometry."""
        try:
            if hasattr(flowline, 'aVertex') and flowline.aVertex:
                # Check vertex coordinates
                for j, vertex in enumerate(flowline.aVertex):
                    if not hasattr(vertex, 'dX') or not hasattr(vertex, 'dY'):
                        result.add_error(
                            f"Flowline {flowline_id} vertex {j} missing coordinates",
                            flowline_id, "flowline"
                        )

                # Check for zero-length segments
                for j in range(len(flowline.aVertex) - 1):
                    v1 = flowline.aVertex[j]
                    v2 = flowline.aVertex[j + 1]

                    if (hasattr(v1, 'dX') and hasattr(v1, 'dY') and
                        hasattr(v2, 'dX') and hasattr(v2, 'dY')):
                        distance = calculate_distance_2d(v1.dX, v1.dY, v2.dX, v2.dY)
                        if distance < self.tolerance:
                            result.add_warning(
                                f"Flowline {flowline_id} has zero-length segment at vertex {j}",
                                flowline_id, "flowline"
                            )

        except Exception as e:
            result.add_error(f"Flowline geometry validation failed: {e}", flowline_id, "flowline")

    def _validate_flowline_attributes(self, flowline: Any, flowline_id: int, result: ValidationResult):
        """Validate flowline attributes."""
        try:
            # Check length attribute
            if hasattr(flowline, 'dLength'):
                if flowline.dLength < 0:
                    result.add_error(f"Flowline {flowline_id} has negative length", flowline_id, "flowline")
                elif flowline.dLength == 0:
                    result.add_warning(f"Flowline {flowline_id} has zero length", flowline_id, "flowline")

            # Check stream order
            if hasattr(flowline, 'iStreamOrder'):
                if flowline.iStreamOrder < 1:
                    result.add_error(f"Flowline {flowline_id} has invalid stream order", flowline_id, "flowline")

        except Exception as e:
            result.add_error(f"Flowline attributes validation failed: {e}", flowline_id, "flowline")

    def _validate_vertex_structure(self, vertex: Any, vertex_id: int, result: ValidationResult):
        """Validate vertex structure."""
        try:
            # Check required attributes
            if not hasattr(vertex, 'lVertexIndex'):
                result.add_error(f"Vertex {vertex_id} missing lVertexIndex", vertex_id, "vertex")

            if not hasattr(vertex, 'dX'):
                result.add_error(f"Vertex {vertex_id} missing dX coordinate", vertex_id, "vertex")

            if not hasattr(vertex, 'dY'):
                result.add_error(f"Vertex {vertex_id} missing dY coordinate", vertex_id, "vertex")

        except Exception as e:
            result.add_error(f"Vertex structure validation failed: {e}", vertex_id, "vertex")

    def _validate_vertex_coordinates(self, vertex: Any, vertex_id: int, result: ValidationResult):
        """Validate vertex coordinates."""
        try:
            if hasattr(vertex, 'dX') and hasattr(vertex, 'dY'):
                # Check for valid coordinate values
                if not isinstance(vertex.dX, (int, float)):
                    result.add_error(f"Vertex {vertex_id} has invalid X coordinate type", vertex_id, "vertex")

                if not isinstance(vertex.dY, (int, float)):
                    result.add_error(f"Vertex {vertex_id} has invalid Y coordinate type", vertex_id, "vertex")

                # Check for reasonable coordinate ranges (assuming geographic coordinates)
                if abs(vertex.dX) > 180:
                    result.add_warning(f"Vertex {vertex_id} X coordinate outside typical range", vertex_id, "vertex")

                if abs(vertex.dY) > 90:
                    result.add_warning(f"Vertex {vertex_id} Y coordinate outside typical range", vertex_id, "vertex")

        except Exception as e:
            result.add_error(f"Vertex coordinates validation failed: {e}", vertex_id, "vertex")

    def _validate_confluence_structure(self, confluence: Any, confluence_id: int, result: ValidationResult):
        """Validate confluence structure."""
        try:
            # Check required attributes
            if not hasattr(confluence, 'lVertexIndex'):
                result.add_error(f"Confluence {confluence_id} missing lVertexIndex", confluence_id, "confluence")

            if hasattr(confluence, 'aFlowline'):
                if not isinstance(confluence.aFlowline, list):
                    result.add_error(f"Confluence {confluence_id} aFlowline is not a list", confluence_id, "confluence")
                elif len(confluence.aFlowline) < 2:
                    result.add_warning(f"Confluence {confluence_id} has fewer than 2 flowlines", confluence_id, "confluence")

        except Exception as e:
            result.add_error(f"Confluence structure validation failed: {e}", confluence_id, "confluence")

    def _find_connected_components(self) -> List[Set[int]]:
        """Find connected components in the network."""
        visited = set()
        components = []

        def dfs(flowline_id: int, component: Set[int]):
            if flowline_id in visited or flowline_id >= len(self.graph.aFlowline):
                return

            flowline = self.graph.aFlowline[flowline_id]
            if flowline is None:
                return

            visited.add(flowline_id)
            component.add(flowline_id)

            # Find connected flowlines through shared vertices
            if hasattr(flowline, 'aVertex') and flowline.aVertex:
                for vertex in flowline.aVertex:
                    if hasattr(vertex, 'lVertexIndex'):
                        vertex_id = vertex.lVertexIndex

                        # Find other flowlines using this vertex
                        for i, other_flowline in enumerate(self.graph.aFlowline):
                            if i != flowline_id and other_flowline is not None and i not in visited:
                                if hasattr(other_flowline, 'aVertex') and other_flowline.aVertex:
                                    for other_vertex in other_flowline.aVertex:
                                        if (hasattr(other_vertex, 'lVertexIndex') and
                                            other_vertex.lVertexIndex == vertex_id):
                                            dfs(i, component)

        for i, flowline in enumerate(self.graph.aFlowline):
            if flowline is not None and i not in visited:
                component = set()
                dfs(i, component)
                if component:
                    components.append(component)

        return components