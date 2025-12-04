"""
Modular pyrivergraph package for river network graph analysis.

This package provides a comprehensive graph-based approach to river network analysis,
including braided channel removal, cycle detection, parallel path resolution, and
intelligent update strategies for dynamic network modifications.

The package is organized into focused modules for better maintainability:
- core: Main pyrivergraph class with simplified public API
- state_management: State tracking, caching, and performance management
- network_analysis: Topology analysis and pattern detection
- network_operations: Network modification and simplification operations
- stream_analysis: Stream ordering and hydrological calculations
- incremental_updates: Efficient incremental network updates
- batch_operations: Batch operation management and execution
- validation: Graph validation and performance analysis
- graph_builders: Graph construction and maintenance utilities
- utils: Shared utilities and helper functions
"""

# Import the main class for backward compatibility
from .core import pyrivergraph

# Import modular components for advanced usage
from .state_management import (
    GraphStateManager,
    GraphChangeType,
    GraphUpdateStrategy,
    GraphState,
    PerformanceMetrics
)

from .utils import RiverGraphUtils

from .batch_operations import BatchOperation

from .incremental_updates import IncrementalUpdates

from .validation import GraphValidator

from .network_analysis import NetworkAnalyzer

from .stream_analysis import StreamAnalyzer

from .network_operations import NetworkOperations

from .graph_builders import GraphBuilders

# Export the main class and key components
__all__ = [
    # Main class (backward compatibility)
    'pyrivergraph',

    # State management components
    'GraphStateManager',
    'GraphChangeType',
    'GraphUpdateStrategy',
    'GraphState',
    'PerformanceMetrics',

    # Modular components
    'RiverGraphUtils',
    'BatchOperation',
    'IncrementalUpdates',
    'GraphValidator',
    'NetworkAnalyzer',
    'StreamAnalyzer',
    'NetworkOperations',
    'GraphBuilders'
]

# Version information
__version__ = '2.0.0'
__author__ = 'PyFlowline Development Team'
__description__ = 'Modular river network graph analysis with hybrid update strategies'

# Backward compatibility aliases
# These ensure that existing code continues to work unchanged
pyRiverGraph = pyrivergraph  # Alternative capitalization
PyRiverGraph = pyrivergraph  # Alternative capitalization

# Add backward compatibility aliases to exports
__all__.extend(['pyRiverGraph', 'PyRiverGraph'])