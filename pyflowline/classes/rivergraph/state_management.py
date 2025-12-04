"""
State management module for pyrivergraph.

This module provides comprehensive state tracking, caching, and performance management
for the river network graph. It includes change detection, update strategy decisions,
and performance monitoring capabilities.
"""

import time
import logging
from enum import Enum
from typing import Dict, Any, Optional, Set, List, Tuple
from dataclasses import dataclass, field
from collections import defaultdict

# Configure logging
logger = logging.getLogger(__name__)


class GraphChangeType(Enum):
    """Types of changes that can occur to the graph."""
    ADD_FLOWLINE = "add_flowline"
    REMOVE_FLOWLINE = "remove_flowline"
    MODIFY_FLOWLINE = "modify_flowline"
    ADD_VERTEX = "add_vertex"
    REMOVE_VERTEX = "remove_vertex"
    MODIFY_VERTEX = "modify_vertex"
    TOPOLOGY_CHANGE = "topology_change"
    BATCH_OPERATION = "batch_operation"
    FULL_REBUILD = "full_rebuild"


class GraphUpdateStrategy(Enum):
    """Strategies for updating the graph after changes."""
    INCREMENTAL = "incremental"
    PARTIAL_REBUILD = "partial_rebuild"
    FULL_REBUILD = "full_rebuild"
    LAZY_UPDATE = "lazy_update"


class GraphState(Enum):
    """Current state of the graph."""
    CLEAN = "clean"
    DIRTY = "dirty"
    UPDATING = "updating"
    INVALID = "invalid"


@dataclass
class PerformanceMetrics:
    """Performance metrics for graph operations."""
    operation_count: int = 0
    total_time: float = 0.0
    average_time: float = 0.0
    last_operation_time: float = 0.0
    cache_hits: int = 0
    cache_misses: int = 0

    def update(self, operation_time: float):
        """Update metrics with new operation time."""
        self.operation_count += 1
        self.total_time += operation_time
        self.average_time = self.total_time / self.operation_count
        self.last_operation_time = operation_time

    def cache_hit(self):
        """Record a cache hit."""
        self.cache_hits += 1

    def cache_miss(self):
        """Record a cache miss."""
        self.cache_misses += 1

    @property
    def cache_hit_rate(self) -> float:
        """Calculate cache hit rate."""
        total = self.cache_hits + self.cache_misses
        return self.cache_hits / total if total > 0 else 0.0


@dataclass
class GraphChange:
    """Represents a change to the graph."""
    change_type: GraphChangeType
    timestamp: float
    affected_elements: Set[int] = field(default_factory=set)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        if self.timestamp == 0:
            self.timestamp = time.time()


class GraphStateManager:
    """
    Manages the state of the river network graph, including change tracking,
    caching, and performance monitoring.
    """

    def __init__(self):
        """Initialize the state manager."""
        self.state = GraphState.CLEAN
        self.changes: List[GraphChange] = []
        self.cache: Dict[str, Any] = {}
        self.cache_timestamps: Dict[str, float] = {}
        self.dirty_properties: Set[str] = set()
        self.performance_metrics: Dict[str, PerformanceMetrics] = defaultdict(PerformanceMetrics)
        self.last_full_rebuild = time.time()
        self.incremental_update_threshold = 10  # Max changes before forcing rebuild
        self.cache_ttl = 300  # Cache time-to-live in seconds

        # Track which properties depend on which graph elements
        self.property_dependencies = {
            'stream_order': {'flowlines', 'topology'},
            'cycles': {'flowlines', 'topology'},
            'confluence_count': {'confluences', 'topology'},
            'vertex_count': {'vertices'},
            'flowline_count': {'flowlines'},
            'topology': {'flowlines', 'vertices', 'confluences'}
        }

    def record_change(self, change_type: GraphChangeType,
                     affected_elements: Optional[Set[int]] = None,
                     metadata: Optional[Dict[str, Any]] = None) -> None:
        """Record a change to the graph."""
        if affected_elements is None:
            affected_elements = set()
        if metadata is None:
            metadata = {}

        change = GraphChange(
            change_type=change_type,
            timestamp=time.time(),
            affected_elements=affected_elements,
            metadata=metadata
        )

        self.changes.append(change)
        self.state = GraphState.DIRTY

        # Invalidate dependent properties
        self._invalidate_dependent_properties(change_type)

        logger.debug(f"Recorded change: {change_type} affecting {len(affected_elements)} elements")

    def _invalidate_dependent_properties(self, change_type: GraphChangeType) -> None:
        """Invalidate cached properties that depend on the changed elements."""
        # Map change types to affected graph elements
        element_mapping = {
            GraphChangeType.ADD_FLOWLINE: {'flowlines', 'topology'},
            GraphChangeType.REMOVE_FLOWLINE: {'flowlines', 'topology'},
            GraphChangeType.MODIFY_FLOWLINE: {'flowlines', 'topology'},
            GraphChangeType.ADD_VERTEX: {'vertices', 'topology'},
            GraphChangeType.REMOVE_VERTEX: {'vertices', 'topology'},
            GraphChangeType.MODIFY_VERTEX: {'vertices'},
            GraphChangeType.TOPOLOGY_CHANGE: {'topology'},
            GraphChangeType.BATCH_OPERATION: {'flowlines', 'vertices', 'confluences', 'topology'},
            GraphChangeType.FULL_REBUILD: {'flowlines', 'vertices', 'confluences', 'topology'}
        }

        affected_elements = element_mapping.get(change_type, set())

        # Find properties that depend on these elements
        for prop, dependencies in self.property_dependencies.items():
            if dependencies.intersection(affected_elements):
                self.dirty_properties.add(prop)
                # Remove from cache
                if prop in self.cache:
                    del self.cache[prop]
                    del self.cache_timestamps[prop]

    def determine_update_strategy(self) -> GraphUpdateStrategy:
        """Determine the best update strategy based on recent changes."""
        if not self.changes:
            return GraphUpdateStrategy.INCREMENTAL

        recent_changes = [c for c in self.changes if time.time() - c.timestamp < 60]  # Last minute

        # Count different types of changes
        change_counts = defaultdict(int)
        total_affected_elements = set()

        for change in recent_changes:
            change_counts[change.change_type] += 1
            total_affected_elements.update(change.affected_elements)

        # Decision logic
        if len(recent_changes) > self.incremental_update_threshold:
            return GraphUpdateStrategy.FULL_REBUILD

        if GraphChangeType.FULL_REBUILD in change_counts:
            return GraphUpdateStrategy.FULL_REBUILD

        if GraphChangeType.BATCH_OPERATION in change_counts:
            return GraphUpdateStrategy.PARTIAL_REBUILD

        if GraphChangeType.TOPOLOGY_CHANGE in change_counts:
            return GraphUpdateStrategy.PARTIAL_REBUILD

        # Check if changes affect critical topology
        topology_changes = (
            change_counts[GraphChangeType.ADD_FLOWLINE] +
            change_counts[GraphChangeType.REMOVE_FLOWLINE] +
            change_counts[GraphChangeType.MODIFY_FLOWLINE]
        )

        if topology_changes > 5:
            return GraphUpdateStrategy.PARTIAL_REBUILD

        return GraphUpdateStrategy.INCREMENTAL

    def get_cached_property(self, property_name: str) -> Tuple[Any, bool]:
        """
        Get a cached property value.

        Returns:
            Tuple of (value, cache_hit) where cache_hit is True if found in cache
        """
        if property_name in self.cache:
            # Check if cache is still valid
            cache_age = time.time() - self.cache_timestamps[property_name]
            if cache_age < self.cache_ttl and property_name not in self.dirty_properties:
                self.performance_metrics[property_name].cache_hit()
                return self.cache[property_name], True
            else:
                # Cache expired or property is dirty
                del self.cache[property_name]
                del self.cache_timestamps[property_name]
                self.dirty_properties.discard(property_name)

        self.performance_metrics[property_name].cache_miss()
        return None, False

    def cache_property(self, property_name: str, value: Any) -> None:
        """Cache a property value."""
        self.cache[property_name] = value
        self.cache_timestamps[property_name] = time.time()
        self.dirty_properties.discard(property_name)

        logger.debug(f"Cached property: {property_name}")

    def clear_cache(self, property_names: Optional[List[str]] = None) -> None:
        """Clear cache for specific properties or all properties."""
        if property_names is None:
            self.cache.clear()
            self.cache_timestamps.clear()
            self.dirty_properties.clear()
            logger.debug("Cleared all cache")
        else:
            for prop in property_names:
                self.cache.pop(prop, None)
                self.cache_timestamps.pop(prop, None)
                self.dirty_properties.discard(prop)
            logger.debug(f"Cleared cache for: {property_names}")

    def start_operation(self, operation_name: str) -> float:
        """Start timing an operation."""
        return time.time()

    def end_operation(self, operation_name: str, start_time: float) -> float:
        """End timing an operation and update metrics."""
        operation_time = time.time() - start_time
        self.performance_metrics[operation_name].update(operation_time)
        return operation_time

    def mark_clean(self) -> None:
        """Mark the graph as clean (up-to-date)."""
        self.state = GraphState.CLEAN
        self.changes.clear()
        logger.debug("Graph marked as clean")

    def mark_updating(self) -> None:
        """Mark the graph as currently being updated."""
        self.state = GraphState.UPDATING

    def mark_invalid(self) -> None:
        """Mark the graph as invalid (needs rebuild)."""
        self.state = GraphState.INVALID
        self.clear_cache()

    def get_performance_summary(self) -> Dict[str, Any]:
        """Get a summary of performance metrics."""
        summary = {}
        for operation, metrics in self.performance_metrics.items():
            summary[operation] = {
                'operation_count': metrics.operation_count,
                'total_time': metrics.total_time,
                'average_time': metrics.average_time,
                'last_operation_time': metrics.last_operation_time,
                'cache_hit_rate': metrics.cache_hit_rate
            }

        summary['overall'] = {
            'total_changes': len(self.changes),
            'current_state': self.state.value,
            'cache_size': len(self.cache),
            'dirty_properties': len(self.dirty_properties),
            'time_since_last_rebuild': time.time() - self.last_full_rebuild
        }

        return summary

    def should_force_rebuild(self) -> bool:
        """Determine if a full rebuild should be forced."""
        # Force rebuild if too many changes have accumulated
        if len(self.changes) > self.incremental_update_threshold * 2:
            return True

        # Force rebuild if it's been too long since last rebuild
        time_since_rebuild = time.time() - self.last_full_rebuild
        if time_since_rebuild > 3600:  # 1 hour
            return True

        # Force rebuild if graph is in invalid state
        if self.state == GraphState.INVALID:
            return True

        return False

    def record_full_rebuild(self) -> None:
        """Record that a full rebuild has occurred."""
        self.last_full_rebuild = time.time()
        self.changes.clear()
        self.clear_cache()
        self.state = GraphState.CLEAN
        logger.info("Full rebuild completed and recorded")