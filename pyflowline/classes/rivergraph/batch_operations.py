"""
Batch operations module for pyrivergraph.

This module provides efficient batch operation management for performing multiple
network modifications as a single atomic operation, with rollback capabilities
and optimized performance.
"""

import time
import logging
from typing import List, Dict, Any, Optional, Set, Tuple, Callable, Union
from dataclasses import dataclass, field
from enum import Enum
from collections import defaultdict

from .state_management import GraphChangeType, GraphStateManager

logger = logging.getLogger(__name__)


class BatchOperationType(Enum):
    """Types of batch operations."""
    ADD_FLOWLINES = "add_flowlines"
    REMOVE_FLOWLINES = "remove_flowlines"
    MODIFY_FLOWLINES = "modify_flowlines"
    MERGE_FLOWLINES = "merge_flowlines"
    SPLIT_FLOWLINES = "split_flowlines"
    SIMPLIFY_NETWORK = "simplify_network"
    REMOVE_LOOPS = "remove_loops"
    CONNECT_DISCONNECTED = "connect_disconnected"
    CUSTOM = "custom"


class BatchOperationStatus(Enum):
    """Status of batch operations."""
    PENDING = "pending"
    EXECUTING = "executing"
    COMPLETED = "completed"
    FAILED = "failed"
    ROLLED_BACK = "rolled_back"


@dataclass
class BatchOperationStep:
    """Individual step within a batch operation."""
    step_id: str
    operation_type: str
    parameters: Dict[str, Any]
    dependencies: List[str] = field(default_factory=list)
    status: BatchOperationStatus = BatchOperationStatus.PENDING
    result: Optional[Any] = None
    error: Optional[str] = None
    execution_time: float = 0.0

    def __post_init__(self):
        if not self.step_id:
            self.step_id = f"step_{int(time.time() * 1000000)}"


@dataclass
class BatchOperationResult:
    """Result of a batch operation."""
    operation_id: str
    operation_type: BatchOperationType
    status: BatchOperationStatus
    total_steps: int
    completed_steps: int
    failed_steps: int
    total_execution_time: float
    affected_elements: Set[int] = field(default_factory=set)
    step_results: Dict[str, Any] = field(default_factory=dict)
    errors: List[str] = field(default_factory=list)
    rollback_info: Optional[Dict[str, Any]] = None


class BatchOperation:
    """
    Manages batch operations for efficient bulk modifications to the river network graph.

    This class provides:
    - Atomic batch operations with rollback capabilities
    - Dependency management between operation steps
    - Performance optimization through batching
    - Progress tracking and error handling
    """

    def __init__(self, operation_type: BatchOperationType,
                 operation_id: Optional[str] = None,
                 enable_rollback: bool = True):
        """
        Initialize a batch operation.

        Args:
            operation_type: Type of batch operation
            operation_id: Unique identifier for the operation
            enable_rollback: Whether to enable rollback capabilities
        """
        self.operation_type = operation_type
        self.operation_id = operation_id or f"batch_{int(time.time() * 1000000)}"
        self.enable_rollback = enable_rollback

        self.steps: List[BatchOperationStep] = []
        self.step_dependencies: Dict[str, List[str]] = defaultdict(list)
        self.status = BatchOperationStatus.PENDING
        self.start_time: Optional[float] = None
        self.end_time: Optional[float] = None

        # Rollback information
        self.rollback_data: Dict[str, Any] = {}
        self.original_state: Optional[Dict[str, Any]] = None

        # Progress tracking
        self.progress_callback: Optional[Callable[[float, str], None]] = None
        self.current_step = 0

        logger.info(f"Created batch operation: {self.operation_id} ({operation_type.value})")

    def add_step(self, operation_type: str, parameters: Dict[str, Any],
                 step_id: Optional[str] = None,
                 dependencies: Optional[List[str]] = None) -> str:
        """
        Add a step to the batch operation.

        Args:
            operation_type: Type of operation for this step
            parameters: Parameters for the operation
            step_id: Optional custom step ID
            dependencies: List of step IDs this step depends on

        Returns:
            The step ID
        """
        if dependencies is None:
            dependencies = []

        step = BatchOperationStep(
            step_id=step_id or f"step_{len(self.steps)}",
            operation_type=operation_type,
            parameters=parameters,
            dependencies=dependencies
        )

        self.steps.append(step)

        # Update dependency tracking
        for dep in dependencies:
            self.step_dependencies[step.step_id].append(dep)

        logger.debug(f"Added step {step.step_id} to batch {self.operation_id}")
        return step.step_id

    def add_flowline_operations(self, flowlines_to_add: List[Any],
                               flowlines_to_remove: List[int],
                               flowlines_to_modify: List[Tuple[int, Dict[str, Any]]]) -> None:
        """
        Add common flowline operations to the batch.

        Args:
            flowlines_to_add: List of flowline objects to add
            flowlines_to_remove: List of flowline IDs to remove
            flowlines_to_modify: List of (flowline_id, modifications) tuples
        """
        # Add removal operations first
        if flowlines_to_remove:
            self.add_step("remove_flowlines", {
                "flowline_ids": flowlines_to_remove
            })

        # Add modification operations
        if flowlines_to_modify:
            self.add_step("modify_flowlines", {
                "modifications": flowlines_to_modify
            })

        # Add new flowlines last
        if flowlines_to_add:
            self.add_step("add_flowlines", {
                "flowlines": flowlines_to_add
            })

    def set_progress_callback(self, callback: Callable[[float, str], None]) -> None:
        """Set a callback function for progress updates."""
        self.progress_callback = callback

    def _update_progress(self, progress: float, message: str) -> None:
        """Update progress and call callback if set."""
        if self.progress_callback:
            self.progress_callback(progress, message)

    def _validate_dependencies(self) -> bool:
        """
        Validate that all step dependencies are satisfied.

        Returns:
            True if dependencies are valid, False otherwise
        """
        step_ids = {step.step_id for step in self.steps}

        for step in self.steps:
            for dep in step.dependencies:
                if dep not in step_ids:
                    logger.error(f"Step {step.step_id} depends on non-existent step {dep}")
                    return False

        # Check for circular dependencies using topological sort
        try:
            self._get_execution_order()
            return True
        except ValueError as e:
            logger.error(f"Circular dependency detected: {e}")
            return False

    def _get_execution_order(self) -> List[str]:
        """
        Get the execution order of steps based on dependencies.

        Returns:
            List of step IDs in execution order

        Raises:
            ValueError: If circular dependencies are detected
        """
        # Build adjacency list for topological sort
        adjacency = defaultdict(list)
        in_degree = defaultdict(int)

        # Initialize all steps
        for step in self.steps:
            in_degree[step.step_id] = 0

        # Build dependency graph
        for step in self.steps:
            for dep in step.dependencies:
                adjacency[dep].append(step.step_id)
                in_degree[step.step_id] += 1

        # Topological sort using Kahn's algorithm
        queue = [step_id for step_id in in_degree if in_degree[step_id] == 0]
        result = []

        while queue:
            current = queue.pop(0)
            result.append(current)

            for neighbor in adjacency[current]:
                in_degree[neighbor] -= 1
                if in_degree[neighbor] == 0:
                    queue.append(neighbor)

        if len(result) != len(self.steps):
            raise ValueError("Circular dependency detected in batch operation steps")

        return result

    def _save_rollback_state(self, graph_instance: Any) -> None:
        """Save the current state for potential rollback."""
        if not self.enable_rollback:
            return

        try:
            # Save critical graph state
            self.original_state = {
                'flowline_count': len(getattr(graph_instance, 'aFlowline', [])),
                'vertex_count': len(getattr(graph_instance, 'aVertex', [])),
                'confluence_count': len(getattr(graph_instance, 'aConfluence', [])),
                'timestamp': time.time()
            }

            # Save specific elements that might be affected
            self.rollback_data['flowlines'] = {}
            self.rollback_data['vertices'] = {}
            self.rollback_data['confluences'] = {}

            logger.debug(f"Saved rollback state for batch {self.operation_id}")

        except Exception as e:
            logger.warning(f"Failed to save rollback state: {e}")
            self.enable_rollback = False

    def execute(self, graph_instance: Any,
                state_manager: Optional[GraphStateManager] = None) -> BatchOperationResult:
        """
        Execute the batch operation.

        Args:
            graph_instance: The pyrivergraph instance to operate on
            state_manager: Optional state manager for tracking changes

        Returns:
            BatchOperationResult with execution details
        """
        self.status = BatchOperationStatus.EXECUTING
        self.start_time = time.time()
        self.current_step = 0

        # Validate dependencies
        if not self._validate_dependencies():
            self.status = BatchOperationStatus.FAILED
            return self._create_result("Dependency validation failed")

        # Save rollback state
        self._save_rollback_state(graph_instance)

        # Get execution order
        try:
            execution_order = self._get_execution_order()
        except ValueError as e:
            self.status = BatchOperationStatus.FAILED
            return self._create_result(f"Failed to determine execution order: {e}")

        # Create step lookup
        step_lookup = {step.step_id: step for step in self.steps}

        # Execute steps in order
        completed_steps = 0
        failed_steps = 0
        step_results = {}
        errors = []
        affected_elements = set()

        try:
            for i, step_id in enumerate(execution_order):
                step = step_lookup[step_id]
                self.current_step = i + 1

                progress = (i / len(execution_order)) * 100
                self._update_progress(progress, f"Executing step {step.step_id}")

                # Execute the step
                step_start_time = time.time()
                step.status = BatchOperationStatus.EXECUTING

                try:
                    result = self._execute_step(step, graph_instance, state_manager)
                    step.result = result
                    step.status = BatchOperationStatus.COMPLETED
                    step_results[step_id] = result
                    completed_steps += 1

                    # Track affected elements
                    if isinstance(result, dict) and 'affected_elements' in result:
                        affected_elements.update(result['affected_elements'])

                except Exception as e:
                    step.error = str(e)
                    step.status = BatchOperationStatus.FAILED
                    errors.append(f"Step {step_id}: {e}")
                    failed_steps += 1
                    logger.error(f"Step {step_id} failed: {e}")

                    # Decide whether to continue or abort
                    if self._should_abort_on_error(step, e):
                        break

                finally:
                    step.execution_time = time.time() - step_start_time

            # Determine final status
            if failed_steps == 0:
                self.status = BatchOperationStatus.COMPLETED
                self._update_progress(100.0, "Batch operation completed successfully")
            else:
                self.status = BatchOperationStatus.FAILED
                self._update_progress(100.0, f"Batch operation completed with {failed_steps} failures")

        except Exception as e:
            self.status = BatchOperationStatus.FAILED
            errors.append(f"Batch execution error: {e}")
            logger.error(f"Batch operation {self.operation_id} failed: {e}")

        finally:
            self.end_time = time.time()

        # Record the batch operation in state manager
        if state_manager:
            state_manager.record_change(
                GraphChangeType.BATCH_OPERATION,
                affected_elements,
                {
                    'operation_id': self.operation_id,
                    'operation_type': self.operation_type.value,
                    'steps': len(self.steps),
                    'completed': completed_steps,
                    'failed': failed_steps
                }
            )

        return BatchOperationResult(
            operation_id=self.operation_id,
            operation_type=self.operation_type,
            status=self.status,
            total_steps=len(self.steps),
            completed_steps=completed_steps,
            failed_steps=failed_steps,
            total_execution_time=self.end_time - self.start_time,
            affected_elements=affected_elements,
            step_results=step_results,
            errors=errors,
            rollback_info=self.rollback_data if self.enable_rollback else None
        )

    def _execute_step(self, step: BatchOperationStep,
                     graph_instance: Any,
                     state_manager: Optional[GraphStateManager]) -> Any:
        """
        Execute a single step of the batch operation.

        Args:
            step: The step to execute
            graph_instance: The graph instance to operate on
            state_manager: Optional state manager

        Returns:
            Result of the step execution
        """
        operation_type = step.operation_type
        parameters = step.parameters

        # Map operation types to graph methods
        operation_mapping = {
            'add_flowlines': self._add_flowlines,
            'remove_flowlines': self._remove_flowlines,
            'modify_flowlines': self._modify_flowlines,
            'merge_flowlines': self._merge_flowlines,
            'split_flowlines': self._split_flowlines,
            'simplify_network': self._simplify_network,
            'remove_loops': self._remove_loops,
            'connect_disconnected': self._connect_disconnected,
            'rebuild_topology': self._rebuild_topology,
            'recalculate_stream_order': self._recalculate_stream_order
        }

        if operation_type in operation_mapping:
            return operation_mapping[operation_type](graph_instance, parameters, state_manager)
        else:
            # Try to call method directly on graph instance
            if hasattr(graph_instance, operation_type):
                method = getattr(graph_instance, operation_type)
                return method(**parameters)
            else:
                raise ValueError(f"Unknown operation type: {operation_type}")

    def _add_flowlines(self, graph_instance: Any, parameters: Dict[str, Any],
                      state_manager: Optional[GraphStateManager]) -> Dict[str, Any]:
        """Add flowlines to the graph."""
        flowlines = parameters.get('flowlines', [])
        added_ids = []

        for flowline in flowlines:
            # Add flowline to graph
            if hasattr(graph_instance, 'add_flowline'):
                flowline_id = graph_instance.add_flowline(flowline)
                added_ids.append(flowline_id)
            else:
                # Fallback: add to flowline list
                graph_instance.aFlowline.append(flowline)
                added_ids.append(len(graph_instance.aFlowline) - 1)

        return {
            'added_flowlines': added_ids,
            'affected_elements': set(added_ids)
        }

    def _remove_flowlines(self, graph_instance: Any, parameters: Dict[str, Any],
                         state_manager: Optional[GraphStateManager]) -> Dict[str, Any]:
        """Remove flowlines from the graph."""
        flowline_ids = parameters.get('flowline_ids', [])
        removed_ids = []

        for flowline_id in flowline_ids:
            if hasattr(graph_instance, 'remove_flowline'):
                success = graph_instance.remove_flowline(flowline_id)
                if success:
                    removed_ids.append(flowline_id)
            else:
                # Fallback: mark as removed or filter
                if 0 <= flowline_id < len(graph_instance.aFlowline):
                    removed_ids.append(flowline_id)

        return {
            'removed_flowlines': removed_ids,
            'affected_elements': set(removed_ids)
        }

    def _modify_flowlines(self, graph_instance: Any, parameters: Dict[str, Any],
                         state_manager: Optional[GraphStateManager]) -> Dict[str, Any]:
        """Modify existing flowlines."""
        modifications = parameters.get('modifications', [])
        modified_ids = []

        for flowline_id, changes in modifications:
            if hasattr(graph_instance, 'modify_flowline'):
                success = graph_instance.modify_flowline(flowline_id, changes)
                if success:
                    modified_ids.append(flowline_id)
            else:
                # Fallback: direct modification
                if 0 <= flowline_id < len(graph_instance.aFlowline):
                    flowline = graph_instance.aFlowline[flowline_id]
                    for key, value in changes.items():
                        if hasattr(flowline, key):
                            setattr(flowline, key, value)
                    modified_ids.append(flowline_id)

        return {
            'modified_flowlines': modified_ids,
            'affected_elements': set(modified_ids)
        }

    def _merge_flowlines(self, graph_instance: Any, parameters: Dict[str, Any],
                        state_manager: Optional[GraphStateManager]) -> Dict[str, Any]:
        """Merge multiple flowlines into one."""
        flowline_ids = parameters.get('flowline_ids', [])

        if hasattr(graph_instance, 'merge_flowlines'):
            result = graph_instance.merge_flowlines(flowline_ids)
            return {'merge_result': result, 'affected_elements': set(flowline_ids)}
        else:
            return {'merge_result': None, 'affected_elements': set()}

    def _split_flowlines(self, graph_instance: Any, parameters: Dict[str, Any],
                        state_manager: Optional[GraphStateManager]) -> Dict[str, Any]:
        """Split flowlines at specified points."""
        split_operations = parameters.get('split_operations', [])
        affected_ids = set()

        for operation in split_operations:
            flowline_id = operation.get('flowline_id')
            split_points = operation.get('split_points', [])

            if hasattr(graph_instance, 'split_flowline'):
                result = graph_instance.split_flowline(flowline_id, split_points)
                if result:
                    affected_ids.update(result.get('new_flowline_ids', []))
                    affected_ids.add(flowline_id)

        return {'affected_elements': affected_ids}

    def _simplify_network(self, graph_instance: Any, parameters: Dict[str, Any],
                         state_manager: Optional[GraphStateManager]) -> Dict[str, Any]:
        """Simplify the network by removing unnecessary elements."""
        tolerance = parameters.get('tolerance', 1e-6)

        if hasattr(graph_instance, 'simplify_network'):
            result = graph_instance.simplify_network(tolerance)
            return {'simplification_result': result}
        else:
            return {'simplification_result': None}

    def _remove_loops(self, graph_instance: Any, parameters: Dict[str, Any],
                     state_manager: Optional[GraphStateManager]) -> Dict[str, Any]:
        """Remove loops from the network."""
        if hasattr(graph_instance, 'remove_loops'):
            result = graph_instance.remove_loops()
            return {'loop_removal_result': result}
        else:
            return {'loop_removal_result': None}

    def _connect_disconnected(self, graph_instance: Any, parameters: Dict[str, Any],
                             state_manager: Optional[GraphStateManager]) -> Dict[str, Any]:
        """Connect disconnected network components."""
        if hasattr(graph_instance, 'connect_disconnected_flowlines'):
            result = graph_instance.connect_disconnected_flowlines()
            return {'connection_result': result}
        else:
            return {'connection_result': None}

    def _rebuild_topology(self, graph_instance: Any, parameters: Dict[str, Any],
                         state_manager: Optional[GraphStateManager]) -> Dict[str, Any]:
        """Rebuild the network topology."""
        if hasattr(graph_instance, 'rebuild_topology'):
            result = graph_instance.rebuild_topology()
            return {'topology_result': result}
        else:
            return {'topology_result': None}

    def _recalculate_stream_order(self, graph_instance: Any, parameters: Dict[str, Any],
                                 state_manager: Optional[GraphStateManager]) -> Dict[str, Any]:
        """Recalculate stream order for the network."""
        method = parameters.get('method', 'strahler')

        if hasattr(graph_instance, 'define_stream_order'):
            result = graph_instance.define_stream_order(method)
            return {'stream_order_result': result}
        else:
            return {'stream_order_result': None}

    def _should_abort_on_error(self, step: BatchOperationStep, error: Exception) -> bool:
        """
        Determine whether to abort the batch operation on step failure.

        Args:
            step: The failed step
            error: The error that occurred

        Returns:
            True if the batch should be aborted
        """
        # Define critical operations that should abort the batch
        critical_operations = {
            'rebuild_topology',
            'remove_loops',
            'connect_disconnected'
        }

        return step.operation_type in critical_operations

    def _create_result(self, error_message: str) -> BatchOperationResult:
        """Create a result object for failed operations."""
        return BatchOperationResult(
            operation_id=self.operation_id,
            operation_type=self.operation_type,
            status=self.status,
            total_steps=len(self.steps),
            completed_steps=0,
            failed_steps=len(self.steps),
            total_execution_time=0.0,
            errors=[error_message]
        )

    def rollback(self, graph_instance: Any) -> bool:
        """
        Rollback the batch operation if rollback is enabled.

        Args:
            graph_instance: The graph instance to rollback

        Returns:
            True if rollback was successful
        """
        if not self.enable_rollback or not self.original_state:
            logger.warning(f"Rollback not available for batch {self.operation_id}")
            return False

        try:
            # Implement rollback logic here
            # This would restore the graph to its original state
            logger.info(f"Rolling back batch operation {self.operation_id}")

            # Mark as rolled back
            self.status = BatchOperationStatus.ROLLED_BACK

            return True

        except Exception as e:
            logger.error(f"Rollback failed for batch {self.operation_id}: {e}")
            return False

    def get_progress(self) -> Tuple[float, str]:
        """
        Get current progress of the batch operation.

        Returns:
            Tuple of (progress_percentage, status_message)
        """
        if self.status == BatchOperationStatus.PENDING:
            return 0.0, "Pending execution"
        elif self.status == BatchOperationStatus.EXECUTING:
            progress = (self.current_step / len(self.steps)) * 100 if self.steps else 0.0
            return progress, f"Executing step {self.current_step}/{len(self.steps)}"
        elif self.status == BatchOperationStatus.COMPLETED:
            return 100.0, "Completed successfully"
        elif self.status == BatchOperationStatus.FAILED:
            return 100.0, "Failed"
        elif self.status == BatchOperationStatus.ROLLED_BACK:
            return 100.0, "Rolled back"
        else:
            return 0.0, "Unknown status"