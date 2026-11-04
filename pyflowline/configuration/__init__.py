"""Configuration management for pyflowline.

Supports YAML (primary, comment-enabled, user-friendly) and JSON (programmatic,
backward-compatible) configuration files for case, basin, and jigsaw setups.

All configuration management, path discovery, workflow compilation, and model
loading functions are centralized in ``config_manager.py`` (matching the pyearthmesh
architecture).
"""

from __future__ import annotations

import sys

from pyflowline.configuration import config_manager
from pyflowline.configuration.config_manager import (
    DEFAULT_BASIN_CONFIG,
    DEFAULT_FLOWLINE_CONFIG,
    BasinConfigManager,
    FlowlineConfigManager,
    JigsawConfigManager,
    WorkflowSpec,
    compile_workflow,
    create_basin_template_configuration_file,
    create_jigsaw_template_configuration_file,
    create_model_configuration,
    create_pyflowline_template_configuration_file,
    create_workflow,
    dump_pyflowline_config,
    dump_workflow_config,
    join_project_path,
    load_configuration,
    load_pyflowline_config,
    load_workflow_config,
    path_manager,
    pyflowline_package_root,
    pyflowline_project_root,
    pyflowline_read_configuration_file,
    root_path_from_path_manager_location,
    root_path_from_pyflowline_package_root,
    root_path_from_setup_file,
)

# Backwards-compatibility module shims for legacy submodule imports
sys.modules.setdefault("pyflowline.configuration.config_manager", config_manager)
sys.modules.setdefault("pyflowline.configuration.path_manager", config_manager)
sys.modules.setdefault("pyflowline.configuration.read_configuration_file", config_manager)
sys.modules.setdefault("pyflowline.configuration.change_config_key_value", config_manager)
sys.modules.setdefault("pyflowline.configuration.change_json_key_value", config_manager)

__all__ = [
    "BasinConfigManager",
    "FlowlineConfigManager",
    "JigsawConfigManager",
    "WorkflowSpec",
    "compile_workflow",
    "create_workflow",
    "load_configuration",
    "load_pyflowline_config",
    "load_workflow_config",
    "dump_pyflowline_config",
    "dump_workflow_config",
    "pyflowline_read_configuration_file",
    "create_basin_template_configuration_file",
    "create_model_configuration",
    "create_pyflowline_template_configuration_file",
    "create_jigsaw_template_configuration_file",
    "DEFAULT_BASIN_CONFIG",
    "DEFAULT_FLOWLINE_CONFIG",
    "config_manager",
    "path_manager",
    "join_project_path",
    "pyflowline_project_root",
    "pyflowline_package_root",
    "root_path_from_setup_file",
    "root_path_from_pyflowline_package_root",
    "root_path_from_path_manager_location",
]
