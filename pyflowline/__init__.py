from __future__ import annotations

import sys

from pyflowline import configuration
from pyflowline import configuration as config
from pyflowline.configuration import (
    create_model_configuration,
    create_workflow,
    load_configuration,
    load_pyflowline_config,
    load_workflow_config,
)

# Cross-project alias for consistency with pyearthmesh.config
sys.modules.setdefault("pyflowline.config", configuration)
sys.modules.setdefault("pyflowline.config.config_manager", configuration.config_manager)

__all__ = [
    "configuration",
    "config",
    "load_configuration",
    "load_workflow_config",
    "load_pyflowline_config",
    "create_model_configuration",
    "create_workflow",
]
