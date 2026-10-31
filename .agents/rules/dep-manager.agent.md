# Rule: Package Dependencies & Build Configuration

Guidelines for managing package configuration, dependencies, and build specifications in Python projects.

## Environment & Tool Execution
- **Dynamic Interpreter Resolution:** Resolve the active Python interpreter from `.vscode/settings.json` (`python.defaultInterpreterPath`) or `${workspaceFolder}/.venv/bin/python`.
- **Run Package Tools via Interpreter:**
  `<python_interpreter> -m pip install -e . --no-deps`
- **Never Run Bare `pip`:** Subshells may point to system Python instead of the project's environment.

## Packaging Standards
- Use modern `pyproject.toml` configuration. Do not introduce legacy `requirements.txt` files unless explicitly requested.
- Validate all newly added dependencies against existing version constraints before adding.
- **Fail-Fast Guard:** If editable install or build fails, STOP immediately. Output the exact build log to the user rather than retrying in an automated loop.
