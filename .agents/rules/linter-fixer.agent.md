# Rule: Code Quality, Typing & Formatting

Guidelines for code style, type annotations, and linting across Python workspaces.

## Scope & Target Bounds
- **Targeted Scope Only:** Inspect and format ONLY the specific files requested by the user or touched during the current task. Never run repository-wide lint or format sweeps.
- **Dynamic Package Structure:** Identify the package root by reading `pyproject.toml` or finding the directory containing `__init__.py`. Tests live in `tests/`.

## Quality Guidelines
- Include type hints for all public functions, arguments, and return types.
- Write clean Google-style docstrings for public classes and functions.
- Preserve all existing comments and docstrings unrelated to your changes.
- Never introduce backward-compatible shims or duplicate alias assignments.

## Tool Execution & Fail-Fast Guard
- **Execute via Resolved Python:**
  `<python_interpreter> -m ruff check <target_file>`
  `<python_interpreter> -m flake8 <target_file>`
- **Missing Tool Guard:** If a quality tool (`ruff`, `mypy`, `flake8`) is not installed in the environment, report that it is missing and complete your turn immediately. Do NOT attempt automated installation or retry in a loop.
- **Single Pass Only:** Never iterate continuously trying to achieve "zero warnings". Apply clean fixes in one pass, report any remaining warnings in your final summary, and stop.
