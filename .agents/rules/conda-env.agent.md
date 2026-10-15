# Rule: Python Environment & Interpreter Resolution

Guidelines for resolving and verifying Python environments across workspaces.

## Saved Workspace Configuration (Default)
- **Primary Resolution:** Check `.vscode/settings.json` for `"python.defaultInterpreterPath"`.
- **Secondary Resolution:** Check for local virtual environments (`${workspaceFolder}/.venv` or `${workspaceFolder}/venv`).
- **Verification (Max 1 Step):**
  `<saved_interpreter> --version`
  If functional, use this interpreter path for all Python tool execution and proceed directly.

## Constraints
- Never execute verbose discovery commands such as `conda config --show` or `conda info` (they waste tens of thousands of tokens).
- Never create, modify, or destroy Conda or virtual environments automatically without explicit user instructions.
- Never run dynamic environment updates (`conda env update`) during task execution.
