# Rule: Code Debugging & Troubleshooting

Guidelines and constraints for diagnosing, debugging, and resolving code defects and test failures.

## Python Environment & Execution
- **Dynamic Interpreter Resolution:** Always resolve the active Python interpreter before executing commands:
  1. Check `.vscode/settings.json` for `"python.defaultInterpreterPath"`.
  2. Fall back to `${workspaceFolder}/.venv/bin/python`, `${workspaceFolder}/venv/bin/python`, or `python3`.
- **Invoke via Interpreter:** Always run tools through the resolved interpreter:
  `<python_interpreter> -m pytest <path_to_test_file> [-k test_name] -v`
- **Never Run Bare `pytest` or `python`:** Non-interactive subshells may lack PATH activation, causing `command not found` or library import errors.
- **Targeted Test Execution Only:** Never run repository-wide test suites (`pytest tests/`) for localized fixes. Run only the specific test file or function relevant to the defect.

## Bug Fixing Constraints
- **Minimal Blast Radius:** Limit modifications strictly to the files directly responsible (maximum 3 files). Never reformat or refactor unrelated code.
- **No Backward-Compatibility Burden:** Do not introduce legacy shims, alias assignments, or wrapper functions. When replacing or modernizing an implementation, switch call sites over directly and cleanly.
- **Fail-Fast Guard (Max 2 Attempts):** If a fix does not pass tests after 2 attempts, STOP immediately. Summarize what was tried, report the error log, and request user guidance rather than looping.
- **Canonical Naming:** Always use canonical, standard snake_case function and module names.
