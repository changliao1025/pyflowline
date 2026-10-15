# Rule: Unit Test Authoring & Repair

Guidelines for writing and repairing pytest unit tests across Python workspaces.

## Execution & Environment
- **Dynamic Interpreter Resolution:** Resolve the active interpreter from `.vscode/settings.json` (`python.defaultInterpreterPath`) or `${workspaceFolder}/.venv/bin/python`.
- **Run Tests via Interpreter:**
  `<python_interpreter> -m pytest <path_to_test_file> -v`
- **Never Run Bare `pytest`:** Subshells lack conda/venv activation and will fail with `command not found`.
- **Targeted Test Execution Only:** Run ONLY the specific test file being authored or repaired. Never run the full test suite.

## Test Design Guidelines
- **Targeted Scope:** Write or repair tests only for the specific module or function requested by the user. Do not initiate open-ended test coverage audits across the entire repository.
- **Placement:** Place test files inside `tests/` following the pattern `tests/test_<module_name>.py`.
- **Coverage Focus:** Validate happy paths and critical edge cases (e.g., empty sequences, boundary values, invalid inputs).
- **Fail-Fast Guard (Max 2 Attempts):** If a new test fails after 2 repair attempts, STOP immediately. Report the failing assertion, expected vs actual values, and ask the user for clarification.
- **Canonical Naming:** Test canonical functions and classes directly; never test or create legacy compatibility aliases.