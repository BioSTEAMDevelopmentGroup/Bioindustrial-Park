# Development Guidelines

## 1. Directory Structure Rules

### Execution Logs
All execution logs, verification scripts, and temporary testing files must be located within the specific version package directory, specifically under `execution_log/`.

**Rule:** `execution_log/` MUST NOT exist in the root `biorefineries/` or project root. It MUST be located at `prefers/v<version>/execution_log/`.

**Example:**
- **Correct:** `biorefineries/prefers/v1/execution_log/test_my_unit.py`
- **Incorrect:** `execution_log/test_my_unit.py`
- **Incorrect:** `biorefineries/execution_log/test_my_unit.py`

This ensures that all testing and execution artifacts are versioned alongside the code they test.
