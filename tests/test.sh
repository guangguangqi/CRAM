#!/bin/bash
set -euo pipefail

# Standard Harbor environment variables with defaults
TESTS_DIR="${TESTS_DIR:-/tests}"
LOG_DIR="${LOG_DIR:-/logs/verifier}"

# Create directory to store testing telemetry logs
mkdir -p "$LOG_DIR"

echo "[TESTING] Starting Pytest verification suite..."

# Execute pytest directly
if python3 -m pytest "$TESTS_DIR/test_outputs.py" -v \
  | tee "$LOG_DIR/pytest.log"
then
    echo "[PASSED] Pipeline successfully verified. Score: 1"
    echo "1" > "$LOG_DIR/reward.txt"
    EXIT_CODE=0
else
    echo "[FAILED] Pipeline failed validation criteria. Score: 0"
    echo "0" > "$LOG_DIR/reward.txt"
    EXIT_CODE=1
fi

exit $EXIT_CODE
