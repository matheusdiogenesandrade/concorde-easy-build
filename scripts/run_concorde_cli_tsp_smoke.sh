#!/usr/bin/env bash
# Smoke-test upstream Concorde CLI on symmetric TSPLIB instances (not ATSP).
# Default: Python runner mirrors discover/sort/10-bucket partition from
# external/concorde_wrapper/tests/test_concorde_symmetric_suite.cpp (111 runs).
set -euo pipefail

_script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_default_tsp="${_script_dir}/../data/tsp"
TSP_DIR="${TSP_DIR:-${_default_tsp}}"

if [[ -z "${CONCORDE_BIN:-}" ]]; then
  echo "run_concorde_cli_tsp_smoke.sh: set CONCORDE_BIN to the concorde_cli executable" >&2
  exit 2
fi
if [[ ! -x "$CONCORDE_BIN" ]]; then
  echo "run_concorde_cli_tsp_smoke.sh: not an executable: $CONCORDE_BIN" >&2
  exit 2
fi
if [[ ! -d "$TSP_DIR" ]]; then
  echo "run_concorde_cli_tsp_smoke.sh: TSP_DIR is not a directory: $TSP_DIR" >&2
  exit 2
fi

PYTHON="${PYTHON:-python3}"
_py_script="${_script_dir}/run_concorde_cli_tsp_buckets.py"
if [[ ! -f "$_py_script" ]]; then
  echo "run_concorde_cli_tsp_smoke.sh: missing $_py_script" >&2
  exit 2
fi
if ! command -v "$PYTHON" >/dev/null 2>&1; then
  echo "run_concorde_cli_tsp_smoke.sh: need python3 (set PYTHON to override)" >&2
  exit 2
fi

_extra=()
if [[ "${QUICK:-}" == "1" ]]; then
  _extra+=(--quick)
fi

# ALL=1 is documented as an alias for the same full run (no separate code path).

exec "$PYTHON" "$_py_script" --tsp-dir "$TSP_DIR" --concorde-bin "$CONCORDE_BIN" "${_extra[@]}"
