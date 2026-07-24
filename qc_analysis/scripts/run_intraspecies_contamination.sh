#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

PYTHON="${PYTHON:-python3}"
CONFIG="${1:-$REPO_ROOT/config/qc_preprocessing.yaml}"

[[ -s "$CONFIG" ]] || {
    echo "ERROR: missing config: $CONFIG" >&2
    exit 2
}

exec "$PYTHON" \
    "$REPO_ROOT/qc_analysis/scripts/run_intraspecies_contamination.py" \
    --config "$CONFIG"
