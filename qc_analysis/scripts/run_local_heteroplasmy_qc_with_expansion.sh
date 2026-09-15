#!/usr/bin/env bash
set -euo pipefail

CONFIG="${1:-config/qc_preprocessing.yaml}"
PYTHON_BIN="${PYTHON:-python3}"

"$PYTHON_BIN" qc_analysis/scripts/run_local_heteroplasmy_qc.py \
  --config "$CONFIG"

"$PYTHON_BIN" qc_analysis/scripts/expand_numt_seed_variants.py \
  --config "$CONFIG"

"$PYTHON_BIN" qc_analysis/scripts/detect_indel_complex_regions.py \
  --config "$CONFIG"

printf '[local_heteroplasmy_qc_with_expansion] complete config=%s\n' "$CONFIG" >&2
