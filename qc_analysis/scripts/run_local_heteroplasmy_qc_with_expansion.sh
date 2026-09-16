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

# Remove dense indel-overlap residual HETs before downstream residual analyses.
"$PYTHON_BIN" qc_analysis/scripts/apply_dense_indel_overlap_filter.py \
  --config "$CONFIG"

# Keep sensitivity diagnostics on the post-dense-indel residual universe.
"$PYTHON_BIN" qc_analysis/scripts/analyze_residual_artifact_sensitivity.py \
  --config "$CONFIG"

# Apply validated indel-seed expansion and recurrent-block production rules.
"$PYTHON_BIN" qc_analysis/scripts/apply_residual_artifact_rules.py \
  --config "$CONFIG"

printf '[local_heteroplasmy_qc_with_expansion] complete config=%s\n' "$CONFIG" >&2
