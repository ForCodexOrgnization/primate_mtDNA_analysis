#!/usr/bin/env bash
set -euo pipefail
"${PYTHON_COMMAND:-python3}" preprocessing/scripts/build_variant_calling_references.py \
  --ref-inputs references/manifests/in_house_score_reference_inputs.tsv \
  --score results/preprocessing/in_house_score/merged_in_house_score.tsv \
  --out-root references/variant_calling \
  --shift 8000 "$@"
