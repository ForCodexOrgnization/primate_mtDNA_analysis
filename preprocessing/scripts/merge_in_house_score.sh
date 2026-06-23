#!/usr/bin/env bash
set -euo pipefail
INDIR=${1:-results/preprocessing/in_house_score}
OUT=${2:-results/preprocessing/in_house_score/merged_in_house_score.tsv}
SOURCE="${INDIR}/all_species.in_house_summary.with_numt_mask.tsv"
mkdir -p "$(dirname "$OUT")"
if [[ -s "$SOURCE" ]]; then
  cp "$SOURCE" "$OUT"
else
  shopt -s nullglob
  files=("$INDIR"/*.summary.tsv)
  if [[ "${#files[@]}" -eq 0 ]]; then
    echo "ERROR: no per-reference in-house score summaries found in ${INDIR}" >&2
    exit 1
  fi
  awk 'FNR==1 && NR!=1{next} {print}' "${files[@]}" > "$OUT"
fi
