#!/usr/bin/env bash
set -euo pipefail
INDIR=${1:-results/preprocessing/in_house_score}
OUT=${2:-results/preprocessing/in_house_score/merged_in_house_score.tsv}
mkdir -p "$(dirname "$OUT")"
awk 'FNR==1 && NR!=1{next} {print}' "$INDIR"/*.tsv > "$OUT"
