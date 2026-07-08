#!/usr/bin/env bash
#SBATCH --job-name=coordinate_liftover
#SBATCH --output=logs/coordinate_liftover/%x_%A_%a.out
#SBATCH --error=logs/coordinate_liftover/%x_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -euo pipefail

REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
CONFIG="${CONFIG:-config/qc_preprocessing.yaml}"
SAMPLE_REF="${SAMPLE_REF:-config/sample_ref_file.tsv}"
PYTHON="${PYTHON:-python}"
SCRIPT="qc_analysis/scripts/run_coordinate_liftover.py"
SAMPLE="${SAMPLE:-}"

cd "$REPO_DIR"
mkdir -p logs/coordinate_liftover

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" && -z "$SAMPLE" ]]; then
  if [[ ! -f "$SAMPLE_REF" ]]; then
    echo "sample_ref_file not found for Slurm array: $SAMPLE_REF" >&2
    exit 1
  fi
  SAMPLE=$(awk -v idx="$SLURM_ARRAY_TASK_ID" '
    BEGIN { FS = "\t" }
    /^[[:space:]]*#/ || NF < 2 { next }
    tolower($1) == "sample" { next }
    { n++; if (n == idx) { print $1; exit } }
  ' "$SAMPLE_REF")
  if [[ -z "$SAMPLE" ]]; then
    echo "No sample found for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID in $SAMPLE_REF" >&2
    exit 1
  fi
fi

cmd=("$PYTHON" "$SCRIPT" --config "$CONFIG")
if [[ -n "$SAMPLE" ]]; then
  cmd+=(--sample "$SAMPLE")
fi

printf 'Running:' >&2
printf ' %q' "${cmd[@]}" >&2
printf '\n' >&2
exec "${cmd[@]}"
