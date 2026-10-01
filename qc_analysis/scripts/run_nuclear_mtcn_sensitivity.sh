#!/usr/bin/env bash
#SBATCH --job-name=qc_nuc_mtcn_sensitivity
#SBATCH --output=logs/qc_threshold_sensitivity/%x_%j.out
#SBATCH --error=logs/qc_threshold_sensitivity/%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -euo pipefail

mkdir -p logs/qc_threshold_sensitivity

PYTHON_BIN="${PYTHON:-python3}"

exec "$PYTHON_BIN" qc_analysis/scripts/run_nuclear_mtcn_sensitivity.py "$@"
