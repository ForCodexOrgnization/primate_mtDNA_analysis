#!/usr/bin/env bash
#SBATCH --job-name=qc_mtcov_percent_sensitivity
#SBATCH --output=logs/qc_threshold_sensitivity/%x_%j.out
#SBATCH --error=logs/qc_threshold_sensitivity/%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G

set -euo pipefail

mkdir -p logs/qc_threshold_sensitivity

PYTHON_BIN="${PYTHON:-python3}"

exec "$PYTHON_BIN" qc_analysis/scripts/run_mtcov_percent_sensitivity.py "$@"
