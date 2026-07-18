#!/usr/bin/env bash
#SBATCH --job-name=global_anchor_discovery
#SBATCH --output=logs/global_anchor_discovery/%x_%j.out
#SBATCH --error=logs/global_anchor_discovery/%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -euo pipefail
CONFIG="${CONFIG:-${1:-config/qc_preprocessing.yaml}}"
REPO_DIR="${REPO_DIR:-$(pwd)}"
PYTHON="${PYTHON:-python3}"
MAFFT_ENV_COMMAND="${MAFFT_ENV_COMMAND:-}"
cd "$REPO_DIR"
mkdir -p logs/global_anchor_discovery
if [[ -n "$MAFFT_ENV_COMMAND" ]]; then
  # shellcheck disable=SC1090
  source <(printf '%s\n' "$MAFFT_ENV_COMMAND")
fi
export MAFFT_N_THREADS="${SLURM_CPUS_PER_TASK:-${CPUS:-8}}"
"$PYTHON" qc_analysis/scripts/discover_global_liftover_anchor.py --config "$CONFIG"
