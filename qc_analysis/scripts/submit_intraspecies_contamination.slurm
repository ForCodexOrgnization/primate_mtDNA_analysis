#!/usr/bin/env bash
#SBATCH --job-name=intraspecies_contamination
#SBATCH --output=logs/intraspecies_contamination/%x_%j.out
#SBATCH --error=logs/intraspecies_contamination/%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
set -euo pipefail
bash qc_analysis/scripts/run_intraspecies_contamination.sh "${1:-config/qc_preprocessing.yaml}"
