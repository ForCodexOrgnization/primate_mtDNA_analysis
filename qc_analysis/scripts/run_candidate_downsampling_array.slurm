#!/usr/bin/env bash
#SBATCH --job-name=primate_candidate_ds
#SBATCH --partition=day
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=results/qc_coverage_only/downsampling_pilot/logs/candidate_ds_%A_%a.out
#SBATCH --error=results/qc_coverage_only/downsampling_pilot/logs/candidate_ds_%A_%a.err
set -euo pipefail
module load SAMtools/1.21-GCC-13.3.0
ROOT=/vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis_coverage_only
cd "$ROOT"
python3 qc_analysis/scripts/candidate_downsampling.py run --task "${SLURM_ARRAY_TASK_ID:?}" --threads "${SLURM_CPUS_PER_TASK:-4}"
