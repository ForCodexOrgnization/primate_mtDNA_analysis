#!/usr/bin/env bash
#SBATCH --job-name=primate_ds_prepare
#SBATCH --partition=day
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=results/qc_coverage_only/downsampling_pilot/logs/candidate_prepare_%j.out
#SBATCH --error=results/qc_coverage_only/downsampling_pilot/logs/candidate_prepare_%j.err
set -euo pipefail
module load SAMtools/1.21-GCC-13.3.0
ROOT=/vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis_coverage_only
cd "$ROOT"
python3 qc_analysis/scripts/candidate_downsampling.py prepare --threads "${SLURM_CPUS_PER_TASK:-4}" --mt-contig "${MT_CONTIG:-chrM}" --candidate-root "${CANDIDATE_ROOT:-$ROOT/results}"
