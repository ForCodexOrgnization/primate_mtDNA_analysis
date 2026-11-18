#!/usr/bin/env bash
#SBATCH --job-name=round2_duplicate_audit
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=round2_duplicate_audit.%j.out
#SBATCH --error=round2_duplicate_audit.%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
python3 qc_analysis/scripts/audit_round2_duplicate_variants.py \
  --vcf-dir results/qc_coverage_only/collected_variant_calling_results/collected_vcf \
  --outdir results/qc_coverage_only_dp50/duplicate_variant_audit \
  --dp-min 50
