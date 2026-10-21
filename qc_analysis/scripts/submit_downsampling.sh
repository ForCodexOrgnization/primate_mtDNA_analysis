#!/usr/bin/env bash
set -euo pipefail
ROOT=/vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis_coverage_only
cd "$ROOT"
BASE=results/qc_coverage_only/downsampling_pilot
MANIFEST="$BASE/downsampling_jobs.tsv"
[[ -s "$MANIFEST" ]] || { echo "Run prepare_downsampling_manifest.sh first" >&2; exit 1; }
N=$(( $(wc -l < "$MANIFEST") - 1 ))
(( N > 0 )) || { echo "No tasks" >&2; exit 1; }
mkdir -p "$BASE/logs"
export MT_CONTIG=${MT_CONTIG:-}
export CRAM_REFERENCE_FASTA=${CRAM_REFERENCE_FASTA:-}
export WGS_INPUT_ROOT=${WGS_INPUT_ROOT:-/home/lt692/ycga_work/primate_downsampling_pilot/wgs_input}
MAX_CONCURRENT=${MAX_CONCURRENT:-4}
[[ "$MAX_CONCURRENT" =~ ^[1-9][0-9]*$ ]] || { echo "Invalid MAX_CONCURRENT=$MAX_CONCURRENT" >&2; exit 2; }
echo "WGS source: $WGS_INPUT_ROOT"
echo "Array concurrency: $MAX_CONCURRENT"
# Optional: --pilot submits task 1 only; otherwise all with max concurrency 6.
if [[ ${1:-} == --pilot ]]; then
  sbatch --export=ALL --array=1 qc_analysis/scripts/run_downsampling_array.slurm
elif [[ $# -eq 0 ]]; then
  sbatch --export=ALL --array="1-${N}%${MAX_CONCURRENT}" qc_analysis/scripts/run_downsampling_array.slurm
else
  echo "Usage: bash $0 [--pilot]" >&2; exit 2
fi