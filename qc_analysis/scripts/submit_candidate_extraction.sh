#!/usr/bin/env bash
set -euo pipefail
ROOT=/vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis_coverage_only
cd "$ROOT"
BASE=results/qc_coverage_only/downsampling_pilot
PILOT=$BASE/pilot_samples.tsv
[[ -s "$PILOT" ]] || { echo "Missing $PILOT" >&2; exit 1; }
N=$(( $(wc -l < "$PILOT") - 1 ))
(( N > 0 )) || { echo "No pilot samples" >&2; exit 1; }
mkdir -p "$BASE/logs"
MAX_CONCURRENT=${MAX_CONCURRENT:-2}
[[ "$MAX_CONCURRENT" =~ ^[1-9][0-9]*$ ]] || exit 2
export CANDIDATE_ROOT=${CANDIDATE_ROOT:-$ROOT/$BASE/reextracted_candidate}
export CANDIDATE_METADATA_ROOT=${CANDIDATE_METADATA_ROOT:-/home/lt692/ycga_work/primate_downsampling_pilot/candidate_input}
export WGS_INPUT_ROOT=${WGS_INPUT_ROOT:-/home/lt692/ycga_work/primate_downsampling_pilot/wgs_input}
export CRAM_REFERENCE_FASTA=${CRAM_REFERENCE_FASTA:-}
if [[ ${1:-} == --pilot ]]; then
  sbatch --export=ALL --array=1 qc_analysis/scripts/extract_candidate_from_cram.slurm
elif [[ ${1:-} == --all ]]; then
  sbatch --export=ALL --array="1-$N%$MAX_CONCURRENT" qc_analysis/scripts/extract_candidate_from_cram.slurm
else
  echo "Usage: bash $0 --pilot|--all" >&2; exit 2
fi
