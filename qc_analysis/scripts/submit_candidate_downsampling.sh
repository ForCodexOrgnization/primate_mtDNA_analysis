#!/usr/bin/env bash
set -euo pipefail
ROOT=/vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis_coverage_only
cd "$ROOT"
BASE=results/qc_coverage_only/downsampling_pilot
mkdir -p "$BASE/logs"
MANIFEST="$BASE/candidate_downsampling_jobs.tsv"
MAX_CONCURRENT=${MAX_CONCURRENT:-4}
[[ "$MAX_CONCURRENT" =~ ^[1-9][0-9]*$ ]] || { echo "Invalid MAX_CONCURRENT" >&2; exit 2; }
case "${1:---help}" in
  --prepare)
    export CANDIDATE_ROOT=${CANDIDATE_ROOT:-$ROOT/results}
    export MT_CONTIG=${MT_CONTIG:-chrM}
    sbatch --export=ALL qc_analysis/scripts/prepare_candidate_downsampling.slurm
    ;;
  --pilot|--all)
    [[ -s "$MANIFEST" ]] || { echo "Missing $MANIFEST; run --prepare then review its log" >&2; exit 1; }
    N=$(( $(wc -l < "$MANIFEST") - 1 ))
    (( N>0 )) || { echo "Manifest has no tasks" >&2; exit 1; }
    if [[ "$1" == --pilot ]]; then
      sbatch --export=ALL --array=1 qc_analysis/scripts/run_candidate_downsampling_array.slurm
    else
      sbatch --export=ALL --array="1-${N}%${MAX_CONCURRENT}" qc_analysis/scripts/run_candidate_downsampling_array.slurm
    fi
    ;;
  *)
    echo "Usage: bash $0 --prepare|--pilot|--all" >&2
    exit 2
    ;;
esac
