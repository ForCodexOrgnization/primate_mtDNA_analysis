#!/bin/bash
#SBATCH --job-name=index_fasta
#SBATCH --output=log/preprocessing/index_fasta_%A_%a.out
#SBATCH --error=log/preprocessing/index_fasta_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=06:00:00
set -euo pipefail

FASTA_DIR="${1:-}"
if [[ -z "$FASTA_DIR" || ! -d "$FASTA_DIR" ]]; then
  echo "Usage: bash preprocessing/scripts/index_all_fasta.sh <fasta_dir>" >&2
  exit 2
fi
SAMTOOLS_COMMAND="${SAMTOOLS_COMMAND:-samtools}"
BWA_COMMAND="${BWA_COMMAND:-bwa}"
GATK_COMMAND="${GATK_COMMAND:-gatk}"
ARRAY_CONCURRENCY="${ARRAY_CONCURRENCY:-50}"
mkdir -p log/preprocessing

mapfile -t FASTAS < <(find "$FASTA_DIR" -maxdepth 1 -type f \( -name '*.fa' -o -name '*.fasta' \) | sort)
N=${#FASTAS[@]}
if (( N == 0 )); then
  echo "No .fa or .fasta files found in ${FASTA_DIR}" >&2
  exit 0
fi

if [[ "${SUBMIT_ARRAY:-0}" == "1" && -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  sbatch --array="1-${N}%${ARRAY_CONCURRENCY}" --export=ALL "$0" "$FASTA_DIR"
  exit 0
fi

idx=0
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
    echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} outside 1-${N}" >&2
    exit 2
  fi
  idx=$((SLURM_ARRAY_TASK_ID - 1))
  FASTAS=("${FASTAS[$idx]}")
fi

index_one() {
  local fa="$1" dict="${fa%.*}.dict"
  local expected=("${fa}.fai" "${fa}.amb" "${fa}.ann" "${fa}.bwt" "${fa}.pac" "${fa}.sa" "$dict")
  local present=0
  for f in "${expected[@]}"; do [[ -s "$f" ]] && present=$((present+1)); done
  if (( present == ${#expected[@]} )); then
    echo "Skipping ${fa}; all indexes exist"
    return 0
  fi
  if (( present > 0 )); then
    echo "Removing incomplete indexes for ${fa}"
    rm -f "${expected[@]}"
  fi
  "$SAMTOOLS_COMMAND" faidx "$fa"
  "$BWA_COMMAND" index "$fa"
  "$GATK_COMMAND" CreateSequenceDictionary -R "$fa" -O "$dict"
}

for fa in "${FASTAS[@]}"; do
  index_one "$fa"
done
