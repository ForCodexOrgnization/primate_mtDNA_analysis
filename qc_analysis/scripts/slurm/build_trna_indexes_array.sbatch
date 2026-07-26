#!/usr/bin/env bash
#SBATCH --job-name=trna_indexes
#SBATCH --output=results/qc/trna_match/logs/%A_%a.out
#SBATCH --error=results/qc/trna_match/logs/%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=04:00:00
set -euo pipefail
MANIFEST=${1:?usage: $0 MANIFEST [CONFIG]}; CONFIG=${2:-config/qc_preprocessing.yaml}
module load miniconda/24.11.3
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate trnascan_env
command -v tRNAscan-SE >/dev/null || { echo 'tRNAscan-SE unavailable after activating trnascan_env' >&2; exit 127; }
line=$(awk -F '\t' -v n="$((SLURM_ARRAY_TASK_ID+1))" 'NR==n{print;exit}' "$MANIFEST")
IFS=$'\t' read -r task_id reference_key fasta output_index <<<"$line"
[[ -n "$reference_key" && -s "$fasta" ]] || { echo "Invalid task $SLURM_ARRAY_TASK_ID" >&2;exit 2; }
python qc_analysis/scripts/build_trna_position_index.py --reference-key "$reference_key" --fasta "$fasta" --run-trnascan --trnascan-bin tRNAscan-SE --trnascan-mode mito_mammal --threads "${SLURM_CPUS_PER_TASK:-4}" --output "$output_index"
