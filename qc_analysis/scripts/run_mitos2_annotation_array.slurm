#!/usr/bin/env bash
#SBATCH --job-name=mitos2_ref
#SBATCH --output=logs/mitos2_ref_%A_%a.out
#SBATCH --error=logs/mitos2_ref_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
# Set N to the number of data rows in mitos2_reference_tasks.tsv when submitting.
#SBATCH --array=1-N%20

set -euo pipefail

module load miniconda/24.11.3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mitos2

MITOS2_PYTHON="${CONDA_PREFIX}/bin/python"

echo "Python: ${MITOS2_PYTHON}"
"${MITOS2_PYTHON}" -c \
'import sys, Bio; print(sys.executable); print(Bio.__version__)'

"${MITOS2_PYTHON}" qc_analysis/scripts/run_mitos2_annotation.py \
  --config config/qc_preprocessing.yaml \
  --task-id "${SLURM_ARRAY_TASK_ID}"
