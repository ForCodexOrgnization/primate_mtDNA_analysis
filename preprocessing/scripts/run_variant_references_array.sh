#!/usr/bin/env bash
set -euo pipefail

REF_INPUTS=${REF_INPUTS:-references/manifests/in_house_score_reference_inputs.tsv}
SCORE=${SCORE:-results/preprocessing/in_house_score/merged_in_house_score.tsv}
OUT_ROOT=${OUT_ROOT:-references/variant_calling}
MASK_REF_TYPES=${MASK_REF_TYPES:-#C-likely_comp,#C-Ambiguous,#A}
PYTHON_COMMAND=${PYTHON_COMMAND:-python3}
SAMTOOLS_COMMAND=${SAMTOOLS_COMMAND:-samtools}
BWA_COMMAND=${BWA_COMMAND:-bwa}
GATK_COMMAND=${GATK_COMMAND:-gatk}
MAX_CONCURRENT=${VARIANT_REFERENCE_THREADS:-${MAX_CONCURRENT:-1}}
VARIANT_REFERENCE_LOG_DIR=${VARIANT_REFERENCE_LOG_DIR:-logs/preprocessing}
BUILD_SCRIPT=${BUILD_SCRIPT:-preprocessing/scripts/build_variant_calling_references.sh}
SHARD_DIR="${OUT_ROOT}/manifest_shards"

if [[ ! -s "$REF_INPUTS" ]]; then
  echo "ERROR: missing or empty reference inputs: ${REF_INPUTS}" >&2
  exit 1
fi
if [[ ! -s "$SCORE" ]]; then
  echo "ERROR: missing or empty merged in-house score table: ${SCORE}" >&2
  exit 1
fi

N=$("$PYTHON_COMMAND" - "$REF_INPUTS" <<'PY'
import csv, sys
seen = set()
with open(sys.argv[1], newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        key = (row.get("target_species", ""), row.get("wg_fasta_path", ""), row.get("chrM_fasta_path", ""))
        if all(key):
            seen.add(key)
print(len(seen))
PY
)

if [[ "$N" -lt 1 ]]; then
  echo "ERROR: no usable reference rows in ${REF_INPUTS}" >&2
  exit 1
fi

merge_shards() {
  mkdir -p "$OUT_ROOT"
  local manifest="${OUT_ROOT}/variant_calling_reference_manifest.tsv"
  local first="${SHARD_DIR}/000001.tsv"
  if [[ ! -s "$first" ]]; then
    echo "ERROR: missing first variant-reference shard: ${first}" >&2
    exit 1
  fi
  head -n 1 "$first" > "$manifest"
  for i in $(seq 1 "$N"); do
    local shard
    shard=$(printf "%s/%06d.tsv" "$SHARD_DIR" "$i")
    if [[ ! -s "$shard" ]]; then
      echo "ERROR: missing variant-reference shard: ${shard}" >&2
      exit 1
    fi
    tail -n +2 "$shard" >> "$manifest"
  done
  "$PYTHON_COMMAND" preprocessing/scripts/build_variant_calling_references.py --out-root "$OUT_ROOT" --validate-only
}

if [[ "${MERGE_ONLY:-0}" == "1" ]]; then
  merge_shards
  exit 0
fi

mkdir -p "$OUT_ROOT" "$SHARD_DIR" "$VARIANT_REFERENCE_LOG_DIR"

run_one() {
  local idx="$1"
  local shard
  shard=$(printf "%s/%06d.tsv" "$SHARD_DIR" "$idx")
  PYTHON_COMMAND="$PYTHON_COMMAND" \
  SAMTOOLS_COMMAND="$SAMTOOLS_COMMAND" \
  BWA_COMMAND="$BWA_COMMAND" \
  GATK_COMMAND="$GATK_COMMAND" \
    bash "$BUILD_SCRIPT" \
      --ref-inputs "$REF_INPUTS" \
      --score "$SCORE" \
      --out-root "$OUT_ROOT" \
      --mask-ref-types "$MASK_REF_TYPES" \
      --job-index "$idx" \
      --manifest-output "$shard"
}

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  run_one "$SLURM_ARRAY_TASK_ID"
  exit 0
fi

if [[ "${RUN_LOCAL:-0}" != "1" ]] && command -v sbatch >/dev/null 2>&1; then
  jid=$(sbatch --parsable --array="1-${N}%${MAX_CONCURRENT}" --export=ALL,REF_INPUTS="$REF_INPUTS",SCORE="$SCORE",OUT_ROOT="$OUT_ROOT",MASK_REF_TYPES="$MASK_REF_TYPES",PYTHON_COMMAND="$PYTHON_COMMAND",SAMTOOLS_COMMAND="$SAMTOOLS_COMMAND",BWA_COMMAND="$BWA_COMMAND",GATK_COMMAND="$GATK_COMMAND",VARIANT_REFERENCE_LOG_DIR="$VARIANT_REFERENCE_LOG_DIR",BUILD_SCRIPT="$BUILD_SCRIPT" "$0")
  echo "Submitted variant-reference array job ${jid} for ${N} reference packages (max concurrent ${MAX_CONCURRENT})." >&2
  merge_jid=$(sbatch --parsable --dependency="afterok:${jid}" --export=ALL,REF_INPUTS="$REF_INPUTS",SCORE="$SCORE",OUT_ROOT="$OUT_ROOT",PYTHON_COMMAND="$PYTHON_COMMAND",VARIANT_REFERENCE_LOG_DIR="$VARIANT_REFERENCE_LOG_DIR",MERGE_ONLY=1 "$0")
  echo "Submitted variant-reference manifest merge job ${merge_jid} after array ${jid}." >&2
else
  echo "sbatch unavailable or RUN_LOCAL=1; building ${N} reference packages locally with ${MAX_CONCURRENT} workers." >&2
  running=0
  for i in $(seq 1 "$N"); do
    SLURM_ARRAY_TASK_ID="$i" RUN_LOCAL=1 "$0" &
    running=$((running + 1))
    if [[ "$running" -ge "$MAX_CONCURRENT" ]]; then
      wait -n
      running=$((running - 1))
    fi
  done
  wait
  merge_shards
fi
