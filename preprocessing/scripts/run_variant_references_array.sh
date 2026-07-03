#!/usr/bin/env bash
set -euo pipefail

REF_INPUTS=${REF_INPUTS:-references/manifests/in_house_score_reference_inputs.tsv}
SCORE=${SCORE:-results/preprocessing/in_house_score/merged_in_house_score.tsv}
OUT_ROOT=${OUT_ROOT:-references/variant_calling}
MASK_REF_TYPES=${MASK_REF_TYPES:-#C-likely_comp,#C-Ambiguous}
PYTHON_COMMAND=${PYTHON_COMMAND:-python3}
SAMTOOLS_COMMAND=${SAMTOOLS_COMMAND:-samtools}
BWA_COMMAND=${BWA_COMMAND:-bwa}
GATK_COMMAND=${GATK_COMMAND:-gatk}
MAX_CONCURRENT=${VARIANT_REFERENCE_THREADS:-${MAX_CONCURRENT:-1}}
VARIANT_REFERENCE_SLURM_TIME=${VARIANT_REFERENCE_SLURM_TIME:-24:00:00}
VARIANT_REFERENCE_SLURM_CPUS=${VARIANT_REFERENCE_SLURM_CPUS:-4}
VARIANT_REFERENCE_SLURM_MEM=${VARIANT_REFERENCE_SLURM_MEM:-16G}
VARIANT_REFERENCE_FORCE=${VARIANT_REFERENCE_FORCE:-${FORCE:-0}}
VARIANT_REFERENCE_SKIP_EXISTING=${VARIANT_REFERENCE_SKIP_EXISTING:-${SKIP_EXISTING:-0}}
VARIANT_REFERENCE_LOG_DIR=${VARIANT_REFERENCE_LOG_DIR:-logs/preprocessing}
BUILD_SCRIPT=${BUILD_SCRIPT:-preprocessing/scripts/build_variant_calling_references.sh}
SHARD_DIR="${OUT_ROOT}/manifest_shards"

is_truthy() {
  case "${1,,}" in
    1|true|yes|y) return 0 ;;
    *) return 1 ;;
  esac
}

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
export VARIANT_REFERENCE_FORCE
export VARIANT_REFERENCE_SKIP_EXISTING
export VARIANT_REFERENCE_THREADS

run_one() {
  echo "[variant_refs] MASK_REF_TYPES=${MASK_REF_TYPES}" >&2
  echo "[variant_refs] VARIANT_REFERENCE_THREADS=${VARIANT_REFERENCE_THREADS:-NA} MAX_CONCURRENT=${MAX_CONCURRENT}" >&2
  local idx="$1"
  local force_args=()
  if is_truthy "$VARIANT_REFERENCE_FORCE"; then
    force_args+=(--force)
  fi
  local skip_existing_args=()
  if is_truthy "$VARIANT_REFERENCE_SKIP_EXISTING"; then
    skip_existing_args+=(--skip-existing)
  fi
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
      "${force_args[@]}" \
      "${skip_existing_args[@]}" \
      --job-index "$idx" \
      --manifest-output "$shard"
}

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  run_one "$SLURM_ARRAY_TASK_ID"
  exit 0
fi

if [[ "${RUN_LOCAL:-0}" != "1" ]] && command -v sbatch >/dev/null 2>&1; then
  # MASK_REF_TYPES may contain commas; exporting it separately avoids Slurm
  # splitting the value inside the comma-delimited --export list.
  export MASK_REF_TYPES
  echo "[variant_refs] submitting array 1-${N}%${MAX_CONCURRENT}" >&2
  echo "[variant_refs] time=${VARIANT_REFERENCE_SLURM_TIME} cpus=${VARIANT_REFERENCE_SLURM_CPUS} mem=${VARIANT_REFERENCE_SLURM_MEM}" >&2
  jid=$(sbatch --parsable \
    --job-name="variant_refs" \
    --output="${VARIANT_REFERENCE_LOG_DIR}/variant_refs_%A_%a.out" \
    --error="${VARIANT_REFERENCE_LOG_DIR}/variant_refs_%A_%a.err" \
    --time="$VARIANT_REFERENCE_SLURM_TIME" \
    --cpus-per-task="$VARIANT_REFERENCE_SLURM_CPUS" \
    --mem="$VARIANT_REFERENCE_SLURM_MEM" \
    --array="1-${N}%${MAX_CONCURRENT}" \
    --export=ALL,REF_INPUTS="$REF_INPUTS",SCORE="$SCORE",OUT_ROOT="$OUT_ROOT",PYTHON_COMMAND="$PYTHON_COMMAND",SAMTOOLS_COMMAND="$SAMTOOLS_COMMAND",BWA_COMMAND="$BWA_COMMAND",GATK_COMMAND="$GATK_COMMAND",VARIANT_REFERENCE_LOG_DIR="$VARIANT_REFERENCE_LOG_DIR",BUILD_SCRIPT="$BUILD_SCRIPT",VARIANT_REFERENCE_SLURM_TIME="$VARIANT_REFERENCE_SLURM_TIME",VARIANT_REFERENCE_SLURM_CPUS="$VARIANT_REFERENCE_SLURM_CPUS",VARIANT_REFERENCE_SLURM_MEM="$VARIANT_REFERENCE_SLURM_MEM",VARIANT_REFERENCE_THREADS="$VARIANT_REFERENCE_THREADS",VARIANT_REFERENCE_FORCE="$VARIANT_REFERENCE_FORCE",VARIANT_REFERENCE_SKIP_EXISTING="$VARIANT_REFERENCE_SKIP_EXISTING" \
    "$0")
  echo "Submitted variant-reference array job ${jid} for ${N} reference packages (max concurrent ${MAX_CONCURRENT})." >&2
  merge_jid=$(sbatch --parsable \
    --job-name="variant_refs_merge" \
    --output="${VARIANT_REFERENCE_LOG_DIR}/variant_refs_merge_%j.out" \
    --error="${VARIANT_REFERENCE_LOG_DIR}/variant_refs_merge_%j.err" \
    --time="$VARIANT_REFERENCE_SLURM_TIME" \
    --cpus-per-task="1" \
    --mem="8G" \
    --dependency="afterok:${jid}" \
    --export=ALL,REF_INPUTS="$REF_INPUTS",SCORE="$SCORE",OUT_ROOT="$OUT_ROOT",PYTHON_COMMAND="$PYTHON_COMMAND",VARIANT_REFERENCE_LOG_DIR="$VARIANT_REFERENCE_LOG_DIR",MERGE_ONLY=1,VARIANT_REFERENCE_FORCE="$VARIANT_REFERENCE_FORCE",VARIANT_REFERENCE_SKIP_EXISTING="$VARIANT_REFERENCE_SKIP_EXISTING" \
    "$0")
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
