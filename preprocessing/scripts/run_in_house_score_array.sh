#!/usr/bin/env bash
set -euo pipefail

REF_INPUTS=${REF_INPUTS:-references/manifests/in_house_score_reference_inputs.tsv}
OUTDIR=${OUTDIR:-results/preprocessing/in_house_score}
MERGED_IN_HOUSE_SCORE=${MERGED_IN_HOUSE_SCORE:-${OUTDIR}/merged_in_house_score.tsv}
IN_HOUSE_SCORE_SCRIPT=${IN_HOUSE_SCORE_SCRIPT:-preprocessing/scripts/in_house_score_with_minimal_numt_mask.sh}
PYTHON_COMMAND=${PYTHON_COMMAND:-python3}
MAX_CONCURRENT=${MAX_CONCURRENT:-50}
IN_HOUSE_SCORE_LOG_DIR=${IN_HOUSE_SCORE_LOG_DIR:-log/preprocessing}

if [[ ! -s "$REF_INPUTS" ]]; then
  echo "ERROR: missing or empty in-house score reference inputs: ${REF_INPUTS}" >&2
  exit 1
fi

N=$(
  "$PYTHON_COMMAND" - "$REF_INPUTS" <<'PY'
import csv, sys
seen = set()
with open(sys.argv[1], newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        species = (row.get("target_species") or "").strip()
        wg = (row.get("wg_fasta_path") or "").strip()
        chrm = (row.get("chrM_fasta_path") or "").strip()
        if species and wg and chrm:
            seen.add((species, wg, chrm))
print(len(seen))
PY
)

if [[ "$N" -lt 1 ]]; then
  echo "ERROR: no usable reference rows in ${REF_INPUTS}; required columns are target_species, wg_fasta_path, chrM_fasta_path." >&2
  exit 1
fi

mkdir -p "$OUTDIR" "$IN_HOUSE_SCORE_LOG_DIR"

echo "Submit/iterate in-house score jobs using ${REF_INPUTS}; NUMT candidates are reference-level chrM-vs-WG BLAST hits, not CRAM/sample-derived calls." >&2

if [[ "${MERGE_ONLY:-0}" == "1" ]]; then
  REF_INPUTS="$REF_INPUTS" OUTDIR="$OUTDIR" MERGED_IN_HOUSE_SCORE="$MERGED_IN_HOUSE_SCORE" IN_HOUSE_SCORE_LOG_DIR="$IN_HOUSE_SCORE_LOG_DIR" MERGE_ONLY=1 bash "$IN_HOUSE_SCORE_SCRIPT"
  bash preprocessing/scripts/merge_in_house_score.sh "$OUTDIR" "$MERGED_IN_HOUSE_SCORE"
  exit 0
fi

if [[ "${RUN_LOCAL:-0}" != "1" ]] && command -v sbatch >/dev/null 2>&1; then
  jid=$(sbatch --parsable --array="1-${N}%${MAX_CONCURRENT}" --export=ALL,REF_INPUTS="$REF_INPUTS",OUTDIR="$OUTDIR",MERGED_IN_HOUSE_SCORE="$MERGED_IN_HOUSE_SCORE",IN_HOUSE_SCORE_LOG_DIR="$IN_HOUSE_SCORE_LOG_DIR" "$IN_HOUSE_SCORE_SCRIPT")
  echo "Submitted in-house score array job ${jid} for ${N} reference rows." >&2
  merge_jid=$(sbatch --parsable --dependency="afterok:${jid}" --export=ALL,REF_INPUTS="$REF_INPUTS",OUTDIR="$OUTDIR",MERGED_IN_HOUSE_SCORE="$MERGED_IN_HOUSE_SCORE",IN_HOUSE_SCORE_LOG_DIR="$IN_HOUSE_SCORE_LOG_DIR",MERGE_ONLY=1 "$IN_HOUSE_SCORE_SCRIPT")
  echo "Submitted merge job ${merge_jid} after array ${jid}." >&2
else
  echo "sbatch unavailable or RUN_LOCAL=1; running ${N} reference rows locally in series." >&2
  for i in $(seq 1 "$N"); do
    SLURM_ARRAY_TASK_ID="$i" REF_INPUTS="$REF_INPUTS" OUTDIR="$OUTDIR" IN_HOUSE_SCORE_LOG_DIR="$IN_HOUSE_SCORE_LOG_DIR" bash "$IN_HOUSE_SCORE_SCRIPT"
  done
  REF_INPUTS="$REF_INPUTS" OUTDIR="$OUTDIR" MERGED_IN_HOUSE_SCORE="$MERGED_IN_HOUSE_SCORE" IN_HOUSE_SCORE_LOG_DIR="$IN_HOUSE_SCORE_LOG_DIR" MERGE_ONLY=1 bash "$IN_HOUSE_SCORE_SCRIPT"
  bash preprocessing/scripts/merge_in_house_score.sh "$OUTDIR" "$MERGED_IN_HOUSE_SCORE"
fi
