#!/usr/bin/env bash
#SBATCH --job-name=qc_preprocessing
#SBATCH --output=logs/qc_preprocessing/%x_%j.out
#SBATCH --error=logs/qc_preprocessing/%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash qc_analysis/scripts/run_qc_preprocessing.sh [--submit] [--sample SAMPLE] <step> [config/qc_preprocessing.yaml]
  sbatch qc_analysis/scripts/run_qc_preprocessing.sh <step> [config/qc_preprocessing.yaml]

Steps:
  collect_variant_calling_results  Collect and standardize variant-calling outputs only.
  discover_global_anchor           Discover reference-level global MSA anchors only.
  coordinate_liftover              Run coordinate liftover only.
  build_primate_codon_table        Build GenBank-first / MITOS2-fallback sample-level codon annotations.
  mitos2_annotation                 Run MITOS2 on unique final chrM references.
  codon_match                      Annotate lifted VCFs with codon matching.
  trna_match                       Annotate VCFs with tRNA matching.
  rrna_match                       Annotate VCFs with rRNA matching.
  all                              Run all preprocessing and downstream annotation steps.

Run modes:
  --submit                         Submit this wrapper to Slurm from a login/frontend node.
                                   Without --submit, bash runs the requested step immediately.

Environment overrides:
  PYTHON                           Python executable (default: python3).
  BIOPYTHON_USE_MODULE             Load the configured Biopython module for build_primate_codon_table (default: 1).
  BIOPYTHON_MODULE                 Biopython module to load (default: Biopython/1.83-foss-2022b).
  SAMPLE                           Optional sample name for coordinate_liftover.
  SLURM_PARTITION                  Optional partition/queue for --submit.
  SLURM_TIME                       Walltime for --submit (default: 24:00:00).
  SLURM_MEM                        Memory for --submit (default: 16G).
  SLURM_CPUS                       CPUs for --submit (default: 4).
  SLURM_LOG_DIR                    Log directory for --submit (default: logs/qc_preprocessing).
  SLURM_JOB_NAME                   Job name for --submit (default: qc_preprocessing_<step>).

Examples:
  bash qc_analysis/scripts/run_qc_preprocessing.sh --submit all config/qc_preprocessing.yaml
  bash qc_analysis/scripts/run_qc_preprocessing.sh --submit collect_variant_calling_results config/qc_preprocessing.yaml
  bash qc_analysis/scripts/run_qc_preprocessing.sh --submit discover_global_anchor config/qc_preprocessing.yaml
  bash qc_analysis/scripts/run_qc_preprocessing.sh --submit coordinate_liftover config/qc_preprocessing.yaml
  SAMPLE=SAMPLE_NAME bash qc_analysis/scripts/run_qc_preprocessing.sh --submit coordinate_liftover config/qc_preprocessing.yaml
  BIOPYTHON_MODULE=Biopython/1.83-foss-2022b bash qc_analysis/scripts/run_qc_preprocessing.sh build_primate_codon_table config/qc_preprocessing.yaml
  sbatch qc_analysis/scripts/run_qc_preprocessing.sh all config/qc_preprocessing.yaml
USAGE
}

SUBMIT_TO_SLURM=0
if [[ "${1:-}" == "--submit" ]]; then
  SUBMIT_TO_SLURM=1
  shift
fi

if [[ "${1:-}" == "--sample" ]]; then
  SAMPLE="$2"
  export SAMPLE
  shift 2
fi

if [[ $# -lt 1 || $# -gt 2 ]]; then
  usage >&2
  exit 2
fi

STEP="$1"
CONFIG="${2:-config/qc_preprocessing.yaml}"

case "$STEP" in
  -h|--help|help)
    usage
    exit 0
    ;;
  collect_variant_calling_results|discover_global_anchor|coordinate_liftover|build_primate_codon_table|mitos2_annotation|codon_match|trna_match|rrna_match|all)
    ;;
  *)
    echo "ERROR: unknown step: $STEP" >&2
    usage >&2
    exit 2
    ;;
esac

if [[ ! -s "$CONFIG" ]]; then
  echo "ERROR: missing or empty config file: $CONFIG" >&2
  exit 1
fi

submit_to_slurm() {
  if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    echo "ERROR: --submit was requested from inside an existing Slurm job (${SLURM_JOB_ID})." >&2
    exit 2
  fi
  if ! command -v sbatch >/dev/null 2>&1; then
    echo "ERROR: --submit requires sbatch on PATH. Run with bash without --submit to execute immediately." >&2
    exit 127
  fi

  local script_path log_dir job_name
  script_path="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
  log_dir="${SLURM_LOG_DIR:-logs/qc_preprocessing}"
  job_name="${SLURM_JOB_NAME:-qc_preprocessing_${STEP}}"
  mkdir -p "$log_dir"

  local sbatch_args=(
    --job-name="$job_name"
    --output="${log_dir}/%x_%j.out"
    --error="${log_dir}/%x_%j.err"
    --time="${SLURM_TIME:-24:00:00}"
    --cpus-per-task="${SLURM_CPUS:-4}"
    --mem="${SLURM_MEM:-16G}"
  )
  if [[ -n "${SLURM_PARTITION:-}" ]]; then
    sbatch_args+=(--partition="$SLURM_PARTITION")
  fi

  echo "[qc_preprocessing] Submitting ${STEP} to Slurm with config: ${CONFIG}" >&2
  sbatch "${sbatch_args[@]}" "$script_path" "$STEP" "$CONFIG"
}

if [[ "$SUBMIT_TO_SLURM" == "1" ]]; then
  submit_to_slurm
  exit 0
fi

PYTHON="${PYTHON:-python3}"
COLLECT_SCRIPT="qc_analysis/scripts/collect_variant_calling_results.py"
LIFTOVER_SCRIPT="qc_analysis/scripts/run_coordinate_liftover.py"
CODON_SCRIPT="qc_analysis/scripts/run_codon_match.py"
CODON_TABLE_SCRIPT="qc_analysis/scripts/build_primate_codon_table.py"
MITOS2_SCRIPT="qc_analysis/scripts/run_mitos2_annotation.py"
TRNA_SCRIPT="qc_analysis/scripts/run_trna_match.py"
RRNA_SCRIPT="qc_analysis/scripts/run_rrna_match.py"
GLOBAL_ANCHOR_SCRIPT="qc_analysis/scripts/discover_global_liftover_anchor.py"

# Read the small, optional environment.biopython section without depending on
# PyYAML (Biopython must be available before the build script can run).
configured_biopython_value() {
  local requested_key="$1"
  awk -v requested_key="$requested_key" '
    function indent(line) { match(line, /^[[:space:]]*/); return RLENGTH }
    function trim(value) { sub(/^[[:space:]]+/, "", value); sub(/[[:space:]]+$/, "", value); return value }
    {
      line = $0
      sub(/[[:space:]]*#.*/, "", line)
      if (line !~ /[^[:space:]]/) next
      level = indent(line)
      content = trim(line)

      if (content == "environment:") { environment_indent = level; in_environment = 1; in_biopython = 0; next }
      if (in_environment && level <= environment_indent) { in_environment = 0; in_biopython = 0 }
      if (in_environment && content == "biopython:") { biopython_indent = level; in_biopython = 1; next }
      if (in_biopython && level <= biopython_indent) in_biopython = 0
      if (in_biopython && content ~ ("^" requested_key ":[[:space:]]*")) {
        sub("^" requested_key ":[[:space:]]*", "", content)
        print trim(content)
        exit
      }
    }
  ' "$CONFIG"
}

if [[ -z "${BIOPYTHON_USE_MODULE+x}" ]]; then
  configured_use_module="$(configured_biopython_value use_module)"
  case "${configured_use_module,,}" in
    0|false|no) BIOPYTHON_USE_MODULE=0 ;;
    *) BIOPYTHON_USE_MODULE=1 ;;
  esac
fi
BIOPYTHON_MODULE="${BIOPYTHON_MODULE:-$(configured_biopython_value module_load)}"
BIOPYTHON_MODULE="${BIOPYTHON_MODULE:-Biopython/1.83-foss-2022b}"

run_collect_variant_calling_results() {
  echo "[qc_preprocessing] Running collect_variant_calling_results with config: ${CONFIG}" >&2
  "$PYTHON" "$COLLECT_SCRIPT" --config "$CONFIG"
}

run_discover_global_anchor() {
  echo "[qc_preprocessing] Running discover_global_anchor with config: ${CONFIG}" >&2
  echo "[qc_preprocessing] MAFFT environment preflight:" >&2
  "$PYTHON" "$GLOBAL_ANCHOR_SCRIPT" --config "$CONFIG" --check-environment | while IFS= read -r line; do
    echo "[qc_preprocessing] ${line}" >&2
  done
  "$PYTHON" "$GLOBAL_ANCHOR_SCRIPT" --config "$CONFIG"
}

run_coordinate_liftover() {
  echo "[qc_preprocessing] Running coordinate_liftover with config: ${CONFIG}" >&2
  local cmd=("$PYTHON" "$LIFTOVER_SCRIPT" --config "$CONFIG")
  if [[ -n "${SAMPLE:-}" ]]; then
    cmd+=(--sample "$SAMPLE")
  fi
  "${cmd[@]}"
}

run_mitos2_annotation() {
  echo "[qc_preprocessing] Running mitos2_annotation with config: ${CONFIG}" >&2
  local cmd=("$PYTHON" "$MITOS2_SCRIPT" --config "$CONFIG")
  [[ -n "${SAMPLE:-}" ]] && cmd+=(--sample "$SAMPLE")
  "${cmd[@]}"
}

run_build_primate_codon_table() {
  echo "[qc_preprocessing] Running build_primate_codon_table with config: ${CONFIG}" >&2
  if [[ "${BIOPYTHON_USE_MODULE}" == "1" ]]; then
    echo "[qc_preprocessing] Loading Biopython module: ${BIOPYTHON_MODULE}" >&2
    if command -v module >/dev/null 2>&1; then
      module load "${BIOPYTHON_MODULE}"
    elif [[ -f /etc/profile.d/modules.sh ]]; then
      # module is commonly initialized only for login shells on HPC systems.
      source /etc/profile.d/modules.sh
      module load "${BIOPYTHON_MODULE}"
    else
      echo "WARNING: BIOPYTHON_USE_MODULE=1 but module command is unavailable." >&2
    fi
  fi

  if ! "$PYTHON" - <<'PY'
from Bio import Entrez, SeqIO
print("Biopython import OK")
PY
  then
    echo "ERROR: Biopython is not importable after loading the configured module." >&2
    echo "Tried module: ${BIOPYTHON_MODULE}" >&2
    echo "Please check the HPC module name or set BIOPYTHON_USE_MODULE=0 if using a Python environment that already has Biopython." >&2
    exit 1
  fi

  local cmd=("$PYTHON" "$CODON_TABLE_SCRIPT" --config "$CONFIG")
  [[ -n "${SAMPLE:-}" ]] && cmd+=(--sample "$SAMPLE")
  "${cmd[@]}"
}

run_annotation() {
  local name="$1" script="$2"
  echo "[qc_preprocessing] Running ${name} with config: ${CONFIG}" >&2
  local cmd=("$PYTHON" "$script" --config "$CONFIG")
  [[ -n "${SAMPLE:-}" ]] && cmd+=(--sample "$SAMPLE")
  "${cmd[@]}"
}

case "$STEP" in
  collect_variant_calling_results)
    run_collect_variant_calling_results
    ;;
  discover_global_anchor)
    run_discover_global_anchor
    ;;
  coordinate_liftover)
    run_coordinate_liftover
    ;;
  mitos2_annotation) run_mitos2_annotation ;;
  build_primate_codon_table) run_build_primate_codon_table ;;
  codon_match) run_annotation codon_match "$CODON_SCRIPT" ;;
  trna_match) run_annotation trna_match "$TRNA_SCRIPT" ;;
  rrna_match) run_annotation rrna_match "$RRNA_SCRIPT" ;;
  all)
    run_collect_variant_calling_results
    run_discover_global_anchor
    run_coordinate_liftover
    run_mitos2_annotation
    run_build_primate_codon_table
    run_annotation codon_match "$CODON_SCRIPT"
    run_annotation trna_match "$TRNA_SCRIPT"
    run_annotation rrna_match "$RRNA_SCRIPT"
    ;;
esac
