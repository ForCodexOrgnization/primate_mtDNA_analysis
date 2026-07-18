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
  bash qc_analysis/scripts/run_qc_preprocessing.sh [--submit] <step> [config/qc_preprocessing.yaml]
  sbatch qc_analysis/scripts/run_qc_preprocessing.sh <step> [config/qc_preprocessing.yaml]

Steps:
  collect_variant_calling_results  Collect and standardize variant-calling outputs only.
  discover_global_anchor           Discover reference-level global MSA anchors only.
  coordinate_liftover              Run coordinate liftover only.
  all                              Run collect_variant_calling_results, discover_global_anchor, then coordinate_liftover.

Run modes:
  --submit                         Submit this wrapper to Slurm from a login/frontend node.
                                   Without --submit, bash runs the requested step immediately.

Environment overrides:
  PYTHON                           Python executable (default: python3).
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
  sbatch qc_analysis/scripts/run_qc_preprocessing.sh all config/qc_preprocessing.yaml
USAGE
}

SUBMIT_TO_SLURM=0
if [[ "${1:-}" == "--submit" ]]; then
  SUBMIT_TO_SLURM=1
  shift
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
  collect_variant_calling_results|discover_global_anchor|coordinate_liftover|all)
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
GLOBAL_ANCHOR_SCRIPT="qc_analysis/scripts/discover_global_liftover_anchor.py"

run_collect_variant_calling_results() {
  echo "[qc_preprocessing] Running collect_variant_calling_results with config: ${CONFIG}" >&2
  "$PYTHON" "$COLLECT_SCRIPT" --config "$CONFIG"
}

run_discover_global_anchor() {
  echo "[qc_preprocessing] Running discover_global_anchor with config: ${CONFIG}" >&2
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
  all)
    run_collect_variant_calling_results
    run_discover_global_anchor
    run_coordinate_liftover
    ;;
esac
