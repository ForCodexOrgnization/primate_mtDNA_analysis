#!/usr/bin/env bash
#SBATCH --job-name=preprocessing
#SBATCH --output=logs/preprocessing/%x_%j.out
#SBATCH --error=logs/preprocessing/%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash preprocessing/scripts/run_preprocessing.sh [--submit] <step> [config/preprocessing_paths.yaml]
  sbatch preprocessing/scripts/run_preprocessing.sh <step> [config/preprocessing_paths.yaml]

Run modes:
  --submit                  Submit this script to Slurm from a login/frontend node.
                            Without --submit, bash runs the step immediately.

Steps:
  reference_discovery        Run reference discovery only.
  reference_materialization  Build/materialize reviewed reference manifests.
  in_house_score             Launch/iterate in-house score jobs.
  variant_references         Build variant-calling reference packages.
  post_reference_review      Run reference_materialization, in_house_score, variant_references.
  all_steps                  Run all preprocessing stages using the raw discovery summary without manual review.
  reports                    Render preprocessing Quarto reports.
  all                        Run reference_discovery only, then stop for manual review.

HPC environment config keys:
  environment_setup_script  Optional shell script to source before running a step.
                            Use it for module load / conda activate commands.
  rscript_command           Rscript executable name/path (default: Rscript).
  python_command            Python executable name/path (default: python3).
  wget_command              wget executable name/path (default: wget).
  samtools_command          samtools executable name/path (default: samtools).
  bwa_command              bwa executable name/path (default: bwa).
  gatk_command             gatk executable name/path (default: gatk).
  variant_reference_threads
                            Reference packages to build in parallel (default: 1).
  curl_command             curl executable name/path (default: curl).
  efetch_command           efetch executable name/path (default: efetch).
  reference_discovery_threads
                            Species-level worker threads for reference discovery (default: 1).

Slurm submit environment overrides for --submit:
  SLURM_PARTITION           Optional partition/queue name.
  SLURM_TIME                Walltime passed to sbatch (default: 24:00:00).
  SLURM_MEM                 Memory passed to sbatch (default: 8G).
  SLURM_CPUS                CPUs passed to sbatch (default: 1).
  SLURM_LOG_DIR             Log directory (default: logs/preprocessing).
  SLURM_JOB_NAME            Job name prefix (default: preprocessing_<step>).

Examples:
  bash preprocessing/scripts/run_preprocessing.sh --submit reference_discovery config/preprocessing_paths.yaml
  sbatch preprocessing/scripts/run_preprocessing.sh post_reference_review config/preprocessing_paths.yaml
  bash preprocessing/scripts/run_preprocessing.sh reports config/preprocessing_paths.yaml

Important:
  The all step intentionally stops after reference discovery because
  species_reference_chrM_summary.tsv must be manually reviewed before downstream steps.
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
CONFIG="${2:-config/preprocessing_paths.yaml}"

case "$STEP" in
  -h|--help|help)
    usage
    exit 0
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
  log_dir="${SLURM_LOG_DIR:-logs/preprocessing}"
  job_name="${SLURM_JOB_NAME:-preprocessing_${STEP}}"
  mkdir -p "$log_dir"

  local sbatch_args=(
    --job-name="$job_name"
    --output="${log_dir}/%x_%j.out"
    --error="${log_dir}/%x_%j.err"
    --time="${SLURM_TIME:-24:00:00}"
    --cpus-per-task="${SLURM_CPUS:-1}"
    --mem="${SLURM_MEM:-8G}"
  )
  if [[ -n "${SLURM_PARTITION:-}" ]]; then
    sbatch_args+=(--partition="$SLURM_PARTITION")
  fi

  echo "[preprocessing] Submitting ${STEP} to Slurm with config: ${CONFIG}" >&2
  sbatch "${sbatch_args[@]}" "$script_path" "$STEP" "$CONFIG"
}

if [[ "$SUBMIT_TO_SLURM" == "1" ]]; then
  submit_to_slurm
  exit 0
fi

config_get() {
  local key="$1" default="${2:-}"
  local value
  value=$(awk -F: -v key="$key" '
    $1 == key {
      sub(/^[[:space:]]+/, "", $2)
      sub(/[[:space:]]+$/, "", $2)
      gsub(/^"|"$/, "", $2)
      print $2
      exit
    }
  ' "$CONFIG")
  if [[ -n "$value" ]]; then
    printf '%s\n' "$value"
  else
    printf '%s\n' "$default"
  fi
}

SPECIES_TABLE=$(config_get species_table "data/metadata/all_species_list.txt")
MITO_FASTA=$(config_get mito_fasta "/path/to/mitochondrion.1.1.genomic.fna.gz")
TREE_NEWICK=$(config_get tree_newick "/path/to/primate_tree.nwk")
REFERENCE_DISCOVERY_OUTDIR=$(config_get reference_discovery_outdir "results/preprocessing/reference_discovery")
EMAIL=$(config_get email "your_email@yale.edu")
MAX_NEAREST=$(config_get max_nearest "200")
DELAY=$(config_get delay "0.34")
REFERENCE_DISCOVERY_THREADS=$(config_get reference_discovery_threads "1")

SPECIES_REFERENCE_SUMMARY=$(config_get species_reference_summary "data/metadata/species_reference_chrM_summary.tsv")
REFERENCE_MATERIALIZATION_RESULTS_MANIFEST=$(config_get reference_materialization_results_manifest "results/preprocessing/reference_materialization/reference_materialization_manifest.tsv")
REFERENCE_MATERIALIZATION_MANIFEST=$(config_get reference_materialization_manifest "references/manifests/reference_materialization_manifest.tsv")
IN_HOUSE_SCORE_REFERENCE_INPUTS=$(config_get in_house_score_reference_inputs "references/manifests/in_house_score_reference_inputs.tsv")
MERGED_IN_HOUSE_SCORE=$(config_get merged_in_house_score "results/preprocessing/in_house_score/merged_in_house_score.tsv")
NUMT_TARGET_CHRM_COV=$(config_get numt_target_chrm_cov "0.95")
MASK_REF_TYPES=$(config_get mask_ref_types "#C-likely_comp,#C-Ambiguous,#A")
A_MASK_MODE=$(config_get a_mask_mode "mask_if_requested")
VARIANT_CALLING_REFERENCE_OUT_ROOT=$(config_get variant_calling_reference_out_root "references/variant_calling")
VARIANT_REFERENCE_THREADS=$(config_get variant_reference_threads "1")

ENVIRONMENT_SETUP_SCRIPT=$(config_get environment_setup_script "")
RSCRIPT_COMMAND=$(config_get rscript_command "Rscript")
PYTHON_COMMAND=$(config_get python_command "python3")
WGET_COMMAND=$(config_get wget_command "wget")
SAMTOOLS_COMMAND=$(config_get samtools_command "samtools")
BWA_COMMAND=$(config_get bwa_command "bwa")
GATK_COMMAND=$(config_get gatk_command "gatk")
CURL_COMMAND=$(config_get curl_command "curl")
EFETCH_COMMAND=$(config_get efetch_command "efetch")

source_environment_setup() {
  if [[ -z "$ENVIRONMENT_SETUP_SCRIPT" ]]; then
    return 0
  fi
  if [[ ! -s "$ENVIRONMENT_SETUP_SCRIPT" ]]; then
    echo "ERROR: environment_setup_script is configured but missing or empty: $ENVIRONMENT_SETUP_SCRIPT" >&2
    exit 1
  fi
  echo "[preprocessing] Sourcing HPC environment setup: $ENVIRONMENT_SETUP_SCRIPT" >&2
  # shellcheck source=/dev/null
  source "$ENVIRONMENT_SETUP_SCRIPT"
}

require_command() {
  local command_name="$1" description="${2:-$1}"
  if ! command -v "$command_name" >/dev/null 2>&1; then
    cat >&2 <<EOFMSG
ERROR: required command not found for preprocessing: ${command_name} (${description})
Configure environment_setup_script in ${CONFIG} to load your HPC modules/conda env,
or set the corresponding *_command config key to an executable path.
EOFMSG
    exit 127
  fi
}


require_nonempty_file() {
  local path="$1" description="${2:-file}"
  if [[ ! -s "$path" ]]; then
    echo "ERROR: missing or empty ${description}: ${path}" >&2
    exit 1
  fi
}

check_reference_discovery_environment() {
  require_command "$PYTHON_COMMAND" "Python interpreter for reference discovery"
}

check_reference_materialization_environment() {
  require_command "$RSCRIPT_COMMAND" "Rscript for manifest generation"
  require_command "$PYTHON_COMMAND" "Python interpreter for materialization helpers"
  require_command "$WGET_COMMAND" "wget for reference downloads"
  require_command "$SAMTOOLS_COMMAND" "samtools for FASTA indexing/extraction"
  if ! command -v "$EFETCH_COMMAND" >/dev/null 2>&1; then
    require_command "$CURL_COMMAND" "curl fallback for NCBI efetch downloads"
  fi
}

check_variant_references_environment() {
  require_command "$PYTHON_COMMAND" "Python interpreter for variant-reference package generation"
  require_command "$SAMTOOLS_COMMAND" "samtools for variant-reference FASTA indexing"
  require_command "$BWA_COMMAND" "bwa for variant-reference BWA indexing"
  require_command "$GATK_COMMAND" "GATK for variant-reference sequence dictionaries"
}

source_environment_setup

run_reference_discovery() {
  check_reference_discovery_environment
  echo "[preprocessing] Running reference discovery with species table: $SPECIES_TABLE" >&2
  SPECIES_TABLE="$SPECIES_TABLE" \
  MITO_FASTA="$MITO_FASTA" \
  TREE_NEWICK="$TREE_NEWICK" \
  OUTDIR="$REFERENCE_DISCOVERY_OUTDIR" \
  EMAIL="$EMAIL" \
  MAX_NEAREST="$MAX_NEAREST" \
  DELAY="$DELAY" \
  REFERENCE_DISCOVERY_THREADS="$REFERENCE_DISCOVERY_THREADS" \
  PYTHON_COMMAND="$PYTHON_COMMAND" \
    bash preprocessing/scripts/run_reference_discovery.sh
  cat >&2 <<EOFMSG

[preprocessing] Reference discovery complete.
[preprocessing] Review: ${REFERENCE_DISCOVERY_OUTDIR}/species_reference_chrM_summary.tsv
[preprocessing] After review, copy or symlink the reviewed manifest to:
  ${SPECIES_REFERENCE_SUMMARY}
[preprocessing] Then run:
  bash preprocessing/scripts/run_preprocessing.sh post_reference_review ${CONFIG}
EOFMSG
}

use_raw_discovery_summary_without_review() {
  local raw_summary="${REFERENCE_DISCOVERY_OUTDIR}/species_reference_chrM_summary.tsv"
  if [[ ! -s "$raw_summary" ]]; then
    echo "ERROR: raw reference-discovery summary is missing or empty: $raw_summary" >&2
    exit 1
  fi
  mkdir -p "$(dirname "$SPECIES_REFERENCE_SUMMARY")"
  cp "$raw_summary" "$SPECIES_REFERENCE_SUMMARY"
  cat >&2 <<EOFMSG
[preprocessing] WARNING: skipped manual reference review.
[preprocessing] Copied raw discovery summary to reviewed-manifest path:
  ${SPECIES_REFERENCE_SUMMARY}
[preprocessing] Downstream reference choices may need later manual correction.
EOFMSG
}

run_reference_materialization() {
  echo "[preprocessing] Building reference materialization manifest from: $SPECIES_REFERENCE_SUMMARY" >&2
  require_nonempty_file "$SPECIES_REFERENCE_SUMMARY" "reviewed species reference summary"
  check_reference_materialization_environment
  "$RSCRIPT_COMMAND" preprocessing/scripts/build_reference_materialization_manifest.R \
    "$SPECIES_REFERENCE_SUMMARY" \
    "$REFERENCE_MATERIALIZATION_RESULTS_MANIFEST" \
    "$REFERENCE_MATERIALIZATION_MANIFEST"
  echo "[preprocessing] Materializing references from: $REFERENCE_MATERIALIZATION_MANIFEST" >&2
  PYTHON_COMMAND="$PYTHON_COMMAND" \
  WGET_COMMAND="$WGET_COMMAND" \
  SAMTOOLS_COMMAND="$SAMTOOLS_COMMAND" \
  CURL_COMMAND="$CURL_COMMAND" \
  EFETCH_COMMAND="$EFETCH_COMMAND" \
    bash preprocessing/scripts/materialize_references.sh "$REFERENCE_MATERIALIZATION_MANIFEST"
}

run_in_house_score() {
  echo "[preprocessing] Running in-house score step with: $IN_HOUSE_SCORE_REFERENCE_INPUTS" >&2
  local score_outdir
  score_outdir="$(dirname "$MERGED_IN_HOUSE_SCORE")"
  REF_INPUTS="$IN_HOUSE_SCORE_REFERENCE_INPUTS" \
  OUTDIR="$score_outdir" \
  MERGED_IN_HOUSE_SCORE="$MERGED_IN_HOUSE_SCORE" \
  NUMT_TARGET_CHRM_COV="$NUMT_TARGET_CHRM_COV" \
  MASK_REF_TYPES="$MASK_REF_TYPES" \
  A_MASK_MODE="$A_MASK_MODE" \
    bash preprocessing/scripts/run_in_house_score_array.sh
}

run_variant_references() {
  echo "[preprocessing] Building variant-calling reference packages under: $VARIANT_CALLING_REFERENCE_OUT_ROOT" >&2
  require_nonempty_file "$IN_HOUSE_SCORE_REFERENCE_INPUTS" "in-house score reference inputs"
  require_nonempty_file "$MERGED_IN_HOUSE_SCORE" "merged in-house score table"
  check_variant_references_environment
  REF_INPUTS="$IN_HOUSE_SCORE_REFERENCE_INPUTS" \
  SCORE="$MERGED_IN_HOUSE_SCORE" \
  OUT_ROOT="$VARIANT_CALLING_REFERENCE_OUT_ROOT" \
  MASK_REF_TYPES="$MASK_REF_TYPES" \
  VARIANT_REFERENCE_THREADS="$VARIANT_REFERENCE_THREADS" \
  PYTHON_COMMAND="$PYTHON_COMMAND" \
  SAMTOOLS_COMMAND="$SAMTOOLS_COMMAND" \
  BWA_COMMAND="$BWA_COMMAND" \
  GATK_COMMAND="$GATK_COMMAND" \
    bash preprocessing/scripts/run_variant_references_array.sh
}

run_reports() {
  bash preprocessing/scripts/render_preprocessing_reports.sh
}

case "$STEP" in
  reference_discovery)
    run_reference_discovery
    ;;
  reference_materialization)
    run_reference_materialization
    ;;
  in_house_score)
    run_in_house_score
    ;;
  variant_references)
    run_variant_references
    ;;
  post_reference_review)
    run_reference_materialization
    run_in_house_score
    run_variant_references
    ;;
  all_steps)
    run_reference_discovery
    use_raw_discovery_summary_without_review
    run_reference_materialization
    run_in_house_score
    run_variant_references
    ;;
  reports)
    run_reports
    ;;
  all)
    run_reference_discovery
    echo "[preprocessing] all stops here by design; manual reference review is required before downstream steps." >&2
    ;;
  -h|--help|help)
    usage
    ;;
  *)
    echo "ERROR: unknown preprocessing step: $STEP" >&2
    usage >&2
    exit 2
    ;;
esac
