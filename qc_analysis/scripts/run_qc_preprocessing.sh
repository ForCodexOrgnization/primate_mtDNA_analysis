#!/usr/bin/env bash
#SBATCH --job-name=qc_preprocessing
#SBATCH --output=logs/qc_preprocessing/%x_%j.out
#SBATCH --error=logs/qc_preprocessing/%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash qc_analysis/scripts/run_qc_preprocessing.sh [--submit] [--sample SAMPLE] <step> [config/qc_preprocessing.yaml]
  sbatch qc_analysis/scripts/run_qc_preprocessing.sh <step> [config/qc_preprocessing.yaml]

Steps:
  collect_variant_calling_results  Collect and standardize variant-calling outputs only.
  intraspecies_contamination       Write original-coordinate sample contamination QC report.
  sample_variant_filtering         Write the five-criterion biological sample QC report.
  discover_global_anchor           Discover reference-level global MSA anchors only.
  coordinate_liftover              Run coordinate liftover only.
  build_primate_homo_background   Consolidate orthology QC and build temporary homoplasmic background.
  human_contamination              Screen annotated human-coordinate alleles with primate background correction.
  build_primate_codon_table        Build the optional independent GenBank benchmark (not production).
  mitos2_prepare_tasks              Write one Slurm-array task per target-species chrM FASTA.
  mitos2_merge                      Merge completed per-reference MITOS2 raw outputs.
  mitos2_annotation                 Run MITOS2 sequentially on target-species chrM FASTAs.
  codon_match                      Annotate lifted VCFs with codon matching.
  codon_match_validate             Validate and index codon inputs without reading VCFs.
  codon_match_merge                Atomically merge completed per-sample summaries.
  build_trna_indexes               Build one tRNAscan index per unique reference.
  trna_match                       Annotate VCFs with tRNA matching.
  trna_match_merge                 Atomically merge per-sample tRNA summaries.
  trna_gene_qc                     Compare lifted source and human tRNA genes.
  rrna_match                       Annotate VCFs with rRNA matching.
  rrna_match_merge                 Atomically merge per-sample rRNA summaries.
  final_filter                     Combine QC reports and materialize final filtered files.
  all                              Run all preprocessing and downstream annotation steps.

Run modes:
  --submit                         Submit this wrapper to Slurm from a login/frontend node.
                                   Without --submit, bash runs the requested step immediately.

Environment overrides:
  PYTHON                           Python executable (default: python3).
  BIOPYTHON_USE_MODULE             Load the configured Biopython module for build_primate_codon_table (default: 1).
  BIOPYTHON_MODULE                 Biopython module to load (default: Biopython/1.83-foss-2022b).
  CODON_TABLE_WORKERS              Positive worker count for internal codon-table preparation.
  CODON_TABLE_DOWNLOAD_WORKERS     Positive concurrent Entrez download count (rate limited).
  SAMPLE                           Optional sample name for sample-level steps.
  SLURM_ARRAY_CONCURRENCY          Maximum concurrent array tasks (default: 20).
  AUTO_SUBMIT_MERGE                Submit codon/tRNA/rRNA/MITOS2 merge afterok (default: true).
  SKIP_COMPLETED / FORCE_RERUN     Resume controls (defaults: true / false).
  SLURM_PARTITION                  Optional partition/queue for --submit.
  SLURM_TIME                       Walltime for --submit (default: 24:00:00).
  SLURM_MEM                        Memory for --submit (default: 16G).
  SLURM_CPUS                       CPUs for --submit (default: 4).
  SLURM_LOG_DIR                    Override the step-specific scheduler log directory.
  SLURM_JOB_NAME                   Job name for --submit (default: qc_preprocessing_<step>).

Options: --array-concurrency N, --task-file PATH, --prepare-only,
         --prepare-retry, --dry-run-submit. Direct mode remains sequential;
         --submit is recommended for HPC execution.

New array metadata is stored in <step output>/job_arrays and logs default to
<step output>/logs/job_arrays. Optional paths.job_array_dir and paths.log_dir
settings override those derived paths. Historical results/qc/job_arrays files
are left untouched; move them manually only after confirming no old job uses them.

Examples:
  bash qc_analysis/scripts/run_qc_preprocessing.sh --submit all config/qc_preprocessing.yaml
  bash qc_analysis/scripts/run_qc_preprocessing.sh --submit collect_variant_calling_results config/qc_preprocessing.yaml
  bash qc_analysis/scripts/run_qc_preprocessing.sh --submit discover_global_anchor config/qc_preprocessing.yaml
  bash qc_analysis/scripts/run_qc_preprocessing.sh --submit coordinate_liftover config/qc_preprocessing.yaml
  bash qc_analysis/scripts/run_qc_preprocessing.sh --submit codon_match_validate config/qc_preprocessing.yaml
  bash qc_analysis/scripts/run_qc_preprocessing.sh --submit codon_match config/qc_preprocessing.yaml
  bash qc_analysis/scripts/run_qc_preprocessing.sh trna_match_merge config/qc_preprocessing.yaml
  bash qc_analysis/scripts/run_qc_preprocessing.sh rrna_match_merge config/qc_preprocessing.yaml
  SLURM_ARRAY_CONCURRENCY=40 bash qc_analysis/scripts/run_qc_preprocessing.sh --submit codon_match config/qc_preprocessing.yaml
  bash qc_analysis/scripts/run_qc_preprocessing.sh --submit --sample SAMPLE_NAME codon_match config/qc_preprocessing.yaml
  bash qc_analysis/scripts/run_qc_preprocessing.sh --dry-run-submit codon_match config/qc_preprocessing.yaml
  SAMPLE=SAMPLE_NAME bash qc_analysis/scripts/run_qc_preprocessing.sh --submit coordinate_liftover config/qc_preprocessing.yaml
  BIOPYTHON_MODULE=Biopython/1.83-foss-2022b bash qc_analysis/scripts/run_qc_preprocessing.sh build_primate_codon_table config/qc_preprocessing.yaml
  SLURM_CPUS=8 CODON_TABLE_WORKERS=8 bash qc_analysis/scripts/run_qc_preprocessing.sh --submit build_primate_codon_table config/qc_preprocessing.yaml
  sbatch qc_analysis/scripts/run_qc_preprocessing.sh all config/qc_preprocessing.yaml
USAGE
}

SUBMIT_TO_SLURM=0; DRY_RUN_SUBMIT=0; PREPARE_ONLY=0; PREPARE_RETRY=0; ARRAY_TASK_MODE=0
SAMPLE="${SAMPLE:-}"; TASK_FILE=""; ARRAY_CONCURRENCY="${SLURM_ARRAY_CONCURRENCY:-20}"
while [[ "${1:-}" == --* ]]; do
  case "$1" in
    --submit) SUBMIT_TO_SLURM=1; shift ;;
    --dry-run-submit) SUBMIT_TO_SLURM=1; DRY_RUN_SUBMIT=1; shift ;;
    --prepare-only) PREPARE_ONLY=1; shift ;;
    --prepare-retry) PREPARE_RETRY=1; shift ;;
    --sample)
      [[ $# -ge 2 ]] || { echo "ERROR: --sample requires a value" >&2; exit 2; }
      SAMPLE="$2"
      shift 2
      ;;
    --task-file)
      [[ $# -ge 2 ]] || { echo "ERROR: --task-file requires a value" >&2; exit 2; }
      TASK_FILE="$2"
      shift 2
      ;;
    --array-concurrency)
      [[ $# -ge 2 ]] || { echo "ERROR: --array-concurrency requires a value" >&2; exit 2; }
      ARRAY_CONCURRENCY="$2"
      shift 2
      ;;
    --array-task)
      ARRAY_TASK_MODE=1
      shift
      ;;
    --help) usage; exit 0 ;;
    *)
      echo "ERROR: unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done
[[ "$ARRAY_CONCURRENCY" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: array concurrency must be a positive integer: $ARRAY_CONCURRENCY" >&2; exit 2; }
[[ $# -ge 1 && $# -le 2 ]] || { usage >&2; exit 2; }
STEP="$1"; CONFIG="${2:-config/qc_preprocessing.yaml}"
case "$STEP" in
 collect_variant_calling_results|discover_global_anchor|coordinate_liftover|build_primate_homo_background|human_contamination|build_primate_codon_table|compare_genbank_mitos2|mitos2_prepare_tasks|mitos2_merge|mitos2_annotation|codon_match|codon_match_validate|codon_match_merge|build_trna_indexes|trna_match|trna_match_merge|trna_gene_qc|rrna_match|rrna_match_merge|intraspecies_contamination|sample_variant_filtering|final_filter|all) ;;
 -h|--help|help) usage; exit 0;; *) echo "ERROR: unknown step: $STEP" >&2; exit 2;; esac
[[ -s "$CONFIG" ]] || { echo "ERROR: missing or empty config file: $CONFIG" >&2; exit 1; }
export SAMPLE

classify_step() { case "$1" in coordinate_liftover|codon_match|trna_match|rrna_match) echo sample;; mitos2_annotation) echo reference;; *) echo singleton;; esac; }
build_array_expression() { local n="$1" kind="$2"; if [[ "$kind" == singleton || "$n" == 1 ]]; then echo 1-1; else echo "1-${n}%${ARRAY_CONCURRENCY}"; fi; }
trna_setting() {
 awk -v key="$1" '
   function indent(s){match(s,/^[[:space:]]*/);return RLENGTH}
   function trim(s){sub(/^[[:space:]]+/,"",s);sub(/[[:space:]]+$/, "",s);return s}
   { line=$0; sub(/[[:space:]]*#.*/,"",line); level=indent(line); value=trim(line)
     if(value=="trna_match:"){in_trna=1;trna_indent=level;next}
     if(in_trna && level<=trna_indent && value!=""){in_trna=0}
     if(in_trna && value=="settings:"){in_settings=1;settings_indent=level;next}
     if(in_settings && level<=settings_indent && value!=""){in_settings=0}
     if(in_settings && value ~ ("^" key ":[[:space:]]*")){sub("^[^:]+:[[:space:]]*","",value);print value;exit}
   }' "$CONFIG"
}
resolve_step_resources() {
 local prefix=""; case "$1" in codon_match*) prefix=CODON_MATCH;; coordinate_liftover) prefix=LIFTOVER;; build_trna_indexes) prefix=TRNA_INDEX_BUILD;; trna_match*) prefix=TRNA_MATCH;; rrna_match*) prefix=RRNA_MATCH;; mitos2*) prefix=MITOS2;; esac
 local specific=""
 [[ -n "$prefix" ]] && { local vn="${prefix}_SLURM_TIME"; specific="${!vn:-}"; }; RES_TIME="${specific:-${SLURM_TIME:-24:00:00}}"
 specific=""; [[ -n "$prefix" ]] && { local vn="${prefix}_SLURM_MEM"; specific="${!vn:-}"; }; RES_MEM="${specific:-${SLURM_MEM:-16G}}"
 specific=""; [[ -n "$prefix" ]] && { local vn="${prefix}_SLURM_CPUS"; specific="${!vn:-}"; }; RES_CPUS="${specific:-${SLURM_CPUS:-4}}"
 if [[ "$1" == build_trna_indexes && -z "$specific" && -z "${SLURM_CPUS:-}" ]]; then
   local workers threads
   workers="$(trna_setting index_build_workers)"; workers="${workers:-1}"
   threads="$(trna_setting trnascan_threads)"; threads="${threads:-1}"
   [[ "$workers" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: index_build_workers must be >= 1: $workers" >&2; exit 2; }
   [[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: trnascan_threads must be >= 1: $threads" >&2; exit 2; }
   RES_CPUS=$((workers * threads))
 fi
}
resolve_array_item() {
 [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || { echo 'ERROR: SLURM_ARRAY_TASK_ID is required in an array task' >&2; exit 2; }
 [[ "$SLURM_ARRAY_TASK_ID" =~ ^[1-9][0-9]*$ && -f "$TASK_FILE" ]] || { echo "ERROR: invalid task index or missing task file: $TASK_FILE" >&2; exit 2; }
 local count; count=$(awk 'NF{n++} END{print n+0}' "$TASK_FILE")
 (( SLURM_ARRAY_TASK_ID <= count )) || { echo "ERROR: task index $SLURM_ARRAY_TASK_ID exceeds $count" >&2; exit 2; }
 ITEM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASK_FILE"); [[ -n "$ITEM" ]] || { echo 'ERROR: selected array item is empty' >&2; exit 2; }
 printf '[qc_preprocessing] step=%s array_job_id=%s array_task_id=%s selected_item=%s config=%s hostname=%s start_time=%s\n' "$STEP" "${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-unknown}}" "$SLURM_ARRAY_TASK_ID" "$ITEM" "$CONFIG" "$(hostname)" "$(date -u +%FT%TZ)" >&2
  # Task-file entries only represent samples for sample-classified steps.
  # Singleton entries are step names and reference entries are MITOS2 worker
  # keys, neither of which may leak into a downstream --sample argument.
  if [[ "$(classify_step "$STEP")" == sample ]]; then
    SAMPLE="$ITEM"
  else
    SAMPLE=""
  fi
}
prepare_task_manifest() {
 local args=(python3 qc_analysis/scripts/qc_array_manifest.py "$1" "$CONFIG")
 [[ -n "$SAMPLE" ]] && args+=(--sample "$SAMPLE")
 [[ "${FORCE_RERUN:-false}" == true || "${SKIP_COMPLETED:-true}" != true ]] && args+=(--force)
 [[ "$PREPARE_RETRY" == 1 ]] && args+=(--retry)
 local output; output=$("${args[@]}"); TASK_FILE=$(printf '%s\n' "$output"|sed -n 's/^TASK_FILE=//p'); MANIFEST=$(printf '%s\n' "$output"|sed -n 's/^MANIFEST=//p'); TASK_COUNT=$(printf '%s\n' "$output"|sed -n 's/^COUNT=//p')
 OUTPUT_DIR=$(printf '%s\n' "$output"|sed -n 's/^OUTPUT_DIR=//p'); JOB_ARRAY_DIR=$(printf '%s\n' "$output"|sed -n 's/^JOB_ARRAY_DIR=//p'); CONFIG_LOG_DIR=$(printf '%s\n' "$output"|sed -n 's/^LOG_DIR=//p')
}
resolve_runtime_paths() {
 local output; output=$(python3 qc_analysis/scripts/qc_array_manifest.py "$1" "$CONFIG" --resolve-paths)
 OUTPUT_DIR=$(printf '%s\n' "$output"|sed -n 's/^OUTPUT_DIR=//p'); JOB_ARRAY_DIR=$(printf '%s\n' "$output"|sed -n 's/^JOB_ARRAY_DIR=//p'); CONFIG_LOG_DIR=$(printf '%s\n' "$output"|sed -n 's/^LOG_DIR=//p')
}
submit_array() {
 local step="$1" dependency="${2:-}"; STEP="$step"
 if [[ -n "$TASK_FILE" ]]; then
   [[ -f "$TASK_FILE" ]] || { echo "ERROR: task file does not exist: $TASK_FILE" >&2; exit 1; }
   TASK_COUNT=$(awk 'NF{n++} END{print n+0}' "$TASK_FILE")
   (( TASK_COUNT > 0 )) || { echo "ERROR: task file is empty: $TASK_FILE" >&2; exit 1; }
 else prepare_task_manifest "$step"; fi
 [[ -n "${OUTPUT_DIR:-}" ]] || resolve_runtime_paths "$step"
 local kind array;kind=$(classify_step "$step");array=$(build_array_expression "$TASK_COUNT" "$kind");resolve_step_resources "$step"
 local logs="${SLURM_LOG_DIR:-$CONFIG_LOG_DIR}";mkdir -p "$logs"
 local cmd=(sbatch --parsable --job-name="qc_preprocessing_${step}" --array="$array" --output="$logs/%A_%a.out" --error="$logs/%A_%a.err" --time="$RES_TIME" --mem="$RES_MEM" --cpus-per-task="$RES_CPUS")
 [[ -n "${SLURM_PARTITION:-}" ]] && cmd+=(--partition="$SLURM_PARTITION"); [[ -n "$dependency" ]] && cmd+=(--dependency="afterok:$dependency")
 cmd+=("$(readlink -f "${BASH_SOURCE[0]}")" --array-task --task-file "$TASK_FILE" "$step" "$CONFIG")
 printf '[qc_preprocessing] step=%s\n[qc_preprocessing] output_dir=%s\n[qc_preprocessing] task_file=%s\n[qc_preprocessing] manifest=%s\n[qc_preprocessing] logs=%s/%%A_%%a.{out,err}\n[qc_preprocessing] task_count=%s\n[qc_preprocessing] array=%s\n' "$step" "$OUTPUT_DIR" "$TASK_FILE" "${MANIFEST:-unknown}" "$logs" "$TASK_COUNT" "$array" >&2
 printf '[qc_preprocessing] concurrency=%s dependency=%s resources=%s,%s,%s\n' "$ARRAY_CONCURRENCY" "${dependency:-none}" "$RES_TIME" "$RES_MEM" "$RES_CPUS" >&2
 if [[ "$PREPARE_ONLY" == 1 ]]; then LAST_JOB_ID="prepared_${step}"; return; fi
 if [[ "$DRY_RUN_SUBMIT" == 1 ]]; then printf 'DRY RUN:';printf ' %q' "${cmd[@]}";printf '\n'; LAST_JOB_ID="dry_${step}"; else command -v sbatch >/dev/null || { echo 'ERROR: --submit requires sbatch on PATH' >&2;exit 127; };LAST_JOB_ID=$("${cmd[@]}");LAST_JOB_ID=${LAST_JOB_ID%%;*};fi
 local submission="${TASK_FILE%.tasks.txt}.submission.tsv" submitted_at; submitted_at=$(date -u +%FT%TZ)
 printf 'step\tjob_id\ttask_file\tmanifest\tlog_dir\tarray\tsubmitted_at\n%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$step" "$LAST_JOB_ID" "$TASK_FILE" "${MANIFEST:-}" "$logs" "$array" "$submitted_at" >"$submission"
}
submit_workflow() {
 local dep=""; SAMPLE=""; local steps=(collect_variant_calling_results intraspecies_contamination sample_variant_filtering discover_global_anchor coordinate_liftover mitos2_prepare_tasks mitos2_annotation mitos2_merge codon_match_validate codon_match codon_match_merge build_trna_indexes trna_match trna_match_merge rrna_match rrna_match_merge build_primate_homo_background human_contamination final_filter)
 for s in "${steps[@]}"; do
   # Automatic merges are explicit graph nodes here, never also submitted by producers.
   TASK_FILE=""; MANIFEST=""; OUTPUT_DIR=""; CONFIG_LOG_DIR=""; submit_array "$s" "$dep";dep="$LAST_JOB_ID"
 done
}
if [[ "$ARRAY_TASK_MODE" == "1" ]]; then
 [[ -n "$TASK_FILE" ]] || { echo "ERROR: --array-task requires --task-file" >&2; exit 2; }
 resolve_array_item
fi
if [[ "$SUBMIT_TO_SLURM" == 1 ]]; then
 [[ -z "${SLURM_JOB_ID:-}" ]] || { echo 'ERROR: cannot submit from a Slurm job' >&2;exit 2; }
 if [[ "$STEP" == all ]];then submit_workflow
 else submit_array "$STEP"; producer="$LAST_JOB_ID"
   if [[ "${AUTO_SUBMIT_MERGE:-true}" == true ]];then case "$STEP" in codon_match) TASK_FILE=""; submit_array codon_match_merge "$producer";echo "Submitted producer=$producer merge=$LAST_JOB_ID";; trna_match) TASK_FILE=""; submit_array trna_match_merge "$producer";echo "Submitted producer=$producer merge=$LAST_JOB_ID";; rrna_match) TASK_FILE=""; submit_array rrna_match_merge "$producer";echo "Submitted producer=$producer merge=$LAST_JOB_ID";; mitos2_annotation) TASK_FILE=""; submit_array mitos2_merge "$producer";echo "Submitted producer=$producer merge=$LAST_JOB_ID";; esac;fi
 fi
 exit 0
fi

# Keep the workflow interpreter stable: activating MITOS2 must not change the
# interpreter used by subsequent preprocessing steps when running `all`.
BASE_PYTHON="${PYTHON:-python3}"
COLLECT_SCRIPT="qc_analysis/scripts/collect_variant_calling_results.py"
LIFTOVER_SCRIPT="qc_analysis/scripts/run_coordinate_liftover.py"
CODON_SCRIPT="qc_analysis/scripts/run_codon_match.py"
CODON_TABLE_SCRIPT="qc_analysis/scripts/build_primate_codon_table.py"
MITOS2_SCRIPT="qc_analysis/scripts/run_mitos2_annotation.py"
COMPARISON_SCRIPT="qc_analysis/scripts/compare_genbank_mitos2_reference_annotations.py"
TRNA_SCRIPT="qc_analysis/scripts/run_trna_match.py"
TRNA_MERGE_SCRIPT="qc_analysis/scripts/merge_trna_match_summaries.py"
TRNA_INDEX_SCRIPT="qc_analysis/scripts/build_all_trna_indexes.py"
RRNA_SCRIPT="qc_analysis/scripts/run_rrna_match.py"
RRNA_MERGE_SCRIPT="qc_analysis/scripts/merge_rrna_match_summaries.py"
INTRASPECIES_SCRIPT="qc_analysis/scripts/run_intraspecies_contamination.py"
FINAL_FILTER_SCRIPT="qc_analysis/scripts/run_final_filter.py"
SAMPLE_FILTER_SCRIPT="qc_analysis/scripts/run_sample_variant_filtering.py"
GLOBAL_ANCHOR_SCRIPT="qc_analysis/scripts/discover_global_liftover_anchor.py"
HUMAN_CONTAMINATION_SCRIPT="qc_analysis/scripts/run_human_contamination.py"
PRIMATE_BACKGROUND_SCRIPT="qc_analysis/scripts/build_primate_homo_background.py"

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

# Read the MITOS2 conda settings without requiring a Python YAML parser before
# the environment that provides Biopython has been activated.
configured_mitos2_value() {
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

      if (content == "mitos2_annotation:") { mitos2_indent = level; in_mitos2 = 1; in_settings = 0; next }
      if (in_mitos2 && level <= mitos2_indent) { in_mitos2 = 0; in_settings = 0 }
      if (in_mitos2 && content == "settings:") { settings_indent = level; in_settings = 1; next }
      if (in_settings && level <= settings_indent) in_settings = 0
      if (in_settings && content ~ ("^" requested_key ":[[:space:]]*")) {
        sub("^" requested_key ":[[:space:]]*", "", content)
        print trim(content)
        exit
      }
    }
  ' "$CONFIG"
}

# Read the independent tRNAscan-SE conda settings without requiring PyYAML.
# This intentionally does not consult the MITOS2 environment configuration.
configured_trnascan_value() {
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

      if (content == "trna_match:") { trna_indent = level; in_trna = 1; in_settings = 0; next }
      if (in_trna && level <= trna_indent) { in_trna = 0; in_settings = 0 }
      if (in_trna && content == "settings:") { settings_indent = level; in_settings = 1; next }
      if (in_settings && level <= settings_indent) in_settings = 0
      if (in_settings && content ~ ("^" requested_key ":[[:space:]]*")) {
        sub("^" requested_key ":[[:space:]]*", "", content)
        print trim(content)
        exit
      }
    }
  ' "$CONFIG"
}

if [[ -z "${BIOPYTHON_USE_MODULE+x}" ]]; then
  configured_use_module="$(configured_biopython_value use_module)"
  case "$configured_use_module" in
    0|false|False|FALSE|no|No|NO) BIOPYTHON_USE_MODULE=0 ;;
    *) BIOPYTHON_USE_MODULE=1 ;;
  esac
fi
BIOPYTHON_MODULE="${BIOPYTHON_MODULE:-$(configured_biopython_value module_load)}"
BIOPYTHON_MODULE="${BIOPYTHON_MODULE:-Biopython/1.83-foss-2022b}"
MITOS2_CONDA_MODULE="$(configured_mitos2_value conda_module)"
MITOS2_CONDA_ENV="$(configured_mitos2_value conda_env)"
TRNASCAN_CONDA_MODULE="$(configured_trnascan_value conda_module)"
TRNASCAN_CONDA_ENV="$(configured_trnascan_value conda_env)"
TRNASCAN_BIN="$(configured_trnascan_value trnascan_bin)"

run_collect_variant_calling_results() {
  echo "[qc_preprocessing] Running collect_variant_calling_results with config: ${CONFIG}" >&2
  "$BASE_PYTHON" "$COLLECT_SCRIPT" --config "$CONFIG"
}

run_discover_global_anchor() {
  echo "[qc_preprocessing] Running discover_global_anchor with config: ${CONFIG}" >&2
  echo "[qc_preprocessing] MAFFT environment preflight:" >&2
  "$BASE_PYTHON" "$GLOBAL_ANCHOR_SCRIPT" --config "$CONFIG" --check-environment | while IFS= read -r line; do
    echo "[qc_preprocessing] ${line}" >&2
  done
  "$BASE_PYTHON" "$GLOBAL_ANCHOR_SCRIPT" --config "$CONFIG"
}

run_coordinate_liftover() {
  echo "[qc_preprocessing] Running coordinate_liftover with config: ${CONFIG}" >&2
  local cmd=("$BASE_PYTHON" "$LIFTOVER_SCRIPT" --config "$CONFIG")
  if [[ -n "${SAMPLE:-}" ]]; then
    cmd+=(--sample "$SAMPLE")
  fi
  "${cmd[@]}"
}

run_mitos2_annotation() {
  local mode="${1:-}"
  activate_mitos2_environment
  echo "[qc_preprocessing] Running mitos2_annotation ${mode} with config: ${CONFIG}" >&2
  local cmd=("$MITOS2_PYTHON" "$MITOS2_SCRIPT" --config "$CONFIG")
  [[ -n "$mode" ]] && cmd+=("$mode")
  if [[ "$STEP" == mitos2_annotation && "${ITEM:-}" == task:* ]]; then cmd+=(--task-id "${ITEM#task:}")
  elif [[ "$STEP" == mitos2_annotation && "${ITEM:-}" == reference:* ]]; then cmd+=(--reference "${ITEM#reference:}")
  elif [[ -n "${SAMPLE:-}" ]]; then cmd+=(--sample "$SAMPLE"); fi
  "${cmd[@]}"
}

activate_mitos2_environment() {
  if [[ -z "$MITOS2_CONDA_MODULE" || -z "$MITOS2_CONDA_ENV" ]]; then
    echo "ERROR: mitos2_annotation.settings must define conda_module and conda_env." >&2
    exit 1
  fi

  echo "[qc_preprocessing] Loading MITOS2 conda module: ${MITOS2_CONDA_MODULE}" >&2
  if ! command -v module >/dev/null 2>&1; then
    if [[ -f /etc/profile.d/modules.sh ]]; then
      # module is commonly initialized only for login shells on HPC systems.
      source /etc/profile.d/modules.sh
    fi
  fi
  if ! command -v module >/dev/null 2>&1; then
    echo "ERROR: module command is unavailable; cannot load ${MITOS2_CONDA_MODULE}." >&2
    exit 1
  fi
  module load "$MITOS2_CONDA_MODULE"

  # shellcheck disable=SC1090
  source "$(conda info --base)/etc/profile.d/conda.sh"
  echo "[qc_preprocessing] Activating MITOS2 conda environment: ${MITOS2_CONDA_ENV}" >&2
  conda activate "$MITOS2_CONDA_ENV"
  hash -r

  if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "ERROR: conda activation did not set CONDA_PREFIX for MITOS2 environment: ${MITOS2_CONDA_ENV}." >&2
    exit 1
  fi
  MITOS2_PYTHON="${CONDA_PREFIX}/bin/python"
  if [[ ! -x "$MITOS2_PYTHON" ]]; then
    echo "ERROR: MITOS2 Python is missing or not executable: ${MITOS2_PYTHON}" >&2
    exit 1
  fi

  echo "[qc_preprocessing] CONDA_DEFAULT_ENV=${CONDA_DEFAULT_ENV:-}" >&2
  echo "[qc_preprocessing] CONDA_PREFIX=${CONDA_PREFIX}" >&2
  echo "[qc_preprocessing] command -v python=$(command -v python || true)" >&2
  echo "[qc_preprocessing] MITOS2_PYTHON=${MITOS2_PYTHON}" >&2
  echo "[qc_preprocessing] MITOS2_PYTHON version=$($MITOS2_PYTHON --version 2>&1)" >&2

  if ! "$MITOS2_PYTHON" -c \
      'import sys, Bio; from Bio import SeqIO; print(sys.executable); print(Bio.__version__)'
  then
    echo "ERROR: Biopython is not importable in the MITOS2 conda environment: ${MITOS2_CONDA_ENV}." >&2
    exit 1
  fi
  echo "[qc_preprocessing] Biopython version=$($MITOS2_PYTHON -c 'import Bio; print(Bio.__version__)')" >&2
  echo "[qc_preprocessing] command -v runmitos=$(command -v runmitos || true)" >&2
}

activate_trnascan_environment() {
  if [[ -z "$TRNASCAN_CONDA_MODULE" || -z "$TRNASCAN_CONDA_ENV" ]]; then
    echo "ERROR: trna_match.settings must define conda_module and conda_env." >&2
    exit 1
  fi
  if [[ -z "$TRNASCAN_BIN" ]]; then
    echo "ERROR: trna_match.settings must define trnascan_bin." >&2
    exit 1
  fi

  echo "[qc_preprocessing] Loading tRNAscan-SE conda module: ${TRNASCAN_CONDA_MODULE}" >&2
  if ! command -v module >/dev/null 2>&1 && [[ -f /etc/profile.d/modules.sh ]]; then
    # module is commonly initialized only for login shells on HPC systems.
    source /etc/profile.d/modules.sh
  fi
  if ! command -v module >/dev/null 2>&1; then
    echo "ERROR: module command is unavailable; cannot load ${TRNASCAN_CONDA_MODULE}." >&2
    echo "ERROR: failed to activate tRNAscan environment." >&2
    exit 1
  fi
  if ! module load "$TRNASCAN_CONDA_MODULE"; then
    echo "ERROR: failed to load tRNAscan conda module: ${TRNASCAN_CONDA_MODULE}." >&2
    echo "ERROR: failed to activate tRNAscan environment." >&2
    exit 1
  fi

  local conda_base
  if ! conda_base="$(conda info --base)" || [[ -z "$conda_base" || ! -f "$conda_base/etc/profile.d/conda.sh" ]]; then
    echo "ERROR: could not locate conda initialization script for tRNAscan-SE." >&2
    echo "ERROR: failed to activate tRNAscan environment." >&2
    exit 1
  fi
  # shellcheck disable=SC1090
  source "$conda_base/etc/profile.d/conda.sh"
  echo "[qc_preprocessing] Activating tRNAscan conda environment: ${TRNASCAN_CONDA_ENV}" >&2
  if ! conda activate "$TRNASCAN_CONDA_ENV"; then
    echo "ERROR: failed to activate tRNAscan environment: ${TRNASCAN_CONDA_ENV}." >&2
    exit 1
  fi
  hash -r

  if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "ERROR: conda activation did not set CONDA_PREFIX for tRNAscan environment: ${TRNASCAN_CONDA_ENV}." >&2
    exit 1
  fi
  local trnascan_executable trnascan_version
  if ! trnascan_executable="$(command -v "$TRNASCAN_BIN")"; then
    echo "ERROR: tRNAscan-SE was not found after activating conda environment: ${TRNASCAN_CONDA_ENV}" >&2
    echo "Configured executable: ${TRNASCAN_BIN}" >&2
    exit 127
  fi
  # tRNAscan-SE releases differ in their handling of --version: some print
  # useful version or help text but still return a non-zero status.  Finding
  # the configured executable is the mandatory preflight; version output is
  # diagnostic only.
  trnascan_version="$("$TRNASCAN_BIN" --version 2>&1 || true)"
  if [[ -z "$trnascan_version" ]]; then
    "$TRNASCAN_BIN" -h >/dev/null 2>&1 || true
  fi

  echo "[qc_preprocessing] CONDA_DEFAULT_ENV=${CONDA_DEFAULT_ENV:-}" >&2
  echo "[qc_preprocessing] CONDA_PREFIX=${CONDA_PREFIX}" >&2
  echo "[qc_preprocessing] tRNAscan-SE executable=${trnascan_executable}" >&2
  echo "[qc_preprocessing] tRNAscan-SE version/info=${trnascan_version}" >&2
}

run_build_trna_indexes() {
  activate_trnascan_environment
  echo "[qc_preprocessing] Running build_trna_indexes with config: ${CONFIG}" >&2
  local workers
  workers="$(trna_setting index_build_workers)"; workers="${workers:-1}"
  [[ "$workers" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: index_build_workers must be >= 1: $workers" >&2; exit 2; }
  "$BASE_PYTHON" "$TRNA_INDEX_SCRIPT" --config "$CONFIG" --workers "$workers"
}

run_trna_match() {
  activate_trnascan_environment
  run_annotation trna_match "$TRNA_SCRIPT"
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

  if ! "$BASE_PYTHON" - <<'PY'
from Bio import Entrez, SeqIO
print("Biopython import OK")
PY
  then
    echo "ERROR: Biopython is not importable after loading the configured module." >&2
    echo "Tried module: ${BIOPYTHON_MODULE}" >&2
    echo "Please check the HPC module name or set BIOPYTHON_USE_MODULE=0 if using a Python environment that already has Biopython." >&2
    exit 1
  fi

  local cmd=("$BASE_PYTHON" "$CODON_TABLE_SCRIPT" --config "$CONFIG")
  [[ -n "${SAMPLE:-}" ]] && cmd+=(--sample "$SAMPLE")
  local workers="${CODON_TABLE_WORKERS:-${SLURM_CPUS_PER_TASK:-1}}"
  if [[ "$workers" =~ ^[1-9][0-9]*$ ]]; then
    cmd+=(--workers "$workers")
  else
    echo "WARNING: ignoring invalid CODON_TABLE_WORKERS value: ${workers}" >&2
  fi
  if [[ -n "${CODON_TABLE_DOWNLOAD_WORKERS:-}" ]]; then
    if [[ "${CODON_TABLE_DOWNLOAD_WORKERS}" =~ ^[1-9][0-9]*$ ]]; then
      cmd+=(--download-workers "${CODON_TABLE_DOWNLOAD_WORKERS}")
    else
      echo "WARNING: ignoring invalid CODON_TABLE_DOWNLOAD_WORKERS value: ${CODON_TABLE_DOWNLOAD_WORKERS}" >&2
    fi
  fi
  "${cmd[@]}"
}

run_compare_genbank_mitos2() {
  echo "[qc_preprocessing] Running compare_genbank_mitos2 with config: ${CONFIG}" >&2
  "$BASE_PYTHON" "$COMPARISON_SCRIPT" --config "$CONFIG"
}

comparison_enabled() {
  awk '
    /^[^[:space:]]/ { in_section = ($0 ~ /^genbank_mitos2_comparison:/) }
    in_section && /^[[:space:]]+enabled:[[:space:]]*/ { value=$0; sub(/.*enabled:[[:space:]]*/, "", value); print value; exit }
  ' "$CONFIG" | grep -Eiq '^(true|yes|1)$'
}

run_annotation() {
  local name="$1" script="$2"
  echo "[qc_preprocessing] Running ${name} with config: ${CONFIG}" >&2
  local cmd=("$BASE_PYTHON" "$script" --config "$CONFIG")
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
  human_contamination) "$BASE_PYTHON" "$HUMAN_CONTAMINATION_SCRIPT" --config "$CONFIG" ;;
  build_primate_homo_background) "$BASE_PYTHON" "$PRIMATE_BACKGROUND_SCRIPT" --config "$CONFIG" ;;
  mitos2_prepare_tasks) run_mitos2_annotation --prepare-tasks ;;
  mitos2_merge) run_mitos2_annotation --merge-only ;;
  mitos2_annotation) run_mitos2_annotation ;;
  build_primate_codon_table) run_build_primate_codon_table ;;
  compare_genbank_mitos2) run_compare_genbank_mitos2 ;;
  codon_match) run_annotation codon_match "$CODON_SCRIPT" ;;
  codon_match_validate) "$BASE_PYTHON" "$CODON_SCRIPT" --config "$CONFIG" --validate-inputs ;;
  codon_match_merge) "$BASE_PYTHON" "$CODON_SCRIPT" --config "$CONFIG" --merge-summaries ;;
  build_trna_indexes) run_build_trna_indexes ;;
  trna_match) run_trna_match ;;
  trna_match_merge) "$BASE_PYTHON" "$TRNA_MERGE_SCRIPT" --config "$CONFIG" ;;
  trna_gene_qc) echo 'Run run_trna_gene_liftover_qc.py with source index, human index, and coordinate map for each sample.' ;;
  rrna_match) run_annotation rrna_match "$RRNA_SCRIPT" ;;
  rrna_match_merge) "$BASE_PYTHON" "$RRNA_MERGE_SCRIPT" --config "$CONFIG" ;;
  intraspecies_contamination) "$BASE_PYTHON" "$INTRASPECIES_SCRIPT" --config "$CONFIG" ;;
  sample_variant_filtering) "$BASE_PYTHON" "$SAMPLE_FILTER_SCRIPT" --config "$CONFIG" ;;
  final_filter) "$BASE_PYTHON" "$FINAL_FILTER_SCRIPT" --config "$CONFIG" ;;
  all)
    run_collect_variant_calling_results
    "$BASE_PYTHON" "$INTRASPECIES_SCRIPT" --config "$CONFIG"
    "$BASE_PYTHON" "$SAMPLE_FILTER_SCRIPT" --config "$CONFIG"
    run_discover_global_anchor
    run_coordinate_liftover
    run_mitos2_annotation
    "$BASE_PYTHON" "$CODON_SCRIPT" --config "$CONFIG" --validate-inputs
    run_annotation codon_match "$CODON_SCRIPT"
    "$BASE_PYTHON" "$CODON_SCRIPT" --config "$CONFIG" --merge-summaries
    run_build_trna_indexes
    run_trna_match
    "$BASE_PYTHON" "$TRNA_MERGE_SCRIPT" --config "$CONFIG"
    run_annotation rrna_match "$RRNA_SCRIPT"
    "$BASE_PYTHON" "$RRNA_MERGE_SCRIPT" --config "$CONFIG"
    "$BASE_PYTHON" "$PRIMATE_BACKGROUND_SCRIPT" --config "$CONFIG"
    "$BASE_PYTHON" "$HUMAN_CONTAMINATION_SCRIPT" --config "$CONFIG"
    "$BASE_PYTHON" "$FINAL_FILTER_SCRIPT" --config "$CONFIG"
    ;;
esac
