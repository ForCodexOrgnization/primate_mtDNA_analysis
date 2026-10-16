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
  sample_variant_filtering         Write the basic three-criterion biological sample QC report.
  pre_liftover_variant_qc          Annotate native VCFs with immutable source-call QC.
  local_heteroplasmy_qc            Run native-coordinate NUMT/indel local-artifact QC and removal-list generation.
  discover_global_anchor           Discover reference-level global MSA anchors only.
  coordinate_liftover              Run coordinate liftover only.
  interspecies_contamination       Report cross-species contamination in lifted alleles.
  build_primate_homo_background   Consolidate orthology QC and build temporary homoplasmic background.
  human_contamination              Screen annotated human-coordinate alleles with primate background correction (standalone; temporarily excluded from all).
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
  final_filter                     Materialize final files and apply SOURCE-keyed NUMT/indel removals.
  all                              Run production preprocessing/annotation steps, currently excluding human_contamination.

Run modes:
  --submit                         Submit this wrapper to Slurm from a login/frontend node.
  --incremental                    Skip validated completed per-sample VCFs when submitting; rerun global reports.
  --start-from STEP                With --submit all, start at STEP and submit following steps.
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
  SLURM_PARTITION                  Partition for --submit (default: ycga; override with env var).
  SLURM_TIME                       Walltime for --submit (default: 24:00:00).
  SLURM_MEM                        Memory for --submit (default: 16G).
  SLURM_CPUS                       CPUs for --submit (default: 4).
  SLURM_LOG_DIR                    Override the step-specific scheduler log directory.
  SLURM_JOB_NAME                   Job name for --submit (default: qc_preprocessing_<step>).

Options: --incremental, --start-from STEP, --array-concurrency N, --task-file PATH, --prepare-only,
         --prepare-retry, --dry-run-submit. Direct mode remains sequential;
         --submit is recommended for HPC execution.

New array metadata is stored in <step output>/job_arrays and logs default to
<step output>/logs/job_arrays. Optional paths.job_array_dir and paths.log_dir
settings override those derived paths. Historical results/qc/job_arrays files
are left untouched; move them manually only after confirming no old job uses them.

Examples:
  bash qc_analysis/scripts/run_qc_preprocessing.sh --submit all config/qc_preprocessing.yaml
  bash qc_analysis/scripts/run_qc_preprocessing.sh --dry-run-submit --incremental --start-from discover_global_anchor all config/qc_preprocessing_incremental_strict.yaml
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

SUBMIT_TO_SLURM=0; DRY_RUN_SUBMIT=0; PREPARE_ONLY=0; PREPARE_RETRY=0; ARRAY_TASK_MODE=0; RUNTIME_ELIGIBILITY_MODE=0; INCREMENTAL_MODE=0; START_FROM=""
SAMPLE="${SAMPLE:-}"; TASK_FILE=""; ARRAY_CONCURRENCY="${SLURM_ARRAY_CONCURRENCY:-20}"
while [[ "${1:-}" == --* ]]; do
  case "$1" in
    --submit) SUBMIT_TO_SLURM=1; shift ;;
    --dry-run-submit) SUBMIT_TO_SLURM=1; DRY_RUN_SUBMIT=1; shift ;;
    --prepare-only) PREPARE_ONLY=1; shift ;;
    --prepare-retry) PREPARE_RETRY=1; shift ;;
    --incremental) INCREMENTAL_MODE=1; shift ;;
    --start-from)
      [[ $# -ge 2 && -n "$2" ]] || { echo 'ERROR: --start-from requires a workflow step' >&2; exit 2; }
      START_FROM="$2"; shift 2 ;;

    --runtime-eligibility) RUNTIME_ELIGIBILITY_MODE=1; shift ;;
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
 collect_variant_calling_results|pre_liftover_variant_qc|pre_liftover_variant_qc_merge|local_heteroplasmy_qc|discover_global_anchor|coordinate_liftover|interspecies_contamination|build_primate_homo_background|human_contamination|build_primate_codon_table|compare_genbank_mitos2|mitos2_prepare_tasks|mitos2_merge|mitos2_annotation|codon_match|codon_match_validate|codon_match_merge|build_trna_indexes|trna_match|trna_match_merge|trna_gene_qc|rrna_match|rrna_match_merge|intraspecies_contamination|sample_variant_filtering|final_filter|all) ;;
 -h|--help|help) usage; exit 0;; *) echo "ERROR: unknown step: $STEP" >&2; exit 2;; esac
[[ -s "$CONFIG" ]] || { echo "ERROR: missing or empty config file: $CONFIG" >&2; exit 1; }
if [[ -n "$START_FROM" && ( "$STEP" != all || "$SUBMIT_TO_SLURM" != 1 || "$ARRAY_TASK_MODE" == 1 ) ]]; then
  echo 'ERROR: --start-from is only supported with --submit all' >&2; exit 2
fi
export SAMPLE

classify_step() { case "$1" in pre_liftover_variant_qc|coordinate_liftover|codon_match|trna_match|rrna_match) echo sample;; mitos2_annotation) echo reference;; *) echo singleton;; esac; }
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
  if [[ "$(classify_step "$STEP")" == sample ]]; then
    SAMPLE="$ITEM"
  else
    SAMPLE=""
  fi
}
prepare_task_manifest() {
 local args=(python3 qc_analysis/scripts/qc_array_manifest.py "$1" "$CONFIG")
 [[ "$INCREMENTAL_MODE" == 1 ]] && args+=(--incremental)
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
 else prepare_task_manifest "$step"; fi
 [[ -n "${OUTPUT_DIR:-}" ]] || resolve_runtime_paths "$step"
 if (( TASK_COUNT == 0 )); then
   printf '[qc_preprocessing] step=%s state=complete_noop dependency_preserved=%s task_file=%s manifest=%s\n' "$step" "${dependency:-none}" "$TASK_FILE" "${MANIFEST:-unknown}" >&2
   LAST_JOB_ID="$dependency"
   return 0
 fi
 local kind array;kind=$(classify_step "$step");array=$(build_array_expression "$TASK_COUNT" "$kind");resolve_step_resources "$step"
 local logs="${SLURM_LOG_DIR:-$CONFIG_LOG_DIR}";mkdir -p "$logs"
 local cmd=(sbatch --parsable --job-name="qc_preprocessing_${step}" --array="$array" --output="$logs/%A_%a.out" --error="$logs/%A_%a.err" --time="$RES_TIME" --mem="$RES_MEM" --cpus-per-task="$RES_CPUS")
 cmd+=(--partition="${SLURM_PARTITION:-ycga}"); [[ -n "$dependency" ]] && cmd+=(--dependency="afterok:$dependency")
 cmd+=("$(readlink -f "${BASH_SOURCE[0]}")" --array-task --runtime-eligibility)
  [[ "$INCREMENTAL_MODE" == 1 ]] && cmd+=(--incremental)
  cmd+=(--task-file "$TASK_FILE" "$step" "$CONFIG")
 printf '[qc_preprocessing] step=%s\n[qc_preprocessing] output_dir=%s\n[qc_preprocessing] task_file=%s\n[qc_preprocessing] manifest=%s\n[qc_preprocessing] logs=%s/%%A_%%a.{out,err}\n[qc_preprocessing] task_count=%s\n[qc_preprocessing] array=%s\n' "$step" "$OUTPUT_DIR" "$TASK_FILE" "${MANIFEST:-unknown}" "$logs" "$TASK_COUNT" "$array" >&2
 printf '[qc_preprocessing] concurrency=%s dependency=%s resources=%s,%s,%s\n' "$ARRAY_CONCURRENCY" "${dependency:-none}" "$RES_TIME" "$RES_MEM" "$RES_CPUS" >&2
 if [[ "$PREPARE_ONLY" == 1 ]]; then LAST_JOB_ID="prepared_${step}"; return; fi
 if [[ "$DRY_RUN_SUBMIT" == 1 ]]; then printf 'DRY RUN:';printf ' %q' "${cmd[@]}";printf '\n'; LAST_JOB_ID="dry_${step}"; else command -v sbatch >/dev/null || { echo 'ERROR: --submit requires sbatch on PATH' >&2;exit 127; };LAST_JOB_ID=$("${cmd[@]}");LAST_JOB_ID=${LAST_JOB_ID%%;*};fi
 local submission="${TASK_FILE%.tasks.txt}.submission.tsv" submitted_at; submitted_at=$(date -u +%FT%TZ)
 printf 'step\tjob_id\ttask_file\tmanifest\tlog_dir\tarray\tsubmitted_at\n%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$step" "$LAST_JOB_ID" "$TASK_FILE" "${MANIFEST:-}" "$logs" "$array" "$submitted_at" >"$submission"
}
contamination_step_enabled() {
  local section="$1"
  awk -v section="$section" '
    /^[^[:space:]#][^:]*:/ { name=$0; sub(/:.*/, "", name); in_section=(name==section) }
    in_section && /^[[:space:]]+enabled:[[:space:]]*/ {
      value=$0; sub(/.*enabled:[[:space:]]*/, "", value);
      sub(/[[:space:]#].*$/, "", value); print tolower(value); exit
    }
  ' "$CONFIG" | grep -Eq '^(true|yes|1)$'
}
