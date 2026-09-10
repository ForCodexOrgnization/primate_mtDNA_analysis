#!/usr/bin/env bash
# Run QC from a canonical biological-sample cohort.
#
# This is the preferred entrypoint when enriched sample metadata is available.
# It materializes one representative accession per biological sample before
# collection, then all downstream scripts consume the generated sample_ref file
# through qc_analysis.lib.simple_yaml.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"

PYTHON="${PYTHON:-python3}"
DEDUP_CONFIG="${DEDUP_CONFIG:-config/sample_deduplication.yaml}"
QC_CONFIG="${QC_CONFIG:-config/qc_preprocessing.yaml}"
DEDUP_SCRIPT="qc_analysis/scripts/run_sample_deduplication.py"
COLLECT_SCRIPT="qc_analysis/scripts/collect_variant_calling_results.py"
BASE_WRAPPER="qc_analysis/scripts/run_qc_preprocessing.sh"
NUMT_PROPAGATION_SCRIPT="qc_analysis/scripts/propagate_species_numt_positions.py"
DEDUP_REF="results/qc/sample_deduplication/reports/deduplicated_sample_ref_file.tsv"

usage() {
  cat <<'USAGE'
Usage:
  bash qc_analysis/scripts/run_qc_preprocessing_deduplicated.sh <step> [config/qc_preprocessing.yaml]

The enriched metadata location is configured in config/sample_deduplication.yaml.
The deduplication preflight always runs first.

Examples:
  bash qc_analysis/scripts/run_qc_preprocessing_deduplicated.sh collect_variant_calling_results
  bash qc_analysis/scripts/run_qc_preprocessing_deduplicated.sh sample_variant_filtering
  bash qc_analysis/scripts/run_qc_preprocessing_deduplicated.sh local_heteroplasmy_qc
  bash qc_analysis/scripts/run_qc_preprocessing_deduplicated.sh all
USAGE
}

[[ $# -ge 1 && $# -le 2 ]] || { usage >&2; exit 2; }
STEP="$1"
QC_CONFIG="${2:-$QC_CONFIG}"

run_dedup() {
  echo "[qc_deduplicated] Building canonical sample cohort from ${DEDUP_CONFIG}" >&2
  "$PYTHON" "$DEDUP_SCRIPT" --config "$DEDUP_CONFIG"
  [[ -s "$DEDUP_REF" ]] || { echo "ERROR: deduplicated sample_ref was not created: $DEDUP_REF" >&2; exit 1; }
}

run_collect_dedup() {
  echo "[qc_deduplicated] Collecting only canonical/unique accessions from ${DEDUP_REF}" >&2
  "$PYTHON" "$COLLECT_SCRIPT" \
    --config "$QC_CONFIG" \
    --metadata "$DEDUP_REF" \
    --metadata-sample-column sample \
    --metadata-species-column species
}

run_step() {
  local step="$1"
  bash "$BASE_WRAPPER" "$step" "$QC_CONFIG"
  if [[ "$step" == "local_heteroplasmy_qc" ]]; then
    echo "[qc_deduplicated] Propagating NUMT-supported positions across same-species samples" >&2
    "$PYTHON" "$NUMT_PROPAGATION_SCRIPT" --config "$QC_CONFIG"
  fi
}

run_dedup

case "$STEP" in
  collect_variant_calling_results)
    run_collect_dedup
    ;;
  all)
    run_collect_dedup
    run_step sample_variant_filtering
    run_step pre_liftover_variant_qc
    run_step intraspecies_contamination
    run_step local_heteroplasmy_qc
    run_step discover_global_anchor
    run_step coordinate_liftover
    run_step interspecies_contamination
    run_step mitos2_annotation
    run_step codon_match_validate
    run_step codon_match
    run_step codon_match_merge
    run_step build_trna_indexes
    run_step trna_match
    run_step trna_match_merge
    run_step rrna_match
    run_step rrna_match_merge
    run_step build_primate_homo_background
    run_step final_filter
    ;;
  *)
    # The shared YAML reader automatically replaces legacy
    # config/sample_ref_file.tsv paths with DEDUP_REF after run_dedup.
    run_step "$STEP"
    ;;
esac
