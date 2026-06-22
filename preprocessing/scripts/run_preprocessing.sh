#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash preprocessing/scripts/run_preprocessing.sh <step> config/preprocessing_paths.yaml

Steps:
  reference_discovery        Run reference discovery only.
  reference_materialization  Build/materialize reviewed reference manifests.
  in_house_score             Launch/iterate in-house score jobs.
  variant_inputs             Prepare variant-calling input table.
  post_reference_review      Run reference_materialization, in_house_score, variant_inputs.
  all_steps                  Run all preprocessing stages using the raw discovery summary without manual review.
  reports                    Render preprocessing Quarto reports.
  all                        Run reference_discovery only, then stop for manual review.

Important:
  The all step intentionally stops after reference discovery because
  species_reference_chrM_summary.tsv must be manually reviewed before downstream steps.
USAGE
}

if [[ $# -lt 1 || $# -gt 2 ]]; then
  usage >&2
  exit 2
fi

STEP="$1"
CONFIG="${2:-config/preprocessing_paths.yaml}"

if [[ ! -s "$CONFIG" ]]; then
  echo "ERROR: missing or empty config file: $CONFIG" >&2
  exit 1
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

SPECIES_REFERENCE_SUMMARY=$(config_get species_reference_summary "data/metadata/species_reference_chrM_summary.tsv")
REFERENCE_MATERIALIZATION_RESULTS_MANIFEST=$(config_get reference_materialization_results_manifest "results/preprocessing/reference_materialization/reference_materialization_manifest.tsv")
REFERENCE_MATERIALIZATION_MANIFEST=$(config_get reference_materialization_manifest "references/manifests/reference_materialization_manifest.tsv")
IN_HOUSE_SCORE_REFERENCE_INPUTS=$(config_get in_house_score_reference_inputs "references/manifests/in_house_score_reference_inputs.tsv")
SAMPLE_METADATA=$(config_get sample_metadata "data/metadata/sample_metadata.tsv")
MERGED_IN_HOUSE_SCORE=$(config_get merged_in_house_score "results/preprocessing/in_house_score/merged_in_house_score.tsv")
VARIANT_CALLING_INPUT_TABLE=$(config_get variant_calling_input_table "results/preprocessing/variant_calling_inputs/variant_calling_input_table.tsv")

run_reference_discovery() {
  echo "[preprocessing] Running reference discovery with species table: $SPECIES_TABLE" >&2
  SPECIES_TABLE="$SPECIES_TABLE" \
  MITO_FASTA="$MITO_FASTA" \
  TREE_NEWICK="$TREE_NEWICK" \
  OUTDIR="$REFERENCE_DISCOVERY_OUTDIR" \
  EMAIL="$EMAIL" \
  MAX_NEAREST="$MAX_NEAREST" \
  DELAY="$DELAY" \
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
  Rscript preprocessing/scripts/build_reference_materialization_manifest.R \
    "$SPECIES_REFERENCE_SUMMARY" \
    "$REFERENCE_MATERIALIZATION_RESULTS_MANIFEST" \
    "$REFERENCE_MATERIALIZATION_MANIFEST"
  echo "[preprocessing] Materializing references from: $REFERENCE_MATERIALIZATION_MANIFEST" >&2
  bash preprocessing/scripts/materialize_references.sh "$REFERENCE_MATERIALIZATION_MANIFEST"
}

run_in_house_score() {
  echo "[preprocessing] Running in-house score step with: $IN_HOUSE_SCORE_REFERENCE_INPUTS" >&2
  REF_INPUTS="$IN_HOUSE_SCORE_REFERENCE_INPUTS" \
    bash preprocessing/scripts/run_in_house_score_array.sh
}

run_variant_inputs() {
  echo "[preprocessing] Preparing variant-calling inputs: $VARIANT_CALLING_INPUT_TABLE" >&2
  Rscript preprocessing/scripts/prepare_variant_calling_inputs.R \
    "$SAMPLE_METADATA" \
    "$IN_HOUSE_SCORE_REFERENCE_INPUTS" \
    "$MERGED_IN_HOUSE_SCORE" \
    "$VARIANT_CALLING_INPUT_TABLE"
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
  variant_inputs)
    run_variant_inputs
    ;;
  post_reference_review)
    run_reference_materialization
    run_in_house_score
    run_variant_inputs
    ;;
  all_steps)
    run_reference_discovery
    use_raw_discovery_summary_without_review
    run_reference_materialization
    run_in_house_score
    run_variant_inputs
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
