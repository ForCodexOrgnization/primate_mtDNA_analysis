#!/usr/bin/env bash
# Initialize an isolated coverage-only QC workspace in a NEW Git clone.
# No Slurm submission. Never modifies source or existing destination.
set -Eeuo pipefail
if [[ $# -ne 1 ]]; then
  echo "Usage: bash qc_analysis/scripts/prepare_coverage_only.sh /path/to/primate_mtDNA_analysis_incremental" >&2
  exit 2
fi
cd "$(dirname "$0")/../.."
target_root="$(pwd -P)"
source_root="$(readlink -f "$1")"
target="$target_root/results/qc_coverage_only"
source="$source_root/results/qc_incremental"
cfg="config/qc_preprocessing_coverage_only.yaml"
[[ -s "$cfg" && -d "$source" ]] || { echo "ERROR: missing config or source QC directory" >&2; exit 1; }
[[ "$source_root" != "$target_root" ]] || { echo "ERROR: use a separate Git clone for coverage-only" >&2; exit 1; }
[[ ! -e "$target" && ! -L "$target" ]] || { echo "ERROR: destination already exists: $target" >&2; exit 1; }
command -v rsync >/dev/null || { echo "ERROR: rsync required" >&2; exit 127; }
if [[ ! -e references && -d "$source_root/references" ]]; then
  ln -s "$source_root/references" references
fi
if [[ ! -d data/reference_tables && -d "$source_root/data/reference_tables" ]]; then
  mkdir -p data/reference_tables
  rsync -a "$source_root/data/reference_tables/" data/reference_tables/
fi
if [[ ! -s config/sample_ref_file.tsv && -s "$source_root/config/sample_ref_file.tsv" ]]; then
  cp "$source_root/config/sample_ref_file.tsv" config/sample_ref_file.tsv
fi
mkdir -p "$target"
echo "[coverage_only] Copying per-sample results and source metadata from $source"
rsync -a --exclude='job_arrays/' --exclude='logs/' --exclude='final_filter/' \
  --exclude='local_heteroplasmy_qc/' --exclude='global_anchor/' \
  --exclude='primate_homo_background/' --exclude='intraspecies_contamination/' \
  --exclude='interspecies_contamination/' --exclude='human_contamination/' \
  "$source/" "$target/"
# Remove stale aggregated reports only in new workspace; keep per-sample annotations.
for p in \
  "$target/pre_liftover_variant_qc/reports/source_variant_qc.tsv" \
  "$target/pre_liftover_variant_qc/reports/source_variant_qc_summary.tsv" \
  "$target/codon_match/reports/codon_match_summary.tsv" \
  "$target/trna_match/reports/trna_match_summary.tsv" \
  "$target/rrna_match/reports/rrna_match_summary.tsv"; do
  [[ ! -f "$p" ]] || mv "$p" "$p.precoverage_snapshot"
done
# Prefer the latest current collection inventory; never rewrite source collection.
summary="$source/collected_variant_calling_results/reports/variant_calling_collection_summary.tsv"
[[ -s "$summary" ]] || { echo "ERROR: missing current source collection summary: $summary" >&2; exit 1; }
cp "$summary" "$target/updated_collection_summary.tsv"
python3 - <<'PY'
import csv
from pathlib import Path
p=Path("results/qc_coverage_only/updated_collection_summary.tsv")
with p.open(newline="") as f:
    r=csv.DictReader(f,delimiter=chr(9))
    assert r.fieldnames and {"sample","status","Percent_100"} <= set(r.fieldnames)
    rows=list(r)
assert len(rows)==len({x["sample"] for x in rows}), "duplicate sample"
def valid(v):
    try: return float(v)>=90
    except (TypeError,ValueError): return False
eligible=[x for x in rows if x["status"]=="OK" and valid(x["Percent_100"])]
print(f"[coverage_only] inventory={len(rows)}; status_OK_and_Percent_100_ge90={len(eligible)}")
PY
echo "[coverage_only] Independent copy prepared at $target"
echo "[coverage_only] Source untouched: $source"
echo "[coverage_only] NEXT: run sample_variant_filtering with $cfg, inspect counts, then dry-run."
