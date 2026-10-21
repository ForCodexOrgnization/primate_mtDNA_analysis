#!/usr/bin/env bash
set -euo pipefail
# Copy CRAM/CRAI to large scratch, verify full-file checksums, keep originals.
ROOT=/vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis_coverage_only
BASE="$ROOT/results/qc_coverage_only/downsampling_pilot"
SRC=/home/lt692/ycga_work/primate_downsampling_pilot/wgs_input
DST="$BASE/wgs_input"
SAMPLES="$BASE/pilot_sample_ids.txt"
REPORT="$BASE/wgs_staging_report.tsv"
[[ -s "$SAMPLES" ]] || { echo "Missing sample list: $SAMPLES" >&2; exit 1; }
command -v rsync >/dev/null || { echo 'rsync required' >&2; exit 1; }
mkdir -p "$DST"
printf 'sample\tcram_bytes\tcrai_bytes\tchecksum_verified\n' > "$REPORT"
while IFS= read -r sample || [[ -n "$sample" ]]; do
  sample=${sample%$'\r'}
  [[ -z "$sample" || "$sample" == \#* ]] && continue
  [[ "$sample" =~ ^[A-Za-z0-9_.-]+$ ]] || { echo "Invalid ID: $sample" >&2; exit 1; }
  from="$SRC/$sample"; to="$DST/$sample"
  [[ -s "$from/$sample.cram" && -s "$from/$sample.cram.crai" ]] || {
    echo "Missing input CRAM/CRAI for $sample" >&2; exit 1;
  }
  mkdir -p "$to"
  rsync -a --partial -- "$from/$sample.cram" "$from/$sample.cram.crai" "$to/"
  diff_list=$(rsync -acn --out-format='%n' -- "$from/$sample.cram" "$from/$sample.cram.crai" "$to/")
  [[ -z "$diff_list" ]] || { echo "Checksum mismatch for $sample: $diff_list" >&2; exit 1; }
  printf '%s\t%s\t%s\tYES\n' "$sample" \
    "$(stat -c %s "$to/$sample.cram")" "$(stat -c %s "$to/$sample.cram.crai")" >> "$REPORT"
  echo "STAGED: $sample"
done < "$SAMPLES"
echo "DONE: $REPORT (source CRAMs retained)"