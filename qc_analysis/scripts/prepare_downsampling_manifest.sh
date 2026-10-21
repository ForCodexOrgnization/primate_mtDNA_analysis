#!/usr/bin/env bash
set -euo pipefail
ROOT=/vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis_coverage_only
BASE="$ROOT/results/qc_coverage_only/downsampling_pilot"
PILOT="$BASE/pilot_samples.tsv"
MANIFEST="$BASE/downsampling_jobs.tsv"
[[ -s "$PILOT" ]] || { echo "Missing $PILOT" >&2; exit 1; }
awk -F '\t' -v OFS='\t' '
NR==1 {
  for(i=1;i<=NF;i++) h[$i]=i
  for(i=1;i<=3;i++) {
    col=(i==1?"sample":i==2?"species":"mt_median_coverage")
    if(!h[col]) { print "Missing column " col > "/dev/stderr"; exit 2 }
  }
  print "sample","species","original_median_est","target_depth","seed","fraction"
  next
}
{
  s=$(h["sample"]); sp=$(h["species"]); original=$(h["mt_median_coverage"])+0
  if(original<=0 || s=="") { print "Invalid pilot row: " NR > "/dev/stderr"; exit 2 }
  for(i=1;i<=3;i++) {
    target=(i==1?30:i==2?50:100)
    for(j=1;j<=2;j++) {
      seed=(j==1?101:202)
      printf "%s\t%s\t%.6f\t%d\t%d\t%.10f\n",s,sp,original,target,seed,target/original
    }
  }
}' "$PILOT" > "$MANIFEST"
echo "Prepared $(( $(wc -l < "$MANIFEST") - 1 )) tasks: $MANIFEST"
head -n 7 "$MANIFEST"
