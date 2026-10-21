#!/usr/bin/env bash
#SBATCH --job-name=primate_hom_ds
#SBATCH --partition=day
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --output=results/qc_coverage_only/downsampling_pilot/logs/hom_ds_%A_%a.out
#SBATCH --error=results/qc_coverage_only/downsampling_pilot/logs/hom_ds_%A_%a.err
set -euo pipefail
module load SAMtools/1.21-GCC-13.3.0
ROOT=/vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis_coverage_only
BASE="$ROOT/results/qc_coverage_only/downsampling_pilot"
MANIFEST="$BASE/downsampling_jobs.tsv"
THREADS=${SLURM_CPUS_PER_TASK:-4}
TASK=${SLURM_ARRAY_TASK_ID:-1}

[[ -s "$MANIFEST" ]] || { echo "Missing manifest $MANIFEST" >&2; exit 1; }
line=$(sed -n "$((TASK+1))p" "$MANIFEST")
[[ -n "$line" ]] || { echo "Task $TASK not found" >&2; exit 1; }
IFS=$'\t' read -r SAMPLE SPECIES ORIGINAL_EST TARGET SEED FRACTION <<< "$line"
[[ "$FRACTION" =~ ^0?\.[0-9]+$|^1(\.0+)?$ ]] || {
  echo "Invalid fraction $FRACTION" >&2; exit 1;
}
INPUT_ROOT=${WGS_INPUT_ROOT:-/home/lt692/ycga_work/primate_downsampling_pilot/wgs_input}
IN="$INPUT_ROOT/$SAMPLE/$SAMPLE.cram"
[[ -s "$IN" && -s "$IN.crai" ]] || { echo "Missing CRAM/CRAI: $IN" >&2; exit 1; }
echo "Input CRAM: $IN (read-only); BAM output on scratch"

REF_ARGS=()
if [[ -n ${CRAM_REFERENCE_FASTA:-} ]]; then
  [[ -s "$CRAM_REFERENCE_FASTA" ]] || { echo "No reference: $CRAM_REFERENCE_FASTA" >&2; exit 1; }
  REF_ARGS=(-T "$CRAM_REFERENCE_FASTA")
fi

# Explicit MT_CONTIG wins. Otherwise require an unambiguous short contig with mapped reads.
if [[ -z ${MT_CONTIG:-} ]]; then
  mapfile -t CANDIDATES < <(samtools idxstats "$IN" | awk -F '\t' '$2>=14000 && $2<=20000 && $3>0 {print $1}')
  if (( ${#CANDIDATES[@]} != 1 )); then
    echo "Ambiguous mt contig (${#CANDIDATES[@]} candidates): ${CANDIDATES[*]:-none}" >&2
    echo "Check the sample reference manifest and set MT_CONTIG explicitly." >&2
    exit 1
  fi
  MT_CONTIG="${CANDIDATES[0]}"
fi

OUTDIR="$BASE/downsampled/$SAMPLE"
mkdir -p "$OUTDIR"
PREFIX="$SAMPLE.ds${TARGET}x.seed${SEED}"
OUT="$OUTDIR/$PREFIX.bam"
METRICS="$OUTDIR/$PREFIX.metrics.tsv"
echo "sample=$SAMPLE species=$SPECIES target=$TARGET seed=$SEED frac=$FRACTION mt=$MT_CONTIG"

samtools quickcheck -v "$IN"
if [[ ! -s "$OUT" ]]; then
  TMP_BAM="$OUT.tmp.$$.bam"
  trap 'rm -f "$TMP_BAM"' EXIT
  samtools view -@ "$THREADS" "${REF_ARGS[@]}" -b \
    --subsample-seed "$SEED" --subsample "$FRACTION" \
    -o "$TMP_BAM" "$IN"
  samtools quickcheck -v "$TMP_BAM"
  mv "$TMP_BAM" "$OUT"
fi
[[ -s "$OUT.bai" ]] || samtools index -@ "$THREADS" "$OUT"
samtools quickcheck -v "$OUT"

# Include zero-depth sites in the denominator.
TMP_DEPTH="$OUTDIR/$PREFIX.depth.$$.tsv"
trap 'rm -f "$TMP_DEPTH"' EXIT
samtools depth -aa -r "$MT_CONTIG" "$OUT" > "$TMP_DEPTH"
[[ -s "$TMP_DEPTH" ]] || { echo "No mt depth output for $MT_CONTIG" >&2; exit 1; }

MEDIAN=$(cut -f3 "$TMP_DEPTH" | sort -n | awk '{a[NR]=$1} END {
 if (NR==0) exit 1;
 if(NR%2) printf "%.3f",a[(NR+1)/2];
 else printf "%.3f",(a[NR/2]+a[NR/2+1])/2;
}')
STATS=$(awk -F '\t' '{
 n++; sum+=$3;
 if($3>=20)p20++; if($3>=30)p30++; if($3>=50)p50++;
 if($3>=75)p75++; if($3>=100)p100++;
} END {
 if(!n)exit 1;
 printf "%.3f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f",
 sum/n,100*p20/n,100*p30/n,100*p50/n,100*p75/n,100*p100/n
}' "$TMP_DEPTH")
printf 'sample\tspecies\ttarget_depth\tseed\tsampling_fraction\toriginal_median_est\tmt_contig\tactual_median_depth\tactual_mean_depth\tPercent_20\tPercent_30\tPercent_50\tPercent_75\tPercent_100\tbam_path\n' > "$METRICS.tmp"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
 "$SAMPLE" "$SPECIES" "$TARGET" "$SEED" "$FRACTION" "$ORIGINAL_EST" \
 "$MT_CONTIG" "$MEDIAN" "$STATS" "$OUT" >> "$METRICS.tmp"
mv "$METRICS.tmp" "$METRICS"
echo "DONE: $METRICS"
