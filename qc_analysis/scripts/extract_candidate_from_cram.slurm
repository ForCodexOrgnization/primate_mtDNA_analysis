#!/usr/bin/env bash
#SBATCH --job-name=primate_extract_candidate
#SBATCH --partition=day
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --output=results/qc_coverage_only/downsampling_pilot/logs/extract_candidate_%A_%a.out
#SBATCH --error=results/qc_coverage_only/downsampling_pilot/logs/extract_candidate_%A_%a.err
set -euo pipefail
module load SAMtools/1.21-GCC-13.3.0
ROOT=/vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis_coverage_only
BASE=$ROOT/results/qc_coverage_only/downsampling_pilot
PILOT=$BASE/pilot_samples.tsv
TASK=$(printf '%s' "$SLURM_ARRAY_TASK_ID")
SAMPLE=$(awk -F '\t' -v n="$TASK" 'NR==n+1 {print $1; exit}' "$PILOT")
[[ -n "$SAMPLE" ]] || { echo "Missing sample for task $TASK" >&2; exit 1; }
WGS_ROOT=${WGS_INPUT_ROOT:-/home/lt692/ycga_work/primate_downsampling_pilot/wgs_input}
META_ROOT=${CANDIDATE_METADATA_ROOT:-/home/lt692/ycga_work/primate_downsampling_pilot/candidate_input}
OUT_ROOT=${CANDIDATE_ROOT:-$BASE/reextracted_candidate}
CRAM=$WGS_ROOT/$SAMPLE/$SAMPLE.cram
BED=$META_ROOT/$SAMPLE/round_1/candidate_reads/$SAMPLE.target_chrM_plus_NUMT.bed
ORIG_SUMMARY=$META_ROOT/$SAMPLE/round_1/candidate_reads/$SAMPLE.candidate_read_selection.summary.tsv
OUTDIR=$OUT_ROOT/$SAMPLE/round_1/candidate_reads
OUT=$OUTDIR/$SAMPLE.with_mates.bam
T=${SLURM_CPUS_PER_TASK:-4}
[[ -s "$CRAM" && -s "$CRAM.crai" && -s "$BED" ]] || { echo "Missing CRAM/CRAI/BED: $CRAM $BED" >&2; exit 1; }
mkdir -p "$OUTDIR"
if [[ -s "$OUT" && -s "$OUT.bai" ]]; then
  samtools quickcheck -v "$OUT"
  echo "Already complete $OUT"
  exit 0
fi
REF=()
if [[ -n ${CRAM_REFERENCE_FASTA:-} ]]; then
  [[ -s "$CRAM_REFERENCE_FASTA" ]] || { echo "Missing CRAM reference FASTA" >&2; exit 1; }
  REF=(-T "$CRAM_REFERENCE_FASTA")
fi
# Avoid competing tasks or accidental stale partial output.
WORK=$OUTDIR/.work_$SAMPLE.$SLURM_JOB_ID
mkdir "$WORK"
trap 'rm -rf "$WORK"' EXIT
echo "SAMPLE=$SAMPLE INPUT=$CRAM BED=$BED OUTPUT=$OUT"
samtools quickcheck -v "$CRAM"
echo "samtools=$(command -v samtools)"
samtools --version | head -n 2
# -M is supported in current samtools but may not be available in the
# login node's older conda binary. Prefer indexed multi-region iteration
# when supported; fall back safely to -L if this binary lacks -M.
EXTRACT_REGION_ARGS=(-L "$BED")
if samtools view --help 2>&1 | grep -Eq -- '(^|[[:space:]])-M([,[:space:]]|$)|--use-index'; then
  EXTRACT_REGION_ARGS=(-M -L "$BED")
  echo "Using indexed -M -L extraction"
else
  echo "WARNING: samtools binary lacks -M; using slower -L BED filtering"
fi
samtools view -@ "$T" "${REF[@]}" -b -F 0x900 "${EXTRACT_REGION_ARGS[@]}" \
  -o "$WORK/target.bam" "$CRAM"
samtools sort -n -@ "$T" -o "$WORK/name.bam" "$WORK/target.bam"
# Reproduce the original Round-1 QNAME-group filter (n >= 2 primary records).
samtools view -h "$WORK/name.bam" |
  awk 'BEGIN{FS=OFS="\t"}
       /^@/{print;next}
       {q=$1
        if(curr==""){curr=q;n=0}
        if(q!=curr){if(n>=2){for(i=1;i<=n;i++)print rec[i]}
                    delete rec;n=0;curr=q}
        rec[++n]=$0
       }
       END{if(n>=2){for(i=1;i<=n;i++)print rec[i]}}' |
  samtools view -@ "$T" -b -o "$WORK/pairs.bam" -
samtools sort -@ "$T" -o "$WORK/final.bam" "$WORK/pairs.bam"
samtools index -@ "$T" "$WORK/final.bam"
samtools quickcheck -v "$WORK/final.bam"
COUNT=$(samtools view -c "$WORK/final.bam")
(( COUNT>0 )) || { echo "ERROR: no paired candidate reads" >&2; exit 1; }
mv "$WORK/final.bam" "$OUT"
mv "$WORK/final.bam.bai" "$OUT.bai"
echo "DONE: $OUT records=$COUNT original_summary=$ORIG_SUMMARY"
echo "NOTE: Input BED is from prior Run 1; reference FASTA must match the CRAM if required for decoding."
