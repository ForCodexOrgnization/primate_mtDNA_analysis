#!/usr/bin/env bash
# Run ONE existing candidate-downsampled BAM through Round 1 of the
# standard chrM+NUMT decoy pipeline. Does not start WGS extraction.
# Example: bash qc_analysis/scripts/submit_hom_round1_pilot.sh ERS12091908 100 101
set -euo pipefail
ROOT=/vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis_coverage_only
cd "$ROOT"
SAMPLE=${1:-ERS12091908}
TARGET=${2:-100}
SEED=${3:-101}
[[ "$SAMPLE" =~ ^[A-Za-z0-9_]+$ && "$TARGET" =~ ^[0-9]+$ && "$SEED" =~ ^[0-9]+$ ]] || { echo "Invalid sample/target/seed" >&2; exit 2; }
BASE="$ROOT/results/qc_coverage_only/downsampling_pilot"
BAM_ROOT="$BASE/candidate_downsampled"
BAM="$BAM_ROOT/$SAMPLE/$SAMPLE.candidate.ds${TARGET}x.seed${SEED}.bam"
[[ -s "$BAM" && -s "$BAM.bai" ]] || { echo "Missing downsampled candidate BAM or index: $BAM" >&2; exit 1; }
SUM="$BASE/pilot_samples.tsv"
SPECIES=$(awk -F '\t' -v s="$SAMPLE" 'NR>1 && $1==s{print $2;exit}' "$SUM")
REF_NAME=$(awk -F '\t' -v s="$SAMPLE" '$1==s{print $2;exit}' config/sample_ref_file.tsv)
[[ -n "$SPECIES" && -n "$REF_NAME" ]] || { echo "No species/ref mapping for $SAMPLE" >&2; exit 1; }
mkdir -p "$BASE/round1_pilot_inputs/$SAMPLE"
PILOT_INPUT="$BASE/round1_pilot_inputs/$SAMPLE"
SAMPLE_TSV="$PILOT_INPUT/sample.tsv"
printf '%s\t%s\t%s\n' "$SAMPLE" "$SPECIES" "$REF_NAME" > "$SAMPLE_TSV"
# Supply the original sample-specific NUMT intervals from the transferred
# compound BED. chrM is removed because PREPARE_DECOY_REFERENCE builds chrM itself.
BED_SRC="/home/lt692/ycga_work/primate_downsampling_pilot/candidate_input/$SAMPLE/round_1/candidate_reads/$SAMPLE.target_chrM_plus_NUMT.bed"
[[ -s "$BED_SRC" ]] || { echo "Missing target BED: $BED_SRC" >&2; exit 1; }
BED_DIR="$PILOT_INPUT/numt_bed"
mkdir -p "$BED_DIR"
awk 'BEGIN{OFS="\t"} $1!="chrM" && NF>=3 && $3>$2 {print $1,$2,$3}' "$BED_SRC" > "$BED_DIR/$SAMPLE.highconf_numt.bed"
# LOCATE_CRAM still supplies source metadata, but no full-WGS extraction occurs.
CRAM_SRC="/home/lt692/ycga_work/primate_downsampling_pilot/wgs_input/$SAMPLE/$SAMPLE.cram"
[[ -s "$CRAM_SRC" && -s "$CRAM_SRC.crai" ]] || { echo "Missing original CRAM/CRAI" >&2; exit 1; }
CRAM_DIR="$PILOT_INPUT/cram_inputs/$SAMPLE/alignment"
mkdir -p "$CRAM_DIR"
ln -sfn "$CRAM_SRC" "$CRAM_DIR/$SAMPLE.cram"
ln -sfn "$CRAM_SRC.crai" "$CRAM_DIR/$SAMPLE.cram.crai"
# Only run if standard refs and applications are accessible in this environment.
GLOBAL_REF=${PILOT_GLOBAL_REF_DIR:-/nfs/roberts/pi/pi_njl27/lt692/primate_ref_files/whole_genomes_final}
MT_REF=${PILOT_MT_REF_DIR:-/nfs/roberts/pi/pi_njl27/lt692/primate_ref_files/Refseq_chrM}
[[ -s "$GLOBAL_REF/$REF_NAME.fasta" ]] || { echo "Missing whole-genome reference: $GLOBAL_REF/$REF_NAME.fasta" >&2; exit 1; }
[[ -s "$MT_REF/$REF_NAME.fasta" ]] || { echo "Missing mitochondrial reference: $MT_REF/$REF_NAME.fasta" >&2; exit 1; }
OUT="$BASE/variant_calling_pilot/ds${TARGET}x_seed${SEED}/round1"
WORK="$BASE/nf_work/ds${TARGET}x_seed${SEED}/round1"
mkdir -p "$OUT" "$WORK" "$BASE/logs"
# Explicit user-controlled consensus threshold; default is ORIGINAL 100.
# Never silently relax this in a scientific comparison.
CONSENSUS_DP=${CONSENSUS_MIN_DP:-100}
[[ "$CONSENSUS_DP" =~ ^[0-9]+$ ]] || exit 2
JOB="$PILOT_INPUT/run_round1.slurm"
cat > "$JOB" <<'SLURM'
#!/usr/bin/env bash
#SBATCH --job-name=hom_ds_round1
#SBATCH --partition=day
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
set -euo pipefail
module load Nextflow/24.10.2
cd "$ROOT"
nextflow run variant_calling/primate_pipeline_numt_decoy_round1.nf -profile cluster \
  -c variant_calling/nextflow.config -w "$WORK" \
  --sample_tsv "$SAMPLE_TSV" --outdir "$OUT" \
  --cram_dirs "$PILOT_INPUT/cram_inputs" \
  --numt_bed_dir "$BED_DIR" \
  --candidate_bam_root "$BAM_ROOT" \
  --candidate_target_depth "$TARGET" \
  --candidate_seed "$SEED" \
  --consensus_min_dp "$CONSENSUS_DP" \
  --global_ref_dir "$GLOBAL_REF" --ref_dir "$MT_REF"
SLURM
# Variables are expanded by a launcher-specific exported environment, not by
# a quoted heredoc; ensure environment reaches sbatch with --export=ALL.
export ROOT WORK SAMPLE_TSV OUT PILOT_INPUT BED_DIR BAM_ROOT TARGET SEED CONSENSUS_DP GLOBAL_REF MT_REF
echo "Pilot sample=$SAMPLE reference=$REF_NAME target=$TARGET seed=$SEED consensus_min_dp=$CONSENSUS_DP"
echo "Output=$OUT"
sbatch --export=ALL --output="$BASE/logs/hom_round1_${SAMPLE}_ds${TARGET}_seed${SEED}_%j.out" --error="$BASE/logs/hom_round1_${SAMPLE}_ds${TARGET}_seed${SEED}_%j.err" "$JOB"
