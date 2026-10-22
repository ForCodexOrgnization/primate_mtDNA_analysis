# Recommended accelerated experiment: proportional mtDNA + NUMT candidate downsampling

The **recommended** pilot now starts from each sample's existing Round-1 `EXTRACT_CANDIDATE_READS` output, `<sample>.with_mates.bam`. The older whole-WGS scripts are retained only as a comparison/control, not replaced. Do not submit their full array by mistake.

**Scientific interpretation:** The candidate BAM already selects primary read pairs from chrM + sample-specific NUMT intervals (both mates required). Subsampling by QNAME with one common fraction reduces chrM and NUMT-associated records proportionally *in expectation*. This does NOT independently keep nuclear depth constant, and raw BAM chrM depth is NOT the final Round-2 callable depth. It also does not regenerate read extraction from a downsampled WGS CRAM; compare the two methods on a small sample if their equivalence matters.

## Accelerated workflow on McCleary

From `/vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis_coverage_only`:

```bash
git pull origin main
BASE=results/qc_coverage_only/downsampling_pilot
# Must exist and include sample, species, mt_median_coverage:
head -n 3 "$BASE/pilot_samples.tsv"
# Inspect location of source candidate BAMs:
find results -name 'ERS12091908.with_mates.bam' -print
```

If candidate BAMs are outside this repository's `results/` folder, set `CANDIDATE_ROOT` to the **Round-1 pipeline outdir that directly contains `<sample>/round_1/candidate_reads`**. Do not set it to the candidate_reads directory itself. Candidate BAMs must already exist; if not, first run `EXTRACT_CANDIDATE_READS` once for each sample using the same reference and NUMT BED as the original run.

**Step 1: calibration once per sample (Slurm job):**

```bash
CANDIDATE_ROOT=/path/to/round1/pipeline/outdir \
  bash qc_analysis/scripts/submit_candidate_downsampling.sh --prepare
```

This computes source median depth using `samtools depth -aa -r chrM`, and writes `$BASE/candidate_downsampling_jobs.tsv` with 60 tasks and measured fractions. It requires `chrM` in every candidate BAM; if a sample uses another MT contig, prepare needs a sample-specific contig feature before including that sample. It does **not** derive fractions from the old QC summary.

Check `$BASE/logs/candidate_prepare_*.out` and `$BASE/logs/candidate_prepare_*.err`, and inspect the manifest.

**Step 2: one-task pilot, then full Slurm array:**

```bash
bash qc_analysis/scripts/submit_candidate_downsampling.sh --pilot
# After confirming output metrics and BAM:
MAX_CONCURRENT=4 bash qc_analysis/scripts/submit_candidate_downsampling.sh --all
```

Each task uses one `samtools view --subsample-seed --subsample` over an already-filtered candidate BAM, preserving read pair consistency. Output is `$BASE/candidate_downsampled/<sample>/` with an indexed BAM and metrics TSV (including actual raw chrM median/mean and mapped chrM/non-chrM counts before and after).

**Important:** Test task 1 is part of the 60-task manifest; full array will encounter it again and reuse an existing BAM. Do not run full array while pilot is still active. The full **variant-calling** repository now implements the Round-1 restart-from-candidate input; the analysis repository only prepares downsampled BAMs. See `ForCodexOrgnization/primate_mt_variant_calling/scripts/launch_hom_downsampling_round1.sh`. Do not use the stale `variant_calling/` mirror in this analysis repo. Run the single-sample calling pilot from the complete variant-calling repository, then wire Round 2 using its authoritative outputs. The initial raw coverage fraction may not lead to exact desired median because sampling is stochastic, and especially not to precise final calling coverage. Always analyze actual measured depth.

## Legacy whole-WGS experiment (optional control)


This is **Stage 1 only**: directly read WGS CRAM/CRAI from `/home`, subsample WGS reads with Slurm arrays, and write BAM/depth metrics to scratch. It does not run Round 1/Round 2 variant calling or establish HOM recall.

Inputs already available:
- `results/qc_coverage_only/downsampling_pilot/pilot_sample_ids.txt`
- `results/qc_coverage_only/downsampling_pilot/pilot_samples.tsv` (must include `sample`, `species`, `mt_median_coverage`)
- Original CRAM and CRAI: `/home/lt692/ycga_work/primate_downsampling_pilot/wgs_input/<sample>/<sample>.cram[.crai]`

BAM and metrics output scratch root:
`/vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis_coverage_only/results/qc_coverage_only/downsampling_pilot/`

## Run

From the **coverage_only** repository root:

```bash
mkdir -p results/qc_coverage_only/downsampling_pilot/logs
bash qc_analysis/scripts/prepare_downsampling_manifest.sh
head results/qc_coverage_only/downsampling_pilot/downsampling_jobs.tsv

# Submit task 1 only before the full array.
bash qc_analysis/scripts/submit_downsampling.sh --pilot
```

If task 1 completes and the `*.metrics.tsv` actual median mt depth is sensible:

```bash
bash qc_analysis/scripts/submit_downsampling.sh
```

The array has 10 samples x 3 depth targets (30,50,100) x 2 random seeds (101,202) = 60 tasks, **4 concurrent by default**, configurable with `MAX_CONCURRENT=6 bash qc_analysis/scripts/submit_downsampling.sh`. Each task rereads the original CRAM; monitor shared filesystem load and start conservatively. The first task may be resubmitted, but a completed BAM is reused.

Monitor:
```bash
squeue -u "$USER"
find results/qc_coverage_only/downsampling_pilot/downsampled -name '*.metrics.tsv'
```

**Critical biological checks:** The initial subsampling fraction uses mt median coverage from the existing pipeline, not necessarily the original CRAM's raw chrM depth; compare the measured depth and recalibrate if needed. A unique short contig candidate is not proof of mitochondrial identity; verify reference manifest. When auto detection is ambiguous, provide `MT_CONTIG` only when it is correct for *every* selected sample (references may have different contig names). When CRAM decoding needs an external reference, export a matching `CRAM_REFERENCE_FASTA`; a single FASTA may not match all species. The script intentionally fails instead of guessing. Never supply human hg38 unless it is truly the encoding reference.

The resulting site coverage is raw BAM depth on one mt contig; it is not equivalent to round2 merged or shifted max-depth coverage. Input CRAM/CRAI are **read directly from `/home`, unchanged**. All generated BAM, BAI and metrics are written to scratch. `stage_pilot_wgs.sh` remains available as an optional, separate backup/staging operation and is not needed for the main workflow. If a large `stage_pilot_wgs.sh` job is currently running, cancel that job or interrupt it before launching this workflow to avoid duplicate I/O. Incomplete scratch copies can be left alone until inspected.

Before using these BAMs for HOM recovery analysis, confirm they are valid inputs to the actual Round 1/Round 2 workflows and relax the appropriate HOM-specific site DP filters in a separate QC config.

## Troubleshooting: many short nuclear contigs

Earlier versions failed on fragmented primate reference assemblies with thousands of 14-20 kb nuclear scaffolds (for example, 6,795 for ERS12091908). The updated Slurm worker ranks 14-20 kb contigs by mapped reads per base and requires the top candidate to have >=100 mapped reads/base and >=20-fold lead over the second candidate. It prints only the top ten. This is a **candidate selection**, not definitive mtDNA identification: verify the selected contig against the original species reference manifest before interpreting biological results. If ambiguous, the worker stops. The safe manual override is to run one species' pilot task with `MT_CONTIG=<correct_contig> bash qc_analysis/scripts/submit_downsampling.sh --pilot`. Do not export a single MT_CONTIG for a multi-species batch unless the same contig identifier is correct for all samples.

Re-run after pulling updated `main`:

```bash
git pull origin main
bash qc_analysis/scripts/submit_downsampling.sh --pilot
```

No CRAM/CRAI need to be copied. The original failed task quit before downsampling.


## Repository ownership for HOM pilot (2026-10-22)

- **Analysis, extraction, calibration, subsampling:** `ForCodexOrgnization/primate_mtDNA_analysis` (`qc_analysis/scripts/*candidate*`). Output: `$BASE/candidate_downsampled/<sample>/<sample>.candidate.ds<target>x.seed<seed>.bam` and its index.
- **Round 1/2 variant calling:** `ForCodexOrgnization/primate_mt_variant_calling` (complete production Nextflow/WDL pipeline). The optional Round 1 parameters `--candidate_bam_root`, `--candidate_target_depth`, and `--candidate_seed` bypass `EXTRACT_CANDIDATE_READS` without changing the normal production path.

Example after generating the indexed **100x** BAM:

```bash
cd /home/lt692/project_pi_njl27/lt692/primate_mt_variant_calling
git pull origin main
bash scripts/launch_hom_downsampling_round1.sh ERS12091908 100 101
```

Use a distinct output directory per target/seed and keep `--consensus_min_dp` fixed/recorded unless intentionally testing altered consensus rules. **100x here means raw candidate depth, not the final Round-2 callable depth.** The Round-2 handoff and HOM comparison are not yet automated; do not interpret Round-1 completion alone as validation. The original analysis repo's `variant_calling/` folder is an incomplete mirror and is not the execution source for this experiment.
