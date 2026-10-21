# Homoplasmic mtDNA downsampling pilot (McCleary)

This is **Stage 1 only**: stage WGS CRAM/CRAI, subsample WGS reads, and assess observed mtDNA depth. It does not run Round 1/Round 2 variant calling or establish HOM recall.

Inputs already available:
- `results/qc_coverage_only/downsampling_pilot/pilot_sample_ids.txt`
- `results/qc_coverage_only/downsampling_pilot/pilot_samples.tsv` (must include `sample`, `species`, `mt_median_coverage`)
- Original CRAM and CRAI: `/home/lt692/ycga_work/primate_downsampling_pilot/wgs_input/<sample>/<sample>.cram[.crai]`

Destination scratch root:
`/vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis_coverage_only/results/qc_coverage_only/downsampling_pilot/`

## Run

From the **coverage_only** repository root:

```bash
mkdir -p results/qc_coverage_only/downsampling_pilot/logs
bash qc_analysis/scripts/stage_pilot_wgs.sh
cat results/qc_coverage_only/downsampling_pilot/wgs_staging_report.tsv

bash qc_analysis/scripts/prepare_downsampling_manifest.sh
head results/qc_coverage_only/downsampling_pilot/downsampling_jobs.tsv

# Submit task 1 only before the full array.
bash qc_analysis/scripts/submit_downsampling.sh --pilot
```

If task 1 completes and the `*.metrics.tsv` actual median mt depth is sensible:

```bash
bash qc_analysis/scripts/submit_downsampling.sh
```

The array has 10 samples x 3 depth targets (30,50,100) x 2 random seeds (101,202) = 60 tasks, 6 concurrent. The first task may be resubmitted, but a completed BAM is reused.

Monitor:
```bash
squeue -u "$USER"
find results/qc_coverage_only/downsampling_pilot/downsampled -name '*.metrics.tsv'
```

**Critical biological checks:** The initial subsampling fraction uses mt median coverage from the existing pipeline, not necessarily the original CRAM's raw chrM depth; compare the measured depth and recalibrate if needed. A unique short contig candidate is not proof of mitochondrial identity; verify reference manifest. When auto detection is ambiguous, provide `MT_CONTIG` only when it is correct for *every* selected sample (references may have different contig names). When CRAM decoding needs an external reference, export a matching `CRAM_REFERENCE_FASTA`; a single FASTA may not match all species. The script intentionally fails instead of guessing. Never supply human hg38 unless it is truly the encoding reference.

The resulting site coverage is raw BAM depth on one mt contig; it is not equivalent to round2 merged or shifted max-depth coverage. CRAM source files under `/home` are **copied and retained**, never removed.

Before using these BAMs for HOM recovery analysis, confirm they are valid inputs to the actual Round 1/Round 2 workflows and relax the appropriate HOM-specific site DP filters in a separate QC config.
