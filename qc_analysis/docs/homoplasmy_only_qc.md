# Homoplasmy-only QC from previously generated results

Use `homoplasmy_qc` **after** `final_filter`. It does not rerun FASTQ alignment, variant calling, liftover, or tRNA annotation, and does not modify existing final-filter files.

## Eligibility

For each sample:
- Adjust `nuclear_median_coverage` from exactly zero to **1**; otherwise retain the original median.
- For original nuclear median **exactly zero**, recalculate `mtcn_median_recalculated = 2 * mt_median_coverage / 1`; for all other samples retain the original `mtcn_median` without recalculating it. Missing mt coverage is not imputed.
- Calculate `Percent_50` from `collected_cov/{sample}.merged.max_depth.per_base_coverage.tsv` (number of covered positions at depth >=50 divided by number of recorded positions, percentage). This uses the existing merged maximum-depth coverage definition. At least 14,000 recorded positions are required by default as a coverage-file completeness guard.
- A homoplasmic candidate passes if `mtCN >=20` **and** `Percent_50 >=50`, both inclusive.
- These new criteria affect **only** previously FINAL PASS records with `source_call_class=HOM`. HET calls and prior FAIL records remain unchanged.

Outputs are written to a new results directory:
- `reports/homoplasmy_sample_qc.tsv` (input medians, corrected nuclear median, recalculated mtCN, Percent_50, sample gate)
- `reports/final_variant_qc.homoplasmy_adjusted.tsv` (all original rows and extra homoplasmy-adjusted status columns)
- Optionally `final_vcf/*.final.vcf.gz` and tabix indexes when using the Python script with `--write-vcf`.

**Important limitation:** this is a *post-final* gate. It cannot restore candidate HOM variants for samples previously excluded by strict sample QC, because those samples are absent from the earlier final report. To genuinely relax prior sample QC, upstream sample eligibility and downstream per-sample work must be revisited.

## Run on McCleary

```bash
cd /vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis_coverage_only
git pull --ff-only
python3 -m py_compile qc_analysis/scripts/run_homoplasmy_qc.py
bash qc_analysis/scripts/run_qc_preprocessing.sh --submit homoplasmy_qc config/qc_preprocessing_coverage_only.yaml
```

Or for incremental strict:

```bash
cd /vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis_incremental
git pull --ff-only
bash qc_analysis/scripts/run_qc_preprocessing.sh --submit homoplasmy_qc config/qc_preprocessing_incremental_strict.yaml
```

To write final VCFs as well (run directly on a compute allocation, not a login node):
```bash
python3 qc_analysis/scripts/run_homoplasmy_qc.py --config config/qc_preprocessing_coverage_only.yaml --write-vcf
```
