# Mean-mtDNA mtCN cohort QC rerun

Both cohorts use the **same** mtCN definition:
`mtcn_qc_value = 2 * mt_mean_coverage / nuclear_median_coverage`.
The original `mtcn_median` field is preserved as historical provenance, NOT used
as the new QC threshold variable.

For nuclear median coverage exactly 0, the existing HOM workflow convention
adjusts the denominator to 1x and reports the adjustment in
`mtcn_qc_nuclear_zero_adjusted`. Missing/negative depths do not pass.

| Rule | HOM10 | HOMHET40 |
| --- | --- | --- |
| mtCN | >=10 | >=40 |
| mtDNA coverage | >=50% of sites at >=50x | >90% of sites at >=100x |
| Nuclear median coverage | no independent minimum | >=5x |
| Variant source DP | >=50 | >=100 |
| source HET AF | >=0.10 | >=0.10 |
| Variant classes retained | HOM + HET | HOM + HET |

The HOM10 name means HOM *sample-selection* rules, not HOM-only variant output.

## Running on McCleary (from repo root)

```bash
git pull --ff-only
python3 -m py_compile qc_analysis/lib/mtcn_qc.py qc_analysis/scripts/run_sample_variant_filtering.py qc_analysis/scripts/run_homoplasmy_qc.py qc_analysis/scripts/run_final_filter_streaming.py qc_analysis/scripts/prepare_mean_mtcn_cohorts.py qc_analysis/scripts/report_mean_mtcn_cohorts.py
bash -n qc_analysis/scripts/submit_mean_mtcn_cohorts.sh
python3 -m unittest qc_analysis.tests.test_mean_mtcn_cohorts
python3 qc_analysis/scripts/prepare_mean_mtcn_cohorts.py --cohort both
bash qc_analysis/scripts/submit_mean_mtcn_cohorts.sh hom10
bash qc_analysis/scripts/submit_mean_mtcn_cohorts.sh homhet40
squeue -u "$USER"
```

The submit scripts only **submit** jobs, not run analyses on the login node.
Use `sacct -j JOBID --format=JobID,State,ExitCode` and inspect each log in the
cohort's `slurm_logs/` before trusting outputs. Dependency chains use
`afterok`, so failed stages prevent downstream finalization.

## Isolation / upstream reuse

Results:
- `results/qc_cohort_hom10_mean_mtcn/`
- `results/qc_cohort_homhet40_mean_mtcn/`

The preparation script links original collected inputs, alignment/liftover, and
functional annotations without overwriting them. It points each cohort's rRNA
VCF input at the independently refreshed files in
`results/qc_coverage_only_dp50/rrna_match/vcf_rrna`, preserving SOURCE tags.

Each run recomputes cohort sample QC, local artifact detection, and final
filtering. HOM10 uses the already completed DP50 native SOURCE QC; HOMHET40
recomputes native SOURCE QC at DP100 (array jobs + merge) before local artifacts.

Reports:
- `sample_variant_filtering/reports/sample_qc.tsv`
- `final_filter_hom_only/reports/final_sample_qc.tsv`
- `final_filter_hom_only/reports/final_variant_qc.tsv`
- HOM10 only: `homoplasmy_qc_only/reports/homoplasmy_sample_qc.tsv`

After **both** final jobs complete:

```bash
python3 qc_analysis/scripts/report_mean_mtcn_cohorts.py
```

The report writes a sample/species/variant comparison TSV, plus a boxplot of
mean mtDNA depth for newly admitted HOM10 samples with 10 <= mtCN < 20 versus
HOM10 samples with mtCN >= 20. Both plotting groups are subject to the other
HOM10 sample QC conditions.

## Important limitations

This is a **selective QC rerun**, not a new read alignment, variant calling,
coordinate liftover, or annotation. It intentionally reuses previously generated
VCF records and coordinates. Existing sample-specific contamination reports
remain subject to the two cohort configs (currently disabled as sample gates).

A fresh cohort sample count may differ from a previous run using
`mtcn_median`. Report counts only after the Slurm tasks complete; none are
predicted in advance.
