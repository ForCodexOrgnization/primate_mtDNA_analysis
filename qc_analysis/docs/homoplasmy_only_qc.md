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

**Standalone post-final mode limitation:** Running only `homoplasmy_qc` cannot restore excluded samples. Use the full-cohort reassessment above for those samples.


## Full-cohort reassessment (recommended for revised criteria)

Both `config/qc_preprocessing_incremental_strict.yaml` and `config/qc_preprocessing_coverage_only.yaml`
enable `homoplasmy_qc.reassess_all_samples: true`.

The sample-QC stage writes **independent** `heteroplasmy_qc_status` and
`homoplasmy_qc_status`. A sample is admitted to downstream annotation if either
class passes; `final_filter` then rejects individual HOM/HET variants when their
respective class gate fails. The old `final_variant_status` in this *new run*
therefore already represents class-specific final QC, and the original
`final_call_status` still describes whether the sample and variant passed.

From a **complete, previously staged collection** (not only old final PASS
samples), back up prior `final_filter` and `sample_variant_filtering` results
before recomputing. The recommended Slurm workflow is:

```bash
cd /vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis_coverage_only
git pull --ff-only
bash qc_analysis/scripts/run_qc_preprocessing.sh --submit sample_variant_filtering config/qc_preprocessing_coverage_only.yaml
# Check the new reports/sample_qc.tsv, then submit the resumed chain:
bash qc_analysis/scripts/run_qc_preprocessing.sh --dry-run-submit --incremental --start-from pre_liftover_variant_qc all config/qc_preprocessing_coverage_only.yaml
bash qc_analysis/scripts/run_qc_preprocessing.sh --submit --incremental --start-from pre_liftover_variant_qc all config/qc_preprocessing_coverage_only.yaml
```

Do **not** independently submit the second command before sample_variant_filtering
has finished. `--incremental` reuses existing per-sample outputs where they
validate; missing results for newly eligible samples are computed. The global
merge stages and the terminal final filter are rerun.

For incremental strict use the corresponding config and project root.

An explicit `homoplasmy_qc` post-final step is still available to produce
separate audit tables. It does not replace upstream eligibility reassessment.

**Caution**: source collection status, liftover completion, reference maps,
missing coverage, and other biological/annotation QC gates still apply. A new
HOM sample eligibility PASS does not guarantee a passing variant or that all
missing inputs can be reconstructed.

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
