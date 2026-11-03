# HOM-selected cohort, both HOM and HET variants

**Clarification (2026-11-03):** HOM QC selects the **samples only**.
The downstream pipeline keeps the original variant-level filtering and emits
both HOM and HET calls for each HOM-eligible sample.
The filename `final_filter_hom_only` is a legacy directory name for this
HOM-selected cohort, **not** a HOM-only variant file.

The configuration `config/qc_preprocessing_coverage_hom_only.yaml` sets
`homoplasmy_qc.pipeline_call_class: HOM`,
`homoplasmy_qc.sample_selection_only: true`, and
`final_filter.biological_call_classes: HET,HOM`.

Eligibility: median mtCN >=20 and >=50% of mtDNA positions covered at >=50x.
For exactly zero nuclear median coverage, use nuclear median 1 when recalculating
median mtCN; otherwise preserve the original mtCN.

The sample cohort is selected from the existing
`sample_variant_filtering/reports/sample_qc.tsv`, where 2,346 samples
previously passed the HOM criteria. The array planner uses HOM-only sample
eligibility, and final_filter also checks HOM sample eligibility.

**Original final filtering:** In sample-selection-only mode, neither
HOM calls nor HET calls receive additional class-specific QC failures.
Existing VCF/DP, annotation, contamination, and other configured final
filters continue to apply. Both classes appear in
`final_filter_hom_only/reports/final_variant_qc.tsv`, including PASS and FAIL
records (for samples reaching final variant reporting). Output VCFs include
PASS calls from both classes.

```bash
cd /vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis_coverage_only
git pull --ff-only
python3 -m py_compile qc_analysis/scripts/run_final_filter_streaming.py
bash qc_analysis/scripts/run_qc_preprocessing.sh \
  --dry-run-submit --incremental --start-from pre_liftover_variant_qc all \
  config/qc_preprocessing_coverage_hom_only.yaml
bash qc_analysis/scripts/run_qc_preprocessing.sh \
  --submit --incremental --start-from pre_liftover_variant_qc all \
  config/qc_preprocessing_coverage_hom_only.yaml
```

If the per-sample annotations have already completed, **only the last step**
needs rerunning to apply the new final variant inclusion behavior:

```bash
bash qc_analysis/scripts/run_qc_preprocessing.sh --submit final_filter \
  config/qc_preprocessing_coverage_hom_only.yaml
```

This final-only command rebuilds the dedicated
`results/qc_coverage_only/final_filter_hom_only/` directory; backup its
reports first if previous output comparisons are needed.
The existing mixed project directory `results/qc_coverage_only/final_filter/`
is never targeted by the dedicated configuration.

Verify:

```bash
python3 - <<'PY'
import csv
p='results/qc_coverage_only/final_filter_hom_only/reports/final_variant_qc.tsv'
with open(p) as h:
    rows=list(csv.DictReader(h,delimiter='\t'))
print('Rows by source class:', {k:sum(r['source_call_class']==k for r in rows) for k in ('HOM','HET')})
print('PASS by source class:', {k:sum(r['source_call_class']==k and r['final_variant_status']=='PASS' and r['final_call_status']=='PASS' for r in rows) for k in ('HOM','HET')})
PY
```

The 2,346 is the number **sample-QC eligible**, not the guaranteed number
with complete liftover/annotations or final PASS calls.
Shared upstream intermediate directories must not be modified by two
concurrent pipelines.
