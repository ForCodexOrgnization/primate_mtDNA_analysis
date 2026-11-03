# HOM-only coverage-only rerun

**Purpose:** Analyze only samples meeting the revised HOM gate, without processing
HET-only samples or changing existing mixed/HET final results.

The dedicated configuration `config/qc_preprocessing_coverage_hom_only.yaml`
inherits the original coverage-only upstream paths, but selects
`homoplasmy_qc.pipeline_call_class: HOM` and sets
`final_filter.biological_call_classes: HOM`.

HOM eligibility is `median mtCN >=20` and `Percent_50 >=50`, where
`Percent_50` denotes positions with depth >=50x. When median nuclear
coverage is exactly zero, the denominator is set to 1 and
`median mtCN = 2 * median mitochondrial coverage / 1`.

This config **reuses** the existing
`results/qc_coverage_only/sample_variant_filtering/reports/sample_qc.tsv`
created from the 3,332 sample inventory (2,346 HOM-pass samples).
Only HOM-eligible samples are eligible for incremental per-sample tasks. Final
filtering restricts samples to HOM PASS and skips all HET calls entirely.

Before running ensure the earlier sample-QC job has finished and the report
has the `homoplasmy_qc_status` column. The `--incremental` task scheduler
additionally requires collection `status=OK`, and later steps require
successful liftover, usable references and available source variants; thus
2,346 is the cohort *eligibility* count, not a guaranteed final output count.

```bash
cd /vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis_coverage_only
git pull --ff-only
bash qc_analysis/scripts/run_qc_preprocessing.sh \
  --dry-run-submit --incremental --start-from pre_liftover_variant_qc all \
  config/qc_preprocessing_coverage_hom_only.yaml
bash qc_analysis/scripts/run_qc_preprocessing.sh \
  --submit --incremental --start-from pre_liftover_variant_qc all \
  config/qc_preprocessing_coverage_hom_only.yaml
```

The new HOM-only output paths are
`results/qc_coverage_only/final_filter_hom_only/reports/final_variant_qc.tsv`,
`results/qc_coverage_only/final_filter_hom_only/final_vcf/`, and
`results/qc_coverage_only/homoplasmy_qc_only/reports/`.

The prior `results/qc_coverage_only/final_filter/` is **not overwritten**.
Upstream annotation, merged reports, and intermediate output paths are shared
with the old workspace; do not run concurrent mixed and HOM-only workflows
against the same directories. The HOM-only command may recompute missing
annotation outputs for newly eligible samples and regenerate global reports.

After completion, check:

```bash
python3 - <<'PY'
import csv
from pathlib import Path
p=Path("results/qc_coverage_only/final_filter_hom_only/reports/final_variant_qc.tsv")
with p.open() as h:
    reader=csv.DictReader(h,delimiter="\t")
    rows=list(reader)
passing=[r for r in rows if r["final_variant_status"]=="PASS"
         and r["final_call_status"]=="PASS"]
print("final passing variant rows:",len(passing))
print("final passing samples:",len({r["sample"] for r in passing}))
print("HET rows:",sum(r["source_call_class"]=="HET" for r in rows))
assert all(r["source_call_class"]=="HOM" for r in rows), "Non-HOM record found"
PY
```
