# HOM-selected DP50 selective QC rerun

This optional workflow is for `config/qc_preprocessing_coverage_hom_only.yaml` after changing
the **variant-level** DP threshold to **DP >= 50**. It keeps sample eligibility
defined by **mtCN >= 20** and **at least 50% of mtDNA positions at >= 50x**
(the percentage comparison in current code is inclusive). It retains both
HOM and HET calls in HOM-eligible samples.

## What changes

- Recompute native-coordinate source QC at DP50.
- Recompute the complete local heteroplasmy / NUMT / indel artifact removal list.
- Refresh `SOURCE_VARIANT_QC` tags in **copies** of already annotated VCFs.
  Every row must match its native `SOURCE_CHROM/POS/REF/ALT` to a newly generated
  source QC row; source call class and DP must agree. Unverifiable rows abort.
- Rerun `final_filter` and post-final `homoplasmy_qc`.

What is reused: collected read/coverage/mtCN inputs, existing HOM sample-QC
eligibility, coordinate liftover/mapping, MITOS2 annotations and codon/tRNA/rRNA
matching decisions. The old outputs remain untouched.

The script does **not** regenerate the old collection summary's DP100 variant
counts. Those historical counts remain in the staged symlinked summary; for
DP50 reporting, use the freshly generated final variant report or perform a
separate collection recount. Do not treat those old count columns as DP50 counts.

## Run via Slurm

This is now a dependent Slurm job chain, with one source-QC array task per sample
(default concurrency 20) and separate merge, local-artifact, refresh, final and
HOM-QC jobs. The merged native report and VCF source lookup use bounded-memory
streaming/SQLite rather than an all-cohort Python dictionary.

From the repository root:

```bash
git pull --ff-only
bash qc_analysis/scripts/submit_hom_only_dp50_selective.sh
```

The launcher calls `--phase prepare` and *resumes* a pre-existing DP50
workspace left by an interrupted run. It then submits the new Slurm chain.
**Do not submit two DP50 chains simultaneously**. The completed array results
may be regenerated, while the original DP100 tree is not overwritten.

Logs: `results/qc_coverage_only_dp50/slurm_logs/`. Use `squeue -u "$USER"`
and `sacct -j JOBID --format=JobID,State,Elapsed,MaxRSS,ReqMem` to diagnose
OOM and job timing. 24h is a per-job walltime limit, not an estimated runtime.

Outputs:
- `results/qc_coverage_only_dp50/pre_liftover_variant_qc/`
- `results/qc_coverage_only_dp50/local_heteroplasmy_qc/`
- `results/qc_coverage_only_dp50/final_filter_hom_only/`
- `results/qc_coverage_only_dp50/homoplasmy_qc_only/`
- `config/qc_preprocessing_coverage_hom_only_dp50_selective.yaml`

Only the historical unchanged data directories are symlinked. Do not move
`results/qc_coverage_only` while the staged analysis exists.

## Important limitations / verification

- An input annotated VCF must retain the source identity tags. Any missing or
  conflicting native key fails closed rather than silently mixing QC versions.
- If the old annotated VCF is already missing DP50–99 variants, refreshing its
  tags cannot recover those variants. In that case a proper upstream rerun of
  liftover/functional annotations is required.
- Additional final filters (codon, tRNA, rRNA, orthology, artifact removal)
  remain active. This workflow is **not** equivalent to the four simplified
  criteria shown in a schematic.
- The current source code uses inclusive AF >= 0.10, DP >= 50, mtCN >= 20,
  per-position coverage >= 50, and percentage >= 50. Change the percentage
  comparison separately if strict >50% is required.

Check the final counts after the run:

```bash
python3 - <<'PY'
import csv
from pathlib import Path
root = Path('results/qc_coverage_only_dp50/final_filter_hom_only/reports')
for name in ('final_sample_qc.tsv', 'final_variant_qc.tsv'):
    path = root / name
    with path.open() as h:
        rows = list(csv.DictReader(h, delimiter='\t'))
    print(name, 'rows:', len(rows))
    if name.startswith('final_variant'):
        print('PASS variant rows:', sum(r.get('final_variant_status') == 'PASS' for r in rows))
PY
```

Compare with the untouched historical `results/qc_coverage_only/final_filter_hom_only/`.
