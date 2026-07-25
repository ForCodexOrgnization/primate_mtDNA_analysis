# QC preprocessing Slurm arrays

`run_qc_preprocessing.sh --submit` uses one submission path for every step.
`coordinate_liftover`, `codon_match`, `trna_match`, and `rrna_match` are
sample arrays; `mitos2_annotation` is a reference array.  Collection, global
anchor discovery, codon-table construction, comparison, validation, merging,
MITOS2 preparation/merge, and contamination QC are singleton arrays because
they own shared outputs.  Codon-table construction deliberately remains a
singleton until its global tables can be split into independently writable
intermediates.

The default is `SLURM_ARRAY_CONCURRENCY=20` (`1-N%20`).  Override it with an
environment variable or `--array-concurrency`; the value must be positive.
`CODON_MATCH_`, `LIFTOVER_`, `TRNA_MATCH_`, `RRNA_MATCH_`, and `MITOS2_`
prefixes override `SLURM_TIME`, `SLURM_MEM`, and `SLURM_CPUS` for their steps.

Each submission gets immutable `.tasks.txt` and `.manifest.tsv` files in
`results/qc/job_arrays`.  Candidate samples come respectively from the
liftover sample/reference table, the codon sample/reference map plus lifted
VCFs, and the configured tRNA/rRNA fallback chains.  Valid completed annotated
VCFs are skipped by default; use `FORCE_RERUN=true` (or
`SKIP_COMPLETED=false`) to include them. `--prepare-retry --prepare-only`
creates a list containing missing or invalid outputs without submission.

Codon and MITOS2 producer arrays automatically submit their singleton merge
with `afterok` unless `AUTO_SUBMIT_MERGE=false`. `--submit all` creates the
ordered dependency graph printed by the wrapper rather than one monolithic
job. Per-task logs are `logs/qc_preprocessing/<step>/%A_%a.{out,err}`.

```bash
# Validate codon inputs (singleton array)
bash qc_analysis/scripts/run_qc_preprocessing.sh --submit codon_match_validate config/qc_preprocessing.yaml
# All eligible samples, default 20-way limit
bash qc_analysis/scripts/run_qc_preprocessing.sh --submit codon_match config/qc_preprocessing.yaml
# Override concurrency
SLURM_ARRAY_CONCURRENCY=40 bash qc_analysis/scripts/run_qc_preprocessing.sh --submit codon_match config/qc_preprocessing.yaml
# Exactly one sample (array 1-1)
bash qc_analysis/scripts/run_qc_preprocessing.sh --submit --sample SAMPLE_NAME codon_match config/qc_preprocessing.yaml
# Inspect without sbatch
bash qc_analysis/scripts/run_qc_preprocessing.sh --dry-run-submit codon_match config/qc_preprocessing.yaml
# Full afterok graph
bash qc_analysis/scripts/run_qc_preprocessing.sh --submit all config/qc_preprocessing.yaml
```

Without `--submit`, existing direct commands still execute immediately. A
direct command without `--sample` may write cohort summaries sequentially;
array workers never write those merged summaries.
