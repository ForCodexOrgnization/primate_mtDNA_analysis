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

Each submission gets immutable, UTC-timestamped `.tasks.txt`, `.manifest.tsv`,
and `.submission.tsv` files under that step's biological output root. Path
resolution uses `paths.job_array_dir` when configured, then
`paths.output_dir/job_arrays`, then the parent of `paths.reports_dir` plus
`job_arrays`. Steps without those configured locations use their documented
step defaults; unknown future steps fall back to
`.workflow/qc_preprocessing/<step>`. Candidate samples come respectively from the
liftover sample/reference table, the codon sample/reference map plus lifted
VCFs, and the configured tRNA/rRNA fallback chains.  Valid completed annotated
VCFs are skipped by default; use `FORCE_RERUN=true` (or
`SKIP_COMPLETED=false`) to include them. `--prepare-retry --prepare-only`
creates a timestamped retry list and manifest containing missing or invalid
outputs without submission. A submitted worker always receives that immutable
timestamped task-file path; `<step>.current.tsv` is only a convenience pointer
and is never read by running jobs.

Codon, tRNA, and MITOS2 producer arrays automatically submit their singleton merge
with `afterok` unless `AUTO_SUBMIT_MERGE=false`. `--submit all` creates the
ordered dependency graph printed by the wrapper rather than one monolithic
job. Logs resolve in this order: `SLURM_LOG_DIR`, configured `paths.log_dir`,
then `<step output>/logs/job_arrays`. `SLURM_LOG_DIR` remains a global override
but is no longer required for a useful layout.

For example, codon matching keeps biological results, workflow metadata, and
scheduler logs separate:

```text
results/qc/codon_match/
├── vcf_codon/
├── reports/
├── job_arrays/
└── logs/job_arrays/
```

The wrapper prints the resolved step, output directory, immutable task file,
manifest, log pattern, task count, and array expression. Paths containing
spaces are passed as individual shell arguments and remain supported.

Historical files under `results/qc/job_arrays/` are intentionally not moved or
deleted because pending or running arrays may still reference them. All new
submissions use the step-specific layout. After confirming with Slurm that no
old arrays are active, administrators may archive those files manually. No
automatic migration is provided because safe reference detection is site- and
scheduler-dependent.

```bash
# Validate codon inputs (singleton array)
bash qc_analysis/scripts/run_qc_preprocessing.sh --submit codon_match_validate config/qc_preprocessing.yaml
# All eligible samples, default 20-way limit
bash qc_analysis/scripts/run_qc_preprocessing.sh --submit codon_match config/qc_preprocessing.yaml
# Explicitly merge per-sample tRNA summaries (automatic after array success by default)
bash qc_analysis/scripts/run_qc_preprocessing.sh trna_match_merge config/qc_preprocessing.yaml
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
