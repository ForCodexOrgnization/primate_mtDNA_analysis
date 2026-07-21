# MITOS2 annotation integration

MITOS2 is invoked only through its conda environment, not as a presumed system command:

```bash
module load miniconda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mitos2

echo "CONDA_PREFIX=$CONDA_PREFIX"
echo "MITOS executable=$(command -v runmitos || true)"
if ! command -v runmitos >/dev/null 2>&1; then
    echo "ERROR: runmitos was not found after activating conda env mitos2." >&2
    echo "CONDA_PREFIX=${CONDA_PREFIX:-not_set}" >&2
    echo "PATH=$PATH" >&2
    exit 1
fi

echo "Using MITOS2 executable: $(command -v runmitos)"
runmitos --help >/dev/null
```

`run_mitos2_annotation.py` activates that environment in a login shell, validates `runmitos`, and records that executable in `results/qc/mitos2_annotation/mitos2_annotation_summary.tsv`. The conda environment name is `mitos2`, the installed package name is `mitos`, and the CLI executable name is `runmitos`.

The workflow runs one MITOS2 task per target species. Its primary input is the exact variant-calling FASTA, `references/variant_calling/Ref_chrM/{target_species}.fa`, whose standardized record header is `>chrM`. This is the coordinate truth for all emitted positions, even when `final_chrM_species` and the manifest accession identify a cross-species nearest reference. The manifest's `chrM_expected_output_fasta` is used only as a fallback when the target-species FASTA is absent, and its header is sanitized to `>chrM`. Manifest coordinate/reference metadata and the actual MITOS2 input path are both recorded in the task and summary tables. MITOS2 supplies CDS, tRNA, and rRNA *intervals*; its tRNA/rRNA output does not provide secondary-structure stem/loop information and does not replace tRNAscan paired-site annotations or human-guided rRNA stem/loop annotation. The interval table may support future fallback region tables.

The first integration target is CDS/codon fallback. `build_primate_codon_table.py` selects a single annotation source for each sample: valid nonzero GenBank CDS rows first, otherwise MITOS2 rows for that final reference. It never combines sources within a sample. Raw/parsed MITOS2 tables remain separate, and the source comparison is written to `results/qc/codon_table_build/genbank_vs_mitos2_cds_comparison.tsv`. The chosen table used downstream is `data/reference_tables/all_primate_position_codon_table.tsv`.

```bash
bash qc_analysis/scripts/run_qc_preprocessing.sh mitos2_annotation config/qc_preprocessing.yaml
bash qc_analysis/scripts/run_qc_preprocessing.sh build_primate_codon_table config/qc_preprocessing.yaml
bash qc_analysis/scripts/run_qc_preprocessing.sh all config/qc_preprocessing.yaml
```

## One-reference smoke test

Run the validated MITOS2 command path for one target-species variant-calling FASTA:

```bash
python qc_analysis/scripts/run_mitos2_annotation.py \
  --config config/qc_preprocessing.yaml \
  --reference Tarsius_lariang \
  --force
```

This writes `results/qc/mitos2_annotation/all_mitos2_features.tsv`,
`results/qc/mitos2_annotation/all_mitos2_position_codon_table.tsv`, and
`results/qc/mitos2_annotation/mitos2_annotation_summary.tsv`. A completed
Tarsius run contains approximately 13 CDS intervals and 11,000 coding-position
rows.

## Slurm array workflow (recommended)

`runmitos` does not provide useful multithreading, so parallelize across target
species FASTAs instead. First create a stable, one-based task list, then
submit one array task for each data row and merge only after the array finishes:

```bash
bash qc_analysis/scripts/run_qc_preprocessing.sh mitos2_prepare_tasks config/qc_preprocessing.yaml
N=$(($(wc -l < results/qc/mitos2_annotation/mitos2_reference_tasks.tsv)-1))
sbatch --array=1-${N}%20 qc_analysis/scripts/run_mitos2_annotation_array.slurm
bash qc_analysis/scripts/run_qc_preprocessing.sh mitos2_merge config/qc_preprocessing.yaml
```

The task list has one row per target-species variant-calling FASTA and records its task ID,
manifest reference metadata, actual MITOS2 input FASTA, number of linked samples, and completion status. Each array
worker writes only `results/qc/mitos2_annotation/raw/{reference_key}/`; successful
workers create `mitos2.completed.ok`, which causes later runs to skip that
reference unless `--force` is supplied. The merge command parses every completed
raw directory and regenerates the three combined tables.

`%20` limits concurrent workers to 20. Adjust this array concurrency to roughly
10–30 according to current cluster load and local scheduler policy.
