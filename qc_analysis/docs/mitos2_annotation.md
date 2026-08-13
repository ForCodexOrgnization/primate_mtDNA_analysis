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

The workflow runs one MITOS2 task per unique normalized sequence SHA256. Its primary input is the exact variant-calling FASTA, `references/variant_calling/Ref_chrM/{target_species}.fa`. Identical sequences shared by multiple species use one task and one biological `reference_key`; different sequences never collapse because of a shared species or accession. Accession remains optional provenance. MITOS2 supplies CDS, tRNA, and rRNA *intervals*; its tRNA/rRNA output does not provide secondary-structure stem/loop information and does not replace tRNAscan paired-site annotations or human-guided rRNA stem/loop annotation.

Task preparation writes the annotation-independent `data/reference_tables/sample_coordinate_reference_map.tsv`, which maps every sample to the exact Ref_chrM used for variant calling regardless of MITOS2 status. The merge writes the production reference codon table and `codon_sample_reference_map.tsv` under `results/qc/mitos2_annotation/`; only references passing strict production QC enter those codon-specific outputs. Failed references remain explicit in the summary and remain eligible for the independent tRNAscan-SE workflow. There is no GenBank codon fallback. `build_primate_codon_table.py` is retained only to build an independent GenBank benchmark for the optional hash-matched comparison.

```bash
bash qc_analysis/scripts/run_qc_preprocessing.sh mitos2_annotation config/qc_preprocessing.yaml
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

This writes `results/qc/mitos2_annotation/all_mitos2_reference_position_codon_table.tsv`,
`results/qc/mitos2_annotation/codon_sample_reference_map.tsv`, and
`results/qc/mitos2_annotation/mitos2_annotation_summary.tsv` (plus the compact
feature diagnostic table). The merge streams one reference at a time and does
not generate the legacy sample-expanded codon table. A completed
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

The task list has one row per unique normalized variant-calling sequence and records its task ID,
manifest reference metadata, actual MITOS2 input FASTA, number of linked samples, and completion status. Each array
worker writes only `results/qc/mitos2_annotation/raw/{task_key}/`; successful
workers create `mitos2.completed.ok`, which causes later runs to skip that
reference unless `--force` is supplied. The merge command parses every completed
raw directory and regenerates the feature table, production reference codon table,
sample reference map, legacy diagnostic sample table, and QC summary.

`%20` limits concurrent workers to 20. Adjust this array concurrency to roughly
10–30 according to current cluster load and local scheduler policy.
