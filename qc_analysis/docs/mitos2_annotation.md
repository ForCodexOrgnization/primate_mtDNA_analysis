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

The workflow runs one MITOS2 task per unique normalized sequence SHA256. Its primary input is the exact variant-calling FASTA, `references/variant_calling/Ref_chrM/{target_species}.fa`. Identical sequences shared by multiple species use one task and one biological `reference_key`; different sequences never collapse because of a shared species or accession. Accession remains optional provenance. MITOS2 supplies final CDS, tRNA, and rRNA intervals. CDS drives production codon matching. MITOS2 tRNA secondary structure is exported only for comparison/QC; production tRNA matching remains exclusively based on tRNAscan-SE and the existing reference position indexes.

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
`results/qc/mitos2_annotation/all_mitos2_reference_trna_structure.tsv`,
`results/qc/mitos2_annotation/all_mitos2_reference_rrna_regions.tsv`,
`results/qc/mitos2_annotation/all_mitos2_reference_rrna_structure.tsv`,
`results/qc/mitos2_annotation/codon_sample_reference_map.tsv`, and
`results/qc/mitos2_annotation/mitos2_annotation_summary.tsv` (plus the compact
feature diagnostic table). The merge streams one reference at a time and does
not generate the legacy sample-expanded codon table. A completed
Tarsius run contains approximately 13 CDS intervals and 11,000 coding-position
rows.

## rRNA structure table

`build_mitos2_rrna_structure_table.py` is used by the merge path and can also be
run directly. Production parsing reads final tRNA/rRNA records from
`<raw-task-directory>/result.mitos`. This is MITOS's custom tabular format, not
BED or GFF: core columns identify the feature, 0-based inclusive start/end,
numeric strand, and score, while the RNA-specific tail contains dot-bracket
structure. The parser explicitly normalizes MITOS start/end by `+1/+1` and only
uses a structure if feature type, gene, strand, and the complete normalized
interval equal exactly one final `result.gff` feature. A disagreement is reported
as `result_mitos_gff_interval_mismatch`; coordinates are never silently shifted.

Dot-bracket pairs are expanded reciprocally in RNA 5'-to-3' orientation. For a
negative-strand feature this means local position 1 maps to the high genomic
coordinate. All genomic bases and paired bases are read from the exact coordinate
reference FASTA. Parentheses, square/curly brackets, and angle brackets are
supported; dots are unpaired. Structure length must equal the final GFF feature
length. The MITOS2 command still uses `--best --noplots`, and SVG geometry is not
parsed for production annotations. Retained Stockholm parsing is legacy-only and
cannot override a present final `result.mitos`.

The rRNA table is reference-level and includes:
`reference_key`, `reference_species`, `coordinate_reference_accession`,
`coordinate_reference_fasta`, `coordinate_reference_sequence_sha256`,
`rrna_gene`, `genomic_pos`, `local_pos`, `base`, `struct_class`,
`paired_genomic_pos`, `paired_local_pos`, `paired_base`, `pair_type`,
`pair_state`, `annotation_source`, and `structure_source`, plus optional model
and strand diagnostics. rRNA structure extraction status is reported separately
from CDS production QC and does not invalidate otherwise valid codon outputs.

The merge also writes `all_mitos2_reference_rrna_regions.tsv`, one row per final
MITOS2 GFF rRNA feature, and `all_mitos2_reference_trna_structure.tsv`. Distinct
MITOS names such as `trnL1`/`trnL2` and `trnS1`/`trnS2` are preserved. RNA
failures are reported independently in `mitos2_annotation_summary.tsv` and never
remove an otherwise valid reference from the production codon table.

## rRNA matching coordinate interface

Primate rRNA lookup is reference-aware: a sample is resolved through
`sample_coordinate_reference_map.tsv` to its `reference_key`, then looked up in
the MITOS2 region and structure tables using the original/source variant
coordinate. Human rRNA interval and structure annotation remain the curated
tables and use the lifted human coordinate. Only after both sides are annotated
independently is a primate pairing partner lifted to human coordinates for
partner-relation comparison. Human local pairing coordinates are never projected
back onto the primate sequence.

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
