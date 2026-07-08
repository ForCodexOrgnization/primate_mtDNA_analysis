# Coordinate liftover QC module

This self-contained module lifts raw primate mitochondrial VCF and coverage
coordinates to canonical human chrM coordinates. It starts from raw, unrotated
species chrM FASTA files and a raw, unrotated human chrM FASTA; it does **not**
require pre-rotated FASTAs, precomputed MAFFT/MUSCLE/PRANK alignments, position
maps, or a mandatory `rotate_pos_file`.

## Inputs

Provide only:

1. raw species chrM FASTA files
2. raw human chrM FASTA
3. original species VCF files compressed as `.vcf.gz`
4. original species COV files
5. a `sample_ref_file` TSV or a sample list in the config

A minimal `sample_ref_file` should include only these columns:

```text
sample	species
```

Headerless two-column TSV files are also accepted, with column 1 interpreted as
`sample` and column 2 interpreted as `species`. Extra columns in headerless files
are ignored.

With the minimal format, the workflow resolves files from configured input
directories:

- species FASTA: `{species_fasta_dir}/{species}.fa` by default, with additional
  extensions from `species_fasta_extensions`
- VCF: the unique gzipped VCF in `vcf_dir` matching `vcf_pattern`, default
  `{sample}*.vcf.gz` (for example,
  `ERS12091861.round2.original_coords.clean.final.split.vcf.gz`)
- COV: the unique file in `cov_dir` matching `cov_pattern`, default
  `{sample}*.tsv,{sample}*.cov` (for example,
  `SAMN01920507.round1_round2.max_stage_coverage.tsv`)

Backward-compatible explicit columns are also accepted: `species_fasta`, `vcf`,
and `cov`. Optional columns are `species_chrom`, `rotate_anchor`, and
`target_sequence`. Alternatively, set `coordinate_liftover.paths.samples: sampleA,sampleB` and
add entries under `coordinate_liftover.samples_by_name` with either `species` or
explicit `species_fasta`, `vcf`, and `cov` values.

## Run

```bash
python qc_analysis/01_coordinate_liftover/run_coordinate_liftover.py \
  --config config/qc_preprocessing.yaml
```

Single-sample mode:

```bash
python qc_analysis/01_coordinate_liftover/run_coordinate_liftover.py \
  --config config/qc_preprocessing.yaml \
  --sample SAMPLE_NAME
```

## Slurm submission

Use the provided Slurm submission script to run the workflow on the cluster:

```bash
sbatch qc_analysis/01_coordinate_liftover/submit_coordinate_liftover.slurm
```

For one sample, pass `SAMPLE` with `--export`:

```bash
sbatch --export=ALL,SAMPLE=SAMPLE_NAME \
  qc_analysis/01_coordinate_liftover/submit_coordinate_liftover.slurm
```

For an array run, submit one task per non-header row in `config/sample_ref_file.tsv`:

```bash
N=$(awk 'BEGIN{FS="\t"} $0 !~ /^[[:space:]]*#/ && NF >= 2 && tolower($1) != "sample" {n++} END{print n}' config/sample_ref_file.tsv)
sbatch --array=1-${N}%20 qc_analysis/01_coordinate_liftover/submit_coordinate_liftover.slurm
```

Override defaults with `--export`, for example `CONFIG=...`, `SAMPLE_REF=...`,
`PYTHON=...`, or `REPO_DIR=...`.

## MAFFT environment

By default, `coordinate_liftover.alignment.use_conda_env: true`, so the alignment step first tries
to load MAFFT through the configured HPC conda module. Configure `module_load`
and `conda_env` under `coordinate_liftover` in `config/qc_preprocessing.yaml` for your cluster.

When `allow_simple_alignment_fallback = true`, missing module, `conda`, or
`mafft` commands fall back to the deterministic simple alignment path without
printing noisy shell setup errors.

## Outputs

The module creates these subdirectories under `coordinate_liftover.paths.output_dir`:

- `prepared_fastas/` normalized single-record FASTAs
- `rotated_fastas/` rotated species and human FASTAs
- `alignments/` pairwise species-to-human FASTA alignments
- `maps/` coordinate maps with original, rotated, and canonical coordinates
- `vcf_lifted_raw/` raw lifted VCFs with source coordinates preserved in INFO
- `cov_lifted/` lifted coverage files with source coordinates preserved
- `reports/` per-sample QC reports and `all_samples.coordinate_liftover_summary.tsv`

Coordinate maps include `sample`, `species_chrom`, `species_pos_original`,
`species_pos_rotated`, `human_pos_rotated`, `human_pos_canonical`,
`species_base`, `human_base`, and `map_status`.

## Scope

This module only produces coordinate-lifted raw VCF/COV files and liftover QC
outputs. It intentionally does not perform codon annotation, tRNA annotation,
codon matching, tRNA matching, or final filtering. Downstream codon/tRNA matching
should consume this module's lifted raw outputs.
