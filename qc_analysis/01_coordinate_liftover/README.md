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
3. original species VCF files
4. original species COV files
5. a `sample_ref_file` TSV or a sample list in the config

A minimal `sample_ref_file` should include only these columns:

```text
sample	species
```

With the minimal format, the workflow resolves files from configured input
directories:

- species FASTA: `{species_fasta_dir}/{species}.fa` by default, with additional
  extensions from `species_fasta_extensions`
- VCF: the unique file in `vcf_dir` matching `vcf_pattern`, default
  `{sample}*.vcf` (for example,
  `ERS12091861.round2.original_coords.clean.final.split.vcf`)
- COV: the unique file in `cov_dir` matching `cov_pattern`, default
  `{sample}*.tsv,{sample}*.cov` (for example,
  `SAMN01920507.round1_round2.max_stage_coverage.tsv`)

Backward-compatible explicit columns are also accepted: `species_fasta`, `vcf`,
and `cov`. Optional columns are `species_chrom`, `rotate_anchor`, and
`target_sequence`. Alternatively, set `[paths] samples = sampleA,sampleB` and
create sections such as `[sample:sampleA]` with either `species` or explicit
`species_fasta`, `vcf`, and `cov` entries.

## Run

```bash
python qc_analysis/01_coordinate_liftover/run_coordinate_liftover.py \
  --config qc_analysis/01_coordinate_liftover/config.coordinate_liftover.ini
```

Single-sample mode:

```bash
python qc_analysis/01_coordinate_liftover/run_coordinate_liftover.py \
  --config qc_analysis/01_coordinate_liftover/config.coordinate_liftover.ini \
  --sample SAMPLE_NAME
```

## MAFFT environment

By default, the alignment step loads MAFFT through the configured HPC conda
environment before running `mafft`:

```bash
module load miniconda/24.11.3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mafft_env
```

These defaults are controlled by `[alignment] use_conda_env`, `module_load`,
and `conda_env` in `config.coordinate_liftover.ini`.

## Outputs

The module creates these subdirectories under `[paths] output_dir`:

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
