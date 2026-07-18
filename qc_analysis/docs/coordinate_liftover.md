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
directories. Set `vcf_dir` and `cov_dir` directly to the standardized collected
input directories produced by `collect_variant_calling.outdir`:

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
python qc_analysis/scripts/run_coordinate_liftover.py \
  --config config/qc_preprocessing.yaml
```

Single-sample mode:

```bash
python qc_analysis/scripts/run_coordinate_liftover.py \
  --config config/qc_preprocessing.yaml \
  --sample SAMPLE_NAME
```


## Run with collected variant-calling inputs

Use the QC preprocessing wrapper to submit collection, coordinate liftover, or both
steps sequentially through Slurm:

```bash
bash qc_analysis/scripts/run_qc_preprocessing.sh --submit collect_variant_calling_results config/qc_preprocessing.yaml
bash qc_analysis/scripts/run_qc_preprocessing.sh --submit coordinate_liftover config/qc_preprocessing.yaml
bash qc_analysis/scripts/run_qc_preprocessing.sh --submit all config/qc_preprocessing.yaml
```

Set `collect_variant_calling.outdir` in `config/qc_preprocessing.yaml`; the
default coordinate-liftover `vcf_dir` and `cov_dir` point at that collected
output layout.

## Slurm submission

Use the provided Slurm submission script to run the workflow on the cluster:

```bash
sbatch qc_analysis/scripts/submit_coordinate_liftover.slurm
```

For one sample, pass `SAMPLE` with `--export`:

```bash
sbatch --export=ALL,SAMPLE=SAMPLE_NAME \
  qc_analysis/scripts/submit_coordinate_liftover.slurm
```

For an array run, submit one task per non-header row in `config/sample_ref_file.tsv`:

```bash
N=$(awk 'BEGIN{FS="\t"} $0 !~ /^[[:space:]]*#/ && NF >= 2 && tolower($1) != "sample" {n++} END{print n}' config/sample_ref_file.tsv)
sbatch --array=1-${N}%20 qc_analysis/scripts/submit_coordinate_liftover.slurm
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

## Allele orientation and ALT/REF flips

Lifted VCF records are classified against the target human chrM reference before
being written:

- `REF_MATCH`: the target human REF equals the source REF. These records keep the
  original REF, ALT, QUAL, FILTER, INFO annotations, FORMAT keys, and sample
  values, except for lifted `CHROM`/`POS`, `SRC_*` provenance fields, and
  `LIFTOVER_ALLELE_STATUS=REF_MATCH`.
- `ALT_REF_FLIP`: the target human REF equals the source ALT and the record is a
  split, biallelic, non-symbolic SNV with one A/C/G/T base in both REF and ALT.
  The workflow swaps REF and ALT and rewrites allele-order-dependent fields.
- `FLIP_UNSUPPORTED_VARIANT_TYPE`: the target human REF equals the source ALT,
  but the variant is not currently safe to flip automatically, such as an indel
  or other non-simple SNV.
- `TARGET_REF_NOT_SOURCE_REF_OR_ALT`: the target human REF matches neither the
  source REF nor the source ALT.
- `UNMAPPED`: the source position has no mapped human coordinate.

Only split biallelic SNVs are flipped in this first implementation. Multiallelic
records, indels, symbolic alleles, `*` alleles, variants spanning alignment gaps,
and records where the human REF is neither the source REF nor source ALT are not
silently emitted with an incorrect REF. They are written to the unresolved report
instead.

Example:

```text
Species-oriented record:
G>A, AF=0.10, AD=90,10

Human-reference-oriented record:
A>G, AF=0.90, AD=10,90
```

`AF` changes to `1-AF` because after the flip the new ALT is the old REF. When a
valid two-value `AD` is available, the new AF is calculated from the original REF
depth as `old_ref_depth / (old_ref_depth + old_alt_depth)`; otherwise a valid
single AF in `[0,1]` is inverted.

For `ALT_REF_FLIP`, the workflow transforms allele-order-dependent FORMAT fields:

- `AD`, `FAD`, `F1R2`, and `F2R1` two-value arrays are swapped from `REF,ALT` to
  `ALT,REF`.
- `SB` is treated as Mutect-style `REF_FWD,REF_REV,ALT_FWD,ALT_REV` and becomes
  `ALT_FWD,ALT_REV,REF_FWD,REF_REV`.
- `GT` allele indexes are swapped (`0 -> 1`, `1 -> 0`) while preserving ploidy,
  missing alleles, and `/` versus `|` separators. The workflow does not normalize
  `1/0` back to `0/1`.
- Header-declared `Number=R` FORMAT fields with two values are swapped.
- Header-declared `Number=G` genotype-likelihood fields are permuted only for
  supported haploid/diploid biallelic shapes; unsupported values are cleared to
  `.` rather than retained with the wrong allele orientation.

INFO fields are handled conservatively. Header-declared `Number=R` INFO fields
with two values are swapped. Unknown `Number=A` and `Number=G` INFO annotations
are removed after a flip because they describe the old ALT allele and usually do
not contain a value for the old REF/new ALT; removed keys are listed in
`LIFTOVER_DROPPED_INFO_FIELDS`.

Before any successful record is written, the final REF allele is checked against
the human FASTA at the lifted position. A failed final check is reported as
`FINAL_REF_MISMATCH` and is not written to the lifted VCF, so the output VCF must
not contain records whose REF disagrees with the human chrM FASTA.

## Unresolved report

Each sample writes unresolved VCF records to:

```text
reports/<sample>.coordinate_liftover_unresolved.tsv
```

The report includes `sample`, source allele fields, target coordinate/reference
fields, and a standard `reason`, including `UNMAPPED`, `MULTIALLELIC`,
`SYMBOLIC_ALLELE`, `FLIP_UNSUPPORTED_VARIANT_TYPE`,
`TARGET_REF_NOT_SOURCE_REF_OR_ALT`, `FINAL_REF_MISMATCH`, and
`MALFORMED_RECORD`.

## Liftover configuration

The `coordinate_liftover.liftover` section supports:

```yaml
liftover:
  check_ref_against_human_fasta: true
  enable_ref_alt_flip: true
  fail_on_unresolvable_target_ref: false
  fail_on_ref_mismatch: false
```

`enable_ref_alt_flip` controls automatic simple biallelic SNV flips.
`fail_on_unresolvable_target_ref` stops the workflow when the target REF is
neither the source REF nor source ALT. `fail_on_ref_mismatch` is a deprecated
compatibility option used only as the fallback value for
`fail_on_unresolvable_target_ref`; it no longer causes resolvable
`ALT_REF_FLIP` records to be discarded.

## QC summary allele-orientation fields

`all_samples.coordinate_liftover_summary.tsv` and per-sample QC reports include:

- `ref_match_count`: source REF equals target human REF.
- `ref_mismatch_count`: source REF differs from target human REF, including both
  resolved flips and unresolved mismatches.
- `alt_ref_flip_count`: simple biallelic SNVs successfully fixed by ALT/REF flip.
- `unresolved_ref_mismatch_count`: target REF matched neither source REF nor ALT.
- `unsupported_flip_count`: target REF matched source ALT, but the variant type is
  not currently flipped automatically.
- `final_ref_mismatch_count`: transformed records that still failed the final REF
  check and were withheld from the lifted VCF.
