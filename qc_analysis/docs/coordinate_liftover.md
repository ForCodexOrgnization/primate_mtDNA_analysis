# Coordinate liftover

Coordinate liftover rotates circular species and human mtDNA references, aligns each species to human with pairwise MAFFT, builds a base-level coordinate map, and then lifts VCF and coverage records. The final pairwise species–human alignment remains responsible for the production coordinate map.

Core principle:

> Global MSA is used to identify a common homologous circular cut point. Pairwise species–human alignment is used to generate the final coordinate map.

## Anchor sources

Runtime anchor priority is:

1. `SAMPLE_OVERRIDE` from a sample-specific `rotate_anchor`;
2. validated `GLOBAL_MSA_ANCHOR` from `reference_anchor_positions.tsv`;
3. validated `FAMILY_MSA_ANCHOR` from the same table when available;
4. `PAIRWISE_FALLBACK`, only when explicitly enabled;
5. failure with manual-review QC.

The global anchor determines only how circular references are cut before pairwise MAFFT. It does not replace the final pairwise alignment.

## Reference-level anchor reuse

Anchors are keyed by stable `reference_id` and `sequence_sha256`, not by sample name. Samples sharing one identical mtDNA reference sequence reuse the same species anchor. Liftover still runs per sample because VCF and coverage files are sample-specific.

## Sequence verification

When `verify_sequence_sha256: true`, runtime liftover hashes the current FASTA sequence and rejects stored anchors if the FASTA version differs. Strict failure reasons include:

- `ANCHOR_REFERENCE_HASH_MISMATCH`
- `ANCHOR_REFERENCE_LENGTH_MISMATCH`
- `ANCHOR_REFERENCE_ID_COLLISION`
- `ANCHOR_NOT_FOUND`
- `ANCHOR_POSITION_OUT_OF_RANGE`

If pairwise fallback is disabled, the affected sample is skipped and reported rather than using an unsafe anchor.

## Metadata

The sample manifest supports the minimal historical format:

```text
sample  species
```

It can also include optional columns: `family`, `reference_id`, `species_fasta`, `vcf`, `cov`, `species_chrom`, `target_sequence`, and `rotate_anchor`.

## Input files

The source VCF may be either an uncompressed `.vcf` file or a compressed
`.vcf.gz` file. When both configured candidates are present, configured pattern
order defines the selection priority; the standard configuration lists `.vcf.gz`
first, so it is preferred. A broken `.vcf.gz` symlink is reported in diagnostics
but does not prevent selection of a valid uncompressed `.vcf` fallback.

Coordinate liftover requires exactly the collection-step merged maximum-depth
coverage input named `{sample}.merged.max_depth.per_base_coverage.tsv`. The
original `{sample}.round2.original_coords.per_base_coverage.tsv` file is not
used directly.

## Running with preprocessing

Run global anchor discovery before coordinate liftover with the preprocessing wrapper:

```bash
bash qc_analysis/scripts/run_qc_preprocessing.sh discover_global_anchor config/qc_preprocessing.yaml
```

For an end-to-end preprocessing run, `all` executes collection, global-anchor discovery, and coordinate liftover in order:

```bash
bash qc_analysis/scripts/run_qc_preprocessing.sh all config/qc_preprocessing.yaml
```

## Configuration

Use `coordinate_liftover.coordinates.anchor_positions_file` to point at the global discovery table. Recommended production anchor settings are:

```yaml
coordinate_liftover:
  coordinates:
    anchor_positions_file: results/qc/global_anchor/reference_anchor_positions.tsv
  anchor:
    require_validated_anchor: true
    verify_sequence_sha256: true
    allow_pairwise_anchor_fallback: false
    allow_anchor_position_one_fallback: false
```

## QC reports

Per-sample QC includes `reference_id`, `reference_sequence_sha256`, `anchor_method`, `species_anchor_position`, `human_anchor_position`, `anchor_alignment_column`, `anchor_qc_status`, and `pairwise_anchor_fallback_used`. The cohort summary reports counts of samples using global, family, or pairwise anchors and counts of failed anchor validation.

## IUPAC ambiguity handling

FASTA validation accepts the complete standard IUPAC DNA alphabet, case-insensitively: `ACGTRYSWKMBDHVN`. The normalized original sequence is retained for reference identity, length, SHA256 validation, and source-coordinate tracking; hashes are never calculated from an alignment-masked copy.

Only temporary MAFFT inputs are masked: every ambiguous symbol (`RYSWKMBDHVN`) becomes `N`, preserving sequence length. Exact shared k-mer rotation anchors require `A`, `C`, `G`, or `T` at every position. Per-sample QC reports the number of ambiguous reference positions and variants overlapping them. Such variants are written to the unresolved report as `SOURCE_REFERENCE_AMBIGUOUS` (including the original FASTA base), while other records in the sample continue through liftover.

## Downstream annotation handoff

Coordinate liftover remains coordinate-only. Its raw lifted VCF handoff files in
`results/qc/coordinate_liftover/vcf_lifted_raw` are consumed by `codon_match`,
`tRNA_match`, and `rRNA_match`; maps in `maps` retain source-to-human relationships.
