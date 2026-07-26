# Reference-level tRNA structural matching

## Data flow

```text
reference FASTA -> tRNAscan-SE -> .trnascan.out + .trnascan.ss
                -> reference-level position index
sample -> sample_reference_map.tsv -> reference_key -> position index
lifted/codon-annotated VCF + human/source indexes -> tRNA match VCF
```

Index construction depends only on the reference FASTA and tRNAscan-SE; it does
not depend on codon matching. Matching prefers a codon-annotated VCF and falls
back to the lifted raw VCF.

## Coordinate and orientation contract (format version 2)

`pos` and `paired_genomic_pos` are 1-based coordinates in the original input
reference. `base_genomic` and `paired_base_genomic` are DNA letters read from
that FASTA. The `.ss` sequence, `base_rna`, and `paired_base_rna` are in mature
5'-to-3' transcript orientation and use U. Thus plus-strand RNA is genomic DNA
with T changed to U; minus-strand RNA is the complement of genomic DNA with T
changed to U. VCF alleles are genomic and are converted exactly once. Pair type
is calculated only from the two RNA-orientation bases. In particular,
`paired_base_rna` must never be complemented again.

Every row contains:

`index_format_version, base_orientation, pair_type_orientation,
coordinate_space, reference_key, chrom, pos, trna_id, trna_begin, trna_end,
strand, aa, anticodon, score, local_pos, base_genomic, base_rna, struct_char,
struct_class, struct_element, paired_local_pos, paired_genomic_pos,
paired_base_genomic, paired_base_rna, pair_bases_rna, pair_type, pair_status,
pair_state, base, paired_base, fasta_sha256`.

The metadata values are `2`, `genomic_and_rna`, `transcript_rna`, and
`original_reference`. Compatibility aliases `base` and `paired_base` are
**genomic**. A legacy alias-only index must explicitly declare its orientation;
the matcher does not guess.

Pairing semantics are deliberately separate: `pair_status` is
`paired|unpaired`; `pair_state` is `WC|non_WC|NA`; and `pair_type` is
`WC|GU_wobble|non_WC|NA`. The SHA256 is calculated over the selected uppercase
mitochondrial sequence (not FASTA wrapping or headers).

Structural elements are inferred from topology, not fixed positions. The
builder groups consecutive antiparallel pairs into stems, identifies the
outermost terminal stem as the acceptor stem, assigns internal stems in
transcript order to D, anticodon, and T regions, labels the enclosed loops,
and labels the anticodon-to-T gap as the variable loop. Canonical mature-tRNA
position ranges are only a fallback when no pairing topology is present. This
also applies on the negative strand because local positions always follow the
mature 5'-to-3' RNA.

Exact duplicate index rows are collapsed. Distinct predictions at the same
chromosome/position are retained as lists and yield
`AMBIGUOUS_SPECIES_TRNA` or `AMBIGUOUS_HUMAN_TRNA`; the matcher never chooses
the last row. Summaries report overlapping and multi-tRNA position counts.

The existing strict stem rule requires region/element/pair-status agreement,
paired-coordinate agreement, allele-effect agreement, and (by default) ALT
pair compatibility. It does **not** require reference `PAIR_TYPE_MATCH` unless
`strict_stem_require_reference_pair_type_match: true` is configured. The
legacy `MTTRNA_COMPENSATED` name is retained for output compatibility, but it
means only that both compared ALT pairs are WC/GU-compatible—not that a
two-site compensatory mutation occurred.

The builder checks the transcript-oriented `.ss` sequence against the
strand-aware FASTA sequence. Its default mismatch threshold is zero and can be
changed with `--max-sequence-mismatch-rate`.

Reference manifests prefer explicit `chrM_fasta_path`/`chrM_expected_output_fasta`
fields and never fall back to `wg_expected_output_fasta`. A multi-record FASTA
is rejected unless `target_sequence_id` is explicit; in that case only that
record is passed to tRNAscan-SE. Configurable mitochondrial length bounds are
validated. Indexes are fully validated at format version 2 and atomically
renamed into place, so truncated interrupted outputs are never accepted.

## Commands

Build human from saved results:

```bash
python qc_analysis/scripts/build_trna_position_index.py --reference-key human \
  --fasta data/reference_tables/human_chrM.fa \
  --trnascan-out human.trnascan.out --trnascan-ss human.trnascan.ss \
  --output data/reference_tables/trna_index/human.trna_position_index.tsv.gz
```

Run tRNAscan and build one primate reference:

```bash
python qc_analysis/scripts/build_trna_position_index.py --reference-key REF \
  --fasta reference.fa --run-trnascan --trnascan-bin tRNAscan-SE \
  --trnascan-mode mito_mammal --threads 4 \
  --output data/reference_tables/trna_index/references/REF.trna_position_index.tsv.gz
```

Build all unique references (shared references run once):

```bash
python qc_analysis/scripts/build_all_trna_indexes.py \
  --config config/qc_preprocessing.yaml --workers 4
```

Run one sample:

```bash
python qc_analysis/scripts/run_trna_match.py --config config/qc_preprocessing.yaml --sample SAMPLE
```

Prepare the unique-reference Slurm manifest and submit its array:

```bash
python qc_analysis/scripts/build_all_trna_indexes.py --config config/qc_preprocessing.yaml \
  --task-manifest results/qc/trna_match/index_build_reports/trna_index_tasks.tsv
sbatch --array=1-N qc_analysis/scripts/slurm/build_trna_indexes_array.sbatch \
  results/qc/trna_match/index_build_reports/trna_index_tasks.tsv config/qc_preprocessing.yaml
```

The array activates `trnascan_env` and validates `tRNAscan-SE`. Standard
`.out`/`.ss` emitted by tRNAscan-SE v1/v2 are supported. Cove-only output,
pseudogene-specific alternate layouts, and structurally malformed or
sequence-less `.ss` records are rejected rather than inferred.
