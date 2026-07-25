# Codon match
`run_codon_match.py` annotates, but never filters, raw lifted VCF records. It reads the reference-level primate codon table and the sample-reference map, resolves `sample -> reference_key`, then finds annotations by `reference_key + pos`. This prevents sample-level duplication while ensuring that coordinate reference identity—not species alone—selects the codon annotation. Plain TSV and gzip-compressed TSV codon tables are streamed while their position indexes are built.

Codon-table positions are one-to-many: every CDS annotation at a position is preserved. In particular, the mitochondrial `MT-ATP8`/`MT-ATP6` and `MT-ND4L`/`MT-ND4` overlaps are not collapsed. Matching evaluates every source-human annotation combination, so a valid match in one overlapping CDS cannot be hidden by another annotation. The legacy single-valued source/human gene, codon, and phase INFO fields contain a deterministic best representative. `MTCODON_N_PRIMATE_ANNOTATIONS`, `MTCODON_N_HUMAN_ANNOTATIONS`, `MTCODON_N_PAIR_CANDIDATES`, `MTCODON_OVERLAPPING_CDS`, `MTCODON_PRIMATE_GENES`, `MTCODON_HUMAN_GENES`, `MTCODON_MATCHING_GENES`, and `MTCODON_AMBIGUOUS_BEST_MATCH` expose the complete overlap context.

Alternate codons are constructed only when `SRC_REF` and `SRC_ALT` are each one of `A`, `C`, `G`, or `T` (case-insensitive), making the source variant a simple biallelic SNV. Insertions, deletions, multi-allelic ALT values, symbolic alleles, missing alleles, and ambiguous bases receive `MTCODON_SUPPORTED_SNV=no`, `MTCODON_STATUS=UNSUPPORTED_NON_SNV`, and no alternate codon or match. Coordinate, noncoding, and missing-human-annotation statuses take precedence, so their useful annotation context is retained.

Reference-level mode requires both a reference codon table and sample-reference map, and every requested sample must have a mapping. A missing mapping or duplicate sample rows with conflicting reference keys is a fatal configuration error; identical duplicate mappings are tolerated. Historical sample-level tables configured with `all_primate_position_codon_table` continue to use the sample name as their lookup key without requiring a map.

The human table requires `pos`, `gene`, `strand`, codon index/phase, codon sequence, and three genomic codon positions. Statuses are `PASS`, `SKIPPED_NONCODING`, `NO_HUMAN_CODON`, `UNSUPPORTED_NON_SNV`, `GENE_MISMATCH`, `PHASE_MISMATCH`, `MISMATCH`, and `MISSING_COORD`. Summary reports additionally count records with each status, overlapping source CDSs, overlapping human CDSs, either type of overlap, ambiguous best matches, and multiple pair candidates.

Outputs are `vcf_codon/{sample}.lifted.codon.vcf` plus per-sample and merged report TSVs. It recognizes both `SRC_*` and `MTLIFT_ORIG_*` INFO coordinate conventions.

```bash
python qc_analysis/scripts/run_codon_match.py --config config/qc_preprocessing.yaml --sample ERS14600320
```
