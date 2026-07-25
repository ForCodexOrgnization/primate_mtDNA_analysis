# Codon match
`run_codon_match.py` annotates, but never filters, raw lifted VCF records. It reads the reference-level primate codon table and the sample-reference map, resolves `sample -> reference_key`, then finds annotations by `reference_key + pos`. This prevents sample-level duplication while ensuring that coordinate reference identity—not species alone—selects the codon annotation. Plain TSV and gzip-compressed TSV codon tables are streamed while their position indexes are built.

Codon-table positions are one-to-many: every CDS annotation at a position is preserved. In particular, the mitochondrial `MT-ATP8`/`MT-ATP6` and `MT-ND4L`/`MT-ND4` overlaps are not collapsed. Matching evaluates every source-human annotation combination, so a valid match in one overlapping CDS cannot be hidden by another annotation. The legacy single-valued source/human gene, codon, and phase INFO fields contain a deterministic best representative. `MTCODON_N_PRIMATE_ANNOTATIONS`, `MTCODON_N_HUMAN_ANNOTATIONS`, `MTCODON_N_PAIR_CANDIDATES`, `MTCODON_OVERLAPPING_CDS`, `MTCODON_PRIMATE_GENES`, `MTCODON_HUMAN_GENES`, `MTCODON_MATCHING_GENES`, and `MTCODON_AMBIGUOUS_BEST_MATCH` expose the complete overlap context.

The human table requires `pos`, `gene`, `strand`, codon index/phase, codon sequence, and three genomic codon positions. Statuses remain `PASS`, `SKIPPED_NONCODING`, `NO_HUMAN_CODON`, `GENE_MISMATCH`, `PHASE_MISMATCH`, `MISMATCH`, and `MISSING_COORD`. Summary reports additionally count records with overlapping source CDSs, overlapping human CDSs, either type of overlap, ambiguous best matches, and multiple pair candidates.

Outputs are `vcf_codon/{sample}.lifted.codon.vcf` plus per-sample and merged report TSVs. It recognizes both `SRC_*` and `MTLIFT_ORIG_*` INFO coordinate conventions.

```bash
python qc_analysis/scripts/run_codon_match.py --config config/qc_preprocessing.yaml --sample ERS14600320
```
