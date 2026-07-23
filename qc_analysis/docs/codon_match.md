# Codon match
`run_codon_match.py` annotates, but never filters, raw lifted VCF records. It reads the reference-level primate codon table and the sample-reference map, resolves `sample -> reference_key`, then finds annotations by `reference_key + pos`. This prevents sample-level duplication while ensuring that coordinate reference identity—not species alone—selects the codon annotation. The human table requires `pos`, `gene`, `strand`, codon index/phase, codon sequence, and three genomic codon positions. INFO annotations are `MTCODON_STATUS`, match/strict/gene/phase flags, and source/human gene, codon, and phase values. Statuses are `PASS`, `SKIPPED_NONCODING`, `NO_HUMAN_CODON`, `GENE_MISMATCH`, `PHASE_MISMATCH`, `MISMATCH`, and `MISSING_COORD`.

Outputs are `vcf_codon/{sample}.lifted.codon.vcf` plus per-sample and merged report TSVs. It recognizes both `SRC_*` and `MTLIFT_ORIG_*` INFO coordinate conventions.

```bash
python qc_analysis/scripts/run_codon_match.py --config config/qc_preprocessing.yaml --sample ERS14600320
```
