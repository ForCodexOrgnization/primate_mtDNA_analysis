# rRNA match
`run_rrna_match.py` prefers tRNA output, then codon output, then lifted raw VCFs. Human interval tables require `chrom`, `start`, `end`, `rrna_gene`, `strand`; species tables additionally require `sample` or `species`. It normalizes `12S`/`RNR1` to `MT-RNR1` and `16S`/`RNR2` to `MT-RNR2`.

It adds `MTRRNA_*` INFO annotations for source/human genes, local positions, lengths, fractions, fraction delta, strand, gene and region matches. Statuses are `OK`, `NO_SPECIES_RRNA`, `NO_HUMAN_RRNA`, `NO_SPECIES_OR_HUMAN_RRNA`, `GENE_MISMATCH`, and `MISSING_COORD`. Outputs are annotated VCFs and per-sample/merged summaries; secondary-structure settings are reserved for future use.

```bash
python qc_analysis/scripts/run_rrna_match.py --config config/qc_preprocessing.yaml --sample ERS14600320
```
