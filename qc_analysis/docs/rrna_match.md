# rRNA match
`run_rrna_match.py` prefers tRNA output, then codon output, then lifted raw VCFs. Human interval tables require `chrom`, `start`, `end`, `rrna_gene`, `strand`; species tables additionally require `sample` or `species`. It normalizes `12S`/`RNR1` to `MT-RNR1` and `16S`/`RNR2` to `MT-RNR2`.

It first performs interval-level rRNA gene matching and adds `MTRRNA_*` INFO annotations for source/human genes, local positions, lengths, fractions, fraction delta, strand, gene and region matches. Statuses are `OK`, `NO_SPECIES_RRNA`, `NO_HUMAN_RRNA`, `NO_SPECIES_OR_HUMAN_RRNA`, `GENE_MISMATCH`, and `MISSING_COORD`.

When `use_rrna_structure_table` is enabled, the required human structure table (`rrna_gene`, `human_pos`, `local_pos`, `struct_class`) adds **human-reference-guided** loop/stem fields, paired human positions, and the effect of the lifted VCF ALT on the human pair. For stem positions, the script infers the source paired genomic position from the human paired local coordinate, lifts it through the per-sample coordinate map, and assigns `HIGH_CONF_STEM` only when that position agrees. A missing map yields `MODERATE_CONF_STEM` rather than a failure. This is not a prediction of species-specific rRNA secondary structure.

If structural annotation is disabled, or the table does not contain the lifted human position, interval annotations remain available and structural fields are `.`/`NA`. Enabling it with a missing table fails clearly.

```bash
python qc_analysis/scripts/run_rrna_match.py --config config/qc_preprocessing.yaml --sample ERS14600320
```
