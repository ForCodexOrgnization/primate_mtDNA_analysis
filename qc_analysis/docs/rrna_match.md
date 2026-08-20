# rRNA match
`run_rrna_match.py` prefers tRNA output, then codon output, then lifted raw VCFs. Human interval tables require `chrom`, `start`, `end`, `rrna_gene`, `strand`; species tables additionally require `sample` or `species`. It normalizes `12S`/`RNR1` to `MT-RNR1` and `16S`/`RNR2` to `MT-RNR2`.

It first performs interval-level rRNA gene matching and adds `MTRRNA_*` INFO annotations for source/human genes, local positions, lengths, fractions, fraction delta, strand, gene and region matches. Statuses are `OK`, `NO_SPECIES_RRNA`, `NO_HUMAN_RRNA`, `NO_SPECIES_OR_HUMAN_RRNA`, `GENE_MISMATCH`, and `MISSING_COORD`.

When `use_rrna_structure_table` is enabled, the matcher reads two independent
structure sources:

* `human_rrna_structure_table`, keyed by `rrna_gene` and `human_pos`
* `species_rrna_structure_table`, keyed by `reference_key`, `rrna_gene`, and
  `genomic_pos`

The species-side table is reference-level, not sample-level. Samples resolve to
that table through `sample_reference_map`, which must identify the exact
coordinate reference FASTA and sequence SHA256 used for variant calling. Species
names are not used as structure identities.

For each rRNA variant, the script annotates independent human and species
classes (`MTRRNA_H_CLASS`, `MTRRNA_S_CLASS`) and paired-site fields
(`MTRRNA_H_PAIR_POS`, `MTRRNA_S_PAIR_POS`, `MTRRNA_H_PAIR_TYPE`,
`MTRRNA_S_PAIR_TYPE`, and pair states). `MTRRNA_STRUCTURE_MATCH` is one of
`STEM_STEM`, `LOOP_LOOP`, `STEM_LOOP`, `LOOP_STEM`, or `UNKNOWN`.

For stem-stem positions, `MTRRNA_PAIR_RELATION_MATCH` compares pairing
relationships across genomes: the species position's independently annotated
species partner is lifted through the existing coordinate map and compared with
the independently annotated human partner. The old human-projected expected pair
logic is not used. `MTRRNA_S_PAIR_EXPECTED_POS` remains as a deprecated
compatibility INFO field and is no longer populated.

`MTRRNA_MATCH_TIER` now reflects two-sided structure: `HIGH_CONF_STEM` requires
`STEM_STEM` plus `MTRRNA_PAIR_RELATION_MATCH=yes`; `HIGH_CONF_LOOP` requires
`LOOP_LOOP`; `STRUCTURE_DISCORDANT` covers `STEM_LOOP` and `LOOP_STEM`; and
`STRUCTURE_UNKNOWN` is used when either side is missing or unresolved. A
compensatory pair-type change such as human `G-C` versus species `A-U` remains
`STEM_STEM` and is not structural discordance.

If structural annotation is disabled, interval annotations remain available and
structural fields are `.`/`NA`. Enabling it with a missing human or species
structure table fails clearly.

```bash
python qc_analysis/scripts/run_rrna_match.py --config config/qc_preprocessing.yaml --sample ERS14600320
```

After a Slurm-array `rrna_match`, merge the one-row per-sample summaries and
inspect the deterministic cohort summary:

```bash
bash qc_analysis/scripts/run_qc_preprocessing.sh \
  rrna_match_merge config/qc_preprocessing.yaml

less results/qc/rrna_match/reports/all_samples.rrna_match_summary.tsv
```

With `--submit rrna_match`, the singleton merge is submitted automatically
with an `afterok` dependency unless `AUTO_SUBMIT_MERGE=false`.
