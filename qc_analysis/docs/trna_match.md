# tRNA match
`run_trna_match.py` prefers codon-annotated VCFs and falls back to lifted raw VCFs. It is annotation-only and adds `MTTRNA_*` INFO fields for IDs, local positions, class/element, pairing, ALT effect, compensation, strict match, and source lookup coordinates. Statuses are `OK`, `NO_SPECIES_TRNA`, `NO_HUMAN_TRNA`, `NO_SPECIES_OR_HUMAN_TRNA`, and `MISSING_SPECIES_COORD`.

Position indexes are TSV/TSV.GZ files requiring `chrom`, `pos`, `trna_id`,
`local_pos`, `struct_class`, `struct_element`, `pair_type`, `pair_state`,
`paired_local_pos`, `paired_genomic_pos`, `paired_base`, and `strand`. Index
`paired_base` and VCF REF/ALT bases are **genomic orientation**; negative-strand
records are complemented and converted to RNA before `pair_type()` is called.
The stored `pair_type` describes the reference pair in transcript/RNA orientation
(pair class is unchanged when both genomic bases are complemented).

Chromosome normalization supports `none`, `strip_chr`, `add_chr`, and
`mitochondrial_alias` (the latter maps chrM/M/MT/mitochondrial names to `MT`).
Index and VCF normalization are configured independently. With
`*_lookup_ignore_chrom: false`, only normalized chromosome-plus-position matches
are accepted. When true, a position-only fallback is allowed, but ambiguous
positions are reported rather than arbitrarily selecting a record.

Stem variants additionally check source paired-site liftover, pair state, ALT effect, and whether both ALT pairs remain compatible (`WC` or `GU_wobble`). The source-side ALT effect uses original `SRC_ALT`; the human-side effect uses the current lifted VCF ALT. Loop strict matching instead requires matching loop class, element, region, and local position. Stem strict matching requires all structural comparisons, including paired-position and allele-effect agreement; `require_compensated_for_strict_stem` controls whether compatibility is also required (default true).

Missing indexes fail clearly. Index generation (`run_trnascan_if_missing`) and
interval gene QC (`run_trna_gene_qc`) are not implemented and are rejected when
enabled; their related executable settings therefore do not imply execution.
`pass_only` filters non-PASS VCF records and `write_summary` controls per-sample
and combined summaries. A missing coordinate map does not fail annotation, but
is counted and prevents strict stem matching. No gene-QC completion report is
written.

```bash
python qc_analysis/scripts/run_trna_match.py --config config/qc_preprocessing.yaml --sample ERS14600320
```
