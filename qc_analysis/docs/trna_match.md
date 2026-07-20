# tRNA match
`run_trna_match.py` prefers codon-annotated VCFs and falls back to lifted raw VCFs. It is annotation-only and adds `MTTRNA_*` INFO fields for IDs, local positions, class/element, pairing, ALT effect, compensation, strict match, and source lookup coordinates. Statuses are `OK`, `NO_SPECIES_TRNA`, `NO_HUMAN_TRNA`, `NO_SPECIES_OR_HUMAN_TRNA`, and `MISSING_SPECIES_COORD`.

Position indexes are TSV/TSV.GZ files with the configured tRNAscan-derived columns, including `chrom`, `pos`, tRNA interval/identity, structural and pairing fields. Stem variants additionally check source paired-site liftover, pair state, ALT pair effect, and whether both ALT pairs remain compatible (`WC` or `GU_wobble`). The source-side ALT effect uses original `SRC_ALT`; the human-side effect uses the current lifted VCF ALT. Loop strict matching instead requires matching loop class, element, region, and local position. Stem strict matching requires all structural comparisons, including paired-position and allele-effect agreement; `require_compensated_for_strict_stem` controls whether compatibility is also required (default true).

Missing indexes fail clearly unless index generation inputs are configured. A missing coordinate map does not fail annotation, but leaves paired-site liftover unavailable and prevents strict stem matching. Reports include a per-sample summary, merged summary, and tRNA gene-liftover QC placeholder/report. tRNAscan executable/environment settings live in the YAML config.

```bash
python qc_analysis/scripts/run_trna_match.py --config config/qc_preprocessing.yaml --sample ERS14600320
```
