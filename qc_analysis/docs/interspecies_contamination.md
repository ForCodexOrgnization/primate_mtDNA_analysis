# Interspecies contamination QC

This report-only cohort step runs immediately after coordinate liftover and before
MITOS2/codon preparation. It reads only
`results/qc/coordinate_liftover/vcf_lifted_raw/{sample}.lifted.raw.vcf[.gz]` and
never rewrites a VCF or removes a sample. The terminal `final_filter` continues to
consume the report's `interspecies_status` column.

## Evidence model

Only biallelic canonical SNVs with `FILTER=PASS`, `DP >= 100`, and a usable AF
are indexed. Recipient alleles have `0.01 <= AF <= 0.20`; potential source
alleles have `AF >= 0.99`. Identity is the exact post-liftover
`CHROM/POS/REF/ALT` tuple. Source/original-coordinate INFO fields are ignored.

Before cross-species evaluation, a recipient's low-VAF allele is removed when
another cohort member of the recipient species carries it homoplasmically. The
implementation queries an allele-to-high-VAF-samples inverted index and aggregates
matches by source species and source sample; it performs neither all-pairs sample
intersections nor phylogenetic-tree analysis.

The report includes the best source, overlap count and fraction (denominator: the
low-VAF set after background removal), median matched recipient VAF, and the
fraction within the configured tolerance of that median. A coherent, unambiguous
signal meeting overlap thresholds is `FAIL`. Incoherent, ambiguous, and singleton
recipient-species signals are `WARN`; other samples are `PASS`.
`interspecies_status` is always exactly `PASS`, `WARN`, or `FAIL`, while
`classification` and `reason` provide detail.

## Run

```bash
python qc_analysis/scripts/run_interspecies_contamination.py --config config/qc_preprocessing.yaml
bash qc_analysis/scripts/run_qc_preprocessing.sh --submit interspecies_contamination config/qc_preprocessing.yaml
```

Use `--overwrite` to replace an existing report. Output is
`results/qc/interspecies_contamination/reports/interspecies_contamination_report.tsv`.
