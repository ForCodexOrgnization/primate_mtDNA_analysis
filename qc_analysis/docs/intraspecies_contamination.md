# Intra-species contamination analysis

This independent pre-liftover QC module detects evidence consistent with
intra-species sample mixture or contamination. It does not prove the biological
or technical source of the mixture and it does not remove individual variants
from the VCF. Coordinates remain in the species reference; neither human
liftover, human-contamination filtering, nor inter-species filtering is done.

## Inputs and method

The builder recursively reads collected `.vcf`/`.vcf.gz` files and sample/species
metadata (a `sample species` table or `config/sample_ref_file.tsv`). It requires
one VCF sample column and uses only sample `FORMAT/AF`, `FORMAT/DP`, and (when
available) `FORMAT/AD`; INFO/DP is never substituted. The core table requires
`Sample, Species, CHROM, POS, REF, ALT, Type, FILTER, DP`, and `VAF` (or `AF`).

Comparisons are strictly within species to avoid interpreting fixed
between-species reference differences as mixture. For each tested A, low variants
are 0.01–0.20 and source B anchors are >=0.99. The best ordered A→B overlap is
selected. Leave-one-out anchors are high variants from **all other** samples of
the species. The mt-high-hets estimate is `1 - mean(VAF)`: use 0.80–0.998 when
at least three anchors are depressed, otherwise 0.80–1.00 when at least one is
available.

Defaults are DP >=100 and SNVs only. Candidate calls require >=5 low-A variants,
overlap >=3, fraction >=0.50, and estimate >=0.036420574377757434.
High-confidence calls use fraction >=0.6213636363636358 and estimate
>=0.07103935483870959. Singleton species are retained as
`insufficient_singleton_species`; samples with no usable variants and samples
without anchors are retained with insufficient-data statuses.

Mirror patterns are supporting evidence only and never set either contamination
flag. If negative-control pairs are absent, calibration is `not_calibrated_no_file`
and normalized support is `NA`, not false evidence. When supplied, the file must
have `Sample_A`, `Sample_B`, and `negative_control_tier`; tier 2 calibration uses
`tier2_location_and_batch_different`. The default mirror tolerance is zero,
meaning an exact complementary VAF sum, retained for compatibility.

## Configuration and usage

Set `intraspecies_contamination.enabled: true`. Choose exactly one mode:

```yaml
intraspecies_contamination:
  enabled: true
  build_variant_table: true
  vcf_dir: results/qc/collected_variant_calling_results/collected_vcf
  metadata: config/sample_ref_file.tsv
  outdir: results/qc/intraspecies_contamination
```

Or set `build_variant_table: false` and provide `variant_table`. Run locally:

```bash
bash qc_analysis/scripts/run_intraspecies_contamination.sh config/qc_preprocessing.yaml
```

Run on Slurm:

```bash
sbatch qc_analysis/scripts/submit_intraspecies_contamination.slurm config/qc_preprocessing.yaml
```

The main results are under `tables/`, six diagnostic plot placeholders/plots are
under `plots/`, and reproducibility logs and `run_parameters.tsv` record paths,
thresholds, timestamp, and R version. Runtime scales approximately with the
number of within-species sample pairs and variants. The Python driver reads configuration with the repository's standard-library-only
restricted YAML parser, consistent with coordinate liftover and codon match. PyYAML
is not required. The R analysis uses base R; all configured thresholds are passed
explicitly by the Python driver.
