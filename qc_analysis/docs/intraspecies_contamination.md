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

Defaults are **DP >= 100** and SNVs only. This inclusive production boundary
matches the validated contamination R implementation; it intentionally differs
from descriptive plotting and other analyses that may use `DP > 100`. Candidate calls require >=5 low-A variants,
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

## Report-only contamination evidence score

The report also includes a **0-10 weighted evidence score** for ranking and
sensitivity review. It does **not** change the validated candidate/high-confidence
flags, `contamination_status`, or `qc_status`.

A score is reported only when the minimum evidence gate passes:
`n_lowA >= min_n_lowA`, `best_overlap >= min_overlap`, and a best same-species
source sample is available.

Weights are:

- donor/source matching: **4.0** total
  - overlap count: up to 1.5
  - overlap fraction: up to 2.5
- mt-high-hets: **2.0**
- genome-wide dispersion: **2.0**
  - occupied 1-kb bins: up to 0.75
  - circular span / configured mtDNA length: up to 0.75
  - max local fraction in a 1-kb window: up to 0.50
- AF coherence (overlap MAD): **1.5**; when overlap <5, this component is capped
  at 1.0 because very small overlap sets can appear artificially coherent
- mirror support: **0.5**
  - use calibrated normalized mirror support when p95/p99 controls are available
  - otherwise use raw mirror fraction as a provisional supporting metric

Initial report-only interpretations are:
`>=7 strong_evidence`, `5-6.99 candidate_evidence`,
`3-4.99 weak_ambiguous_evidence`, and `<3 little_evidence`.
These bins are exploratory and are not production FAIL/PASS cutoffs.

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

The production implementation is entirely Python and has no R runtime dependency.
Its sole biological result is
`reports/intraspecies_contamination_report.tsv`, with one row per collected sample;
`logs/` and `run_parameters.tsv` provide diagnostics and reproducibility metadata.
The stable report includes the low/high overlap, leave-one-out anchor and mirror
counts, categorical `contamination_status`, boolean candidate/high-confidence
flags, and machine-readable `qc_status` (`PASS`, `WARN`, or `FAIL`). Runtime scales
approximately with the number of within-species sample pairs and variants. The
driver uses the repository's standard-library-only restricted YAML parser, so
neither PyYAML nor R is required.

## Validated production semantics

The R program in `qc_analysis/validation/contamination_reference.R` is the
algorithmic validation truth and Python is the production implementation. Both
operate before liftover, in original species coordinates. In particular,
mirror pairs are a cross product of all low- and high-VAF variants in one sample;
alleles and positions need not match, and the tolerance boundary is inclusive.
`n_mirror_pairs` counts qualifying cross-product rows, while
`n_low_variants_with_mirror` counts distinct low rows with at least one match.
Tier-2 negative-control pairs establish type-7 p95/p99 thresholds for normalized
mirror evidence; unavailable/insufficient calibration is explicitly reported,
never replaced by raw-count thresholds.

The contamination and five-criterion sample-QC jobs only emit evidence reports.
The terminal `final_filter` requires both reports, treats candidate and
insufficient contamination calls as warnings, removes high-confidence calls or
five-criterion QC failures, and selects rRNA → tRNA → codon → raw-lifted VCFs.
