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

The report also includes a **0-1 weighted evidence score** for ranking and
sensitivity review. It does **not** change the validated candidate/high-confidence
flags, `contamination_status`, or `qc_status`.

A score is reported only when the minimum evidence gate passes:
`n_lowA >= min_n_lowA`, `best_overlap >= min_overlap`, and a best same-species
source sample is available.

Weights are:

- donor/source matching: **up to 4.0**
  - each best-source matched allele is first down-weighted by how common that
    high-A allele is among same-species source candidates:
    - <=5% frequency: weight 1.00
    - 5-10%: 0.75
    - 10-25%: 0.50
    - 25-50%: 0.25
    - >=50%: 0.00
  - the resulting background-adjusted overlap count and adjusted overlap fraction
    are converted to 0-1 strengths and combined by their geometric mean; they are
    not added as independent evidence
  - if donor background cannot be estimated because only one source candidate is
    available, the score falls back to raw overlap/fraction and reports this basis
  - the resulting source component is then adjusted by target-source provenance:
    - same cohort: x1.00
    - same project with cohort unknown: x0.90
    - different cohort (without a known project mismatch): x0.75
    - different project: x0.60
    - missing project/cohort metadata: x1.00 (neutral; missing metadata is not penalized)
  - provenance therefore changes only how plausible the nominated donor/source is;
    it does not independently create contamination evidence
- mt-high-hets: **2.5**
  - scored only when `mt_high_hets_mode == depressed_anchors` with at least
    3 depressed-high anchors
  - fallback estimates from fewer than 3 depressed anchors remain reported
    diagnostically but contribute 0 points
- genome-wide dispersion: **2.0**
  - scored from normalized Shannon entropy of best-source-overlap variants across
    fixed native-coordinate mtDNA bins
  - entropy is normalized by log(total bins), so 0 indicates concentration in
    one bin and 1 indicates an even distribution across the full mtDNA
  - scoring:
    - >=0.90: 2.0
    - 0.85-0.90: 1.5
    - 0.75-0.85: 1.0
    - 0.65-0.75: 0.5
    - <0.65: 0
  - occupied bins, circular span, and max-local fraction remain in the report
    for diagnostics only and no longer add score points
- AF coherence (overlap MAD): **1.0**
  - MAD <=0.005: 1.0
  - 0.005-0.01: 0.67
  - 0.01-0.02: 0.33
  - >0.02: 0
  - when overlap <5, this component is capped at 0.67 because very small overlap
    sets can appear artificially coherent
- mirror support: **0.5**
  - use calibrated normalized mirror support only when p95/p99 negative-control
    thresholds are available
  - uncalibrated raw mirror fraction remains reported diagnostically but
    contributes **0 points** to the score

`n_lowA` is retained as an evidence-sufficiency gate rather than an additive
score component. The report records `contamination_score_version =
v7_project_cohort_adjusted` and includes the source and dispersion composite
indices plus project/cohort provenance fields for auditability.

Initial report-only interpretations on the normalized 0-1 scale are:
`>=0.70 strong_evidence`, `0.50-0.699 candidate_evidence`,
`0.30-0.499 weak_ambiguous_evidence`, and `<0.30 little_evidence`.
These bins are exploratory and are not production FAIL/PASS cutoffs.

## Same-species donor specificity diagnostics

The report also quantifies whether best-source-matched low-A alleles are
specific to the nominated donor or common across same-species source candidates.
These diagnostics remain visible for audit, and the same-species background
frequency is now used to background-adjust the source-matching component of the
report-only 0-1 score.

For every allele in the best-source overlap, the analysis counts how many
same-species source candidates carry that allele at source-high AF. It then
reports:

- `best_overlap_n_unique_to_best_source`: number of overlap alleles present at
  source-high AF in only one candidate source sample.
- `best_overlap_unique_to_best_source_fraction`: unique-to-best-source alleles
  divided by total best overlap.
- `best_overlap_mean_source_high_fraction` and
  `best_overlap_median_source_high_fraction`: how common the matched alleles are
  across same-species source candidates.
- `best_overlap_mean_donor_specificity`: mean normalized specificity, where an
  allele carried by only one source has specificity 1 and an allele carried by
  all source candidates has specificity 0.
- `best_overlap_effective_specific_overlap`: sum of per-marker specificity
  across the best overlap; this acts like an overlap count discounted for
  same-species background frequency.
- `best_overlap_fraction_common_ge50`: fraction of best-overlap alleles carried
  at source-high AF by at least half of same-species source candidates.
- `donor_specificity_assessable`: false when there is only one source candidate,
  because donor specificity cannot be estimated relative to background.

The main report also records:
- `target_project`, `target_cohort`
- `best_source_project`, `best_source_cohort`
- `target_source_same_project`, `target_source_same_cohort`
- `provenance_relationship_basis`, `provenance_source_factor`
- `contamination_score_source_total_pre_provenance` and the provenance-adjusted
  `contamination_score_source_total`

Per-variant details are written to
`reports/donor_specificity_variant_detail.tsv`, including the matched allele,
target low VAF, number/fraction of high-A source carriers, normalized donor
specificity, the background weight used for scoring, and the source samples
carrying that allele.

The main report additionally includes
`best_overlap_background_adjusted`,
`best_frac_lowA_in_highB_background_adjusted`, and
`best_overlap_background_adjustment_basis`. The score records
`contamination_score_source_overlap_input`,
`contamination_score_source_fraction_input`, and
`contamination_score_source_basis`.

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
