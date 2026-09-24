# Interspecies contamination QC

This cohort-level step runs after coordinate liftover and before MITOS2/codon
annotation. It reads
`results/qc/coordinate_liftover/vcf_lifted_raw/{sample}.lifted.raw.vcf[.gz]`.
It never rewrites a lifted VCF.

The historical `PASS/WARN/FAIL` classification remains the production
classification consumed by `final_filter`. A separate 0-1 evidence score is
reported for ranking, calibration, and sensitivity review; the score does not
currently change production classification.

## Production classification

Only biallelic canonical SNVs with `FILTER=PASS`, `DP >= 100`, and a usable
AF are indexed. Recipient alleles use `0.01 <= AF <= 0.20`; potential source
alleles use `AF >= 0.99`. Allele identity is exact post-liftover
`CHROM/POS/REF/ALT`.

Before cross-species matching, a recipient low-VAF allele is removed if another
sample of the recipient species carries the same allele at high AF. This
recipient-species background removal reduces false signals caused by
species-level polymorphism/reference divergence.

The remaining low-VAF alleles are matched to high-AF alleles in other species.
The report identifies the best source species and the best individual source
sample. Historical production thresholds are preserved:

- at least 5 informative low-VAF alleles after recipient-species background
  removal;
- best source-species overlap >= 3 and overlap fraction >= 0.50;
- best source-sample overlap >= 3 and overlap fraction >= 0.50;
- at least 70% of matched recipient AFs within median AF +/- 0.03;
- source species must be unambiguous;
- recipient species must have more than one cohort sample for recipient
  background to be estimable.

A coherent, unambiguous signal meeting these criteria is
`INTERSPECIES_CONTAMINATION / FAIL`. Insufficient, ambiguous, incoherent, or
source-sample-unsupported signals remain `WARN`; signals below the configured
cross-species thresholds remain `PASS`.

## Report-only evidence score

The 0-1 score is `v2_cross_species_separation_project_cohort`. It is
calculated only when the minimum evidence gate passes
(`n_lowA_after_species_background >= 5`, best-source overlap >= 3, and a best
source species exists).

The score contains five components totaling 10 points before normalization:

1. **Source matching: up to 4.0 points**
   - Each matched allele is down-weighted according to how common that allele is
     among eligible non-recipient species.
   - Species, rather than sample, is the frequency unit so densely sampled taxa
     do not dominate the background.
   - Default source-specificity weights are:
     - <=5% of eligible species: 1.00
     - >5-10%: 0.75
     - >10-25%: 0.50
     - >25-50%: 0.25
     - >=50%: 0
   - Adjusted overlap and adjusted overlap fraction are converted to strengths
     and combined by geometric mean.
   - If fewer than two other species are available, specificity is reported as
     unassessable and raw overlap is used rather than inventing specificity.

2. **Genome-wide dispersion: up to 2.0 points**
   - Uses normalized Shannon entropy of best-source matched positions across
     1-kb human chrM bins.
   - Higher entropy supports a genome-wide donor signal rather than a local
     alignment/indel artifact.
   - Circular span and maximum local fraction are also reported diagnostically.

3. **AF coherence: up to 1.5 points**
   - Based on MAD of matched recipient low-VAF alleles.
   - MAD <=0.005: 1.5
   - <=0.01: 1.0
   - <=0.02: 0.5
   - >0.02: 0
   - With fewer than five overlaps this component is capped at 1.0.

4. **Best-source-sample concentration: up to 1.5 points**
   - Defined from specificity-adjusted overlap in the best individual source
     sample divided by specificity-adjusted overlap in the best source species.
   - This distinguishes a plausible individual donor from a signal broadly
     distributed across the source species.

5. **Best-vs-runner-up source-species separation: up to 1.0 point**
   - Rank source species by specificity-adjusted overlap.
   - Define separation as `(best - second_best) / best`.
   - This avoids denominator inflation from counting the same shared allele in
     many candidate species.
   - Separation >=0.60: 1.0 point; >=0.40: 0.67; >=0.20: 0.33; <0.20: 0.

The normalized score is `raw_10 / 10`. Initial report-only interpretation bins
mirror the intraspecies framework:

- >=0.70: `strong_evidence`
- 0.50-0.699: `candidate_evidence`
- 0.30-0.499: `weak_ambiguous_evidence`
- <0.30: `little_evidence`

These bins are exploratory and are not production FAIL/PASS thresholds.

## Project/cohort provenance

When project/cohort metadata are available, target-source provenance modifies
only the **source-matching** component. It does not independently create
contamination evidence.

Default factors are:

- same cohort: x1.00
- same project, cohort unknown: x0.90
- different cohort without a known project mismatch: x0.75
- different project: x0.60
- missing metadata: x1.00 (neutral)

The preferred provenance source for the deduplicated workflow is
`results/qc/sample_deduplication/reports/deduplicated_sample_ref_file.tsv`.
If project/cohort columns are already present in `sample_ref_file`, those
values are also used. Missing provenance is deliberately neutral.

## Key output fields

In addition to the historical columns, the report now contains:

- `target_project`, `target_cohort`, `best_source_project`,
  `best_source_cohort`;
- `target_source_same_project`, `target_source_same_cohort`,
  `provenance_relationship_basis`, `provenance_source_factor`;
- specificity-adjusted source overlap/fraction and cross-species allele
  prevalence;
- best-source-sample concentration, best/second adjusted overlaps, and
  best-vs-runner-up source-species separation;
- dispersion and AF-MAD diagnostics;
- score component fields, raw 10-point score, normalized score, and
  interpretation.

## Run

```bash
python qc_analysis/scripts/run_interspecies_contamination.py \
  --config config/qc_preprocessing.yaml

bash qc_analysis/scripts/run_qc_preprocessing.sh \
  --submit interspecies_contamination config/qc_preprocessing.yaml
```

Output:

`results/qc/interspecies_contamination/reports/interspecies_contamination_report.tsv`
