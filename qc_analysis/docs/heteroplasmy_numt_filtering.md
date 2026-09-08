# Local heteroplasmy / NUMT analysis

This module runs in the **native primate mitochondrial coordinate system before liftover**.  It is intentionally report-only until the terminal filtering step.

## Production order

The intended biological order is:

1. collect variant-calling outputs
2. sample-level QC report
3. pre-liftover immutable source-call QC
4. **local heteroplasmy analysis**
5. intra-species contamination report
6. reference anchor discovery and coordinate liftover
7. human contamination report
8. inter-species contamination report
9. codon / tRNA / rRNA matching
10. terminal sample/variant filtering
11. downstream analysis

`human_contamination` is currently a standalone step in the legacy `all` wrapper; this feature branch documents it explicitly between liftover and inter-species contamination.

## Heteroplasmy definition

Cluster discovery uses native-coordinate variants satisfying:

- source call class `HET`
- VCF `FILTER=PASS`
- source variant QC `PASS`
- source DP > 100
- 0.10 < source AF < 0.95

## Local cluster detection

Default parameters are defined in `config/heteroplasmy_numt.yaml`:

- circular 250-bp window
- whole-cluster AF span <= 0.06
- minimum seed size 3
- 2,000 permutations
- empirical P <= 0.01

The permutation null randomizes positions uniformly without replacement while retaining the observed AF values.  A sample-specific critical cluster count is derived from the empirical null.  Independent seeds are selected greedily, then expanded post hoc by adding unassigned variants within 250 bp of a member while preserving the whole-cluster AF span <= 0.06.  Expansion does not alter the seed P value.

## NUMT annotation

Cluster discovery and NUMT annotation are separate stages.

### Sample-level evidence

Each cluster is compared with the target sample's upstream NUMT besthit intervals.  Upstream `highconf_numt.bed` overlap defines `HIGH_CONF_NUMT`; other upstream besthits are `BESTHIT_ONLY_NUMT`.  No new downstream threshold is imposed on pident, alignment length, read count, or MAPQ.

A regional NUMT overlap requires, by default, at least 2 cluster variants and at least 50% of cluster variants to lie in the NUMT -> chrM interval.

### Species-level evidence

All sample NUMT intervals are pooled by `species + reference_key` and merged into a species catalogue.  A target sample can receive species-level NUMT support only from **other samples** using the same species/reference coordinate system.  The output retains the number and identities of supporting samples.

This species-level layer is intended to rescue sample-specific NUMT-discovery false negatives; it does not replace direct sample-level evidence.

## Same-species recurrence

For every cluster, the same native-coordinate `POS + REF + ALT` alleles are compared with other samples from the same species.  The default recurrent-pattern rule is:

- >=2 shared variants
- >=50% of cluster variants shared
- median absolute AF difference <=0.05

NUMT evidence and recurrence are computed independently.

## Four cluster classes

| NUMT evidence | recurrence | class |
|---|---|---|
| yes | yes | `NUMT_RECURRENT` |
| yes | no | `NUMT_ONLY` |
| no | yes | `RECURRENT_ONLY` |
| no | no | `UNRESOLVED` |

`NUMT evidence=yes` includes either direct sample-level evidence or species-level evidence.  The detailed source/tier remains in the cluster report.

## Conservative removal policy

Default actions are:

- sample-level NUMT overlap -> `REMOVE`
- species-level NUMT + recurrence -> `REMOVE`
- species-level NUMT only -> `FLAG`
- recurrence only -> `FLAG`
- unresolved -> `KEEP`

NUMT-associated samples are **not removed wholesale** by default.  The module marks `numt_sample=YES` and emits the specific variants to remove.  Sample exclusion remains the responsibility of sample QC and contamination reports.

## Outputs

`results/qc/local_heteroplasmy_qc/reports/` contains:

- `local_heteroplasmy_sample_summary.tsv`
- `local_heteroplasmy_cluster_summary.tsv`
- `local_heteroplasmy_variant_detail.tsv`
- `species_numt_intervals.tsv`
- `numt_samples.tsv`
- `numt_variants_to_remove.tsv`
- `heteroplasmy_analysis_summary.tsv`

The removal table uses the immutable source key:

`sample + source_chrom + source_pos + source_ref + source_alt`

## Terminal filtering

Run the normal terminal final filter first, then apply the heteroplasmy blacklist using the immutable `SOURCE_*` INFO annotations that survive coordinate liftover:

```bash
python qc_analysis/scripts/run_final_filter_with_heteroplasmy.py \
  --config config/qc_preprocessing.yaml
```

The wrapper:

1. runs `run_final_filter_streaming.py`;
2. reads `numt_variants_to_remove.tsv`;
3. marks matching rows in `reports/final_variant_qc.tsv` as `FAIL`;
4. annotates `reports/final_sample_qc.tsv` with NUMT sample/variant counts;
5. removes matching records from `final_vcf/` by `SOURCE_*` identity;
6. writes `reports/heteroplasmy_final_filter_summary.tsv`.

## Running the analysis while this feature is under review

```bash
# prerequisite native-coordinate source QC
bash qc_analysis/scripts/run_qc_preprocessing.sh \
  pre_liftover_variant_qc config/qc_preprocessing.yaml

# heteroplasmy/NUMT analysis
python qc_analysis/scripts/run_local_heteroplasmy_qc.py \
  --config config/heteroplasmy_numt.yaml

# continue biological QC in the intended order
bash qc_analysis/scripts/run_qc_preprocessing.sh \
  intraspecies_contamination config/qc_preprocessing.yaml
bash qc_analysis/scripts/run_qc_preprocessing.sh \
  discover_global_anchor config/qc_preprocessing.yaml
bash qc_analysis/scripts/run_qc_preprocessing.sh \
  coordinate_liftover config/qc_preprocessing.yaml
bash qc_analysis/scripts/run_qc_preprocessing.sh \
  human_contamination config/qc_preprocessing.yaml
bash qc_analysis/scripts/run_qc_preprocessing.sh \
  interspecies_contamination config/qc_preprocessing.yaml

# after codon/tRNA/rRNA annotation, run the terminal filter with NUMT removal
python qc_analysis/scripts/run_final_filter_with_heteroplasmy.py \
  --config config/qc_preprocessing.yaml
```

Before production use, verify the two external NUMT directories in `config/heteroplasmy_numt.yaml` against the current HPC NUMT-discovery outputs.
