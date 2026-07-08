# Collect variant calling results

This preprocessing QC step standardizes variant calling outputs before coordinate
liftover. It does **not** perform coordinate liftover, modify VCF records, or edit
mtCN files. It collects each sample's mtCN and VCF outputs by symlink (default) or
copy and creates a merged max-depth per-base coverage file.

## Expected input layout

`--input-root` must contain one directory per sample. The script searches
recursively inside each sample directory, so files can be under subdirectories such
as `alignment/`, `round_1/`, `round_1_variant_calling_decoy/`, `round_2/`, or
`round_2_variant_calling_original_coords/`.

Required file patterns for sample `SAMPLE` are:

* `SAMPLE.round2.mtcn.tsv`
* `SAMPLE.round2.original_coords.clean.final.split.vcf.gz`
  * If the gzipped VCF is absent, `SAMPLE.round2.original_coords.clean.final.split.vcf` is accepted.
* `SAMPLE.numt_decoy.clean.realigned.per_base_coverage.tsv`
* `SAMPLE.round2.original_coords.per_base_coverage.tsv`

## Outputs

Given `--outdir OUT`, the step creates:

* `OUT/collected_mtcn/SAMPLE.round2.mtcn.tsv`
* `OUT/collected_vcf/SAMPLE.round2.original_coords.clean.final.split.vcf.gz`
  * If only an uncompressed input VCF exists, the collected symlink/copy keeps the `.vcf` extension so the file content is not modified or mislabeled.
* `OUT/collected_cov/SAMPLE.merged.max_depth.per_base_coverage.tsv`
* `OUT/reports/variant_calling_collection_summary.tsv`
* `OUT/logs/collection_warnings.log`

The merged coverage file always has four tab-delimited columns:

```text
chrom	pos	target	coverage
```

Coverage files are joined by `chrom`, `pos`, and `target`; the output coverage is
the maximum of the decoy and round-2 original-coordinate coverage values for each
position. If one coverage file is missing, the sample fails by default. Use
`--allow-single-cov` to merge/write the available coverage file and add a warning
note.

## Summary report

The summary table contains one row for every sample directory found under
`--input-root`, including samples with missing inputs. Columns are:

```text
sample species mt_median_coverage nuclear_median_coverage mtcn_median Percent_100 MAD n_hetero n_homo vcf_file cov_file mtcn_file status missing_files notes
```

Species are read from an optional metadata TSV/CSV configured with `--metadata`,
`--metadata-sample-column`, and `--metadata-species-column`. The same two-column
`config/sample_ref_file.tsv` used by coordinate liftover can be supplied here;
headerless files are interpreted as `sample` and `species`. If no species can be
found, `NA` is written.

mtCN columns can be configured explicitly with `--mtcn-mt-column`,
`--mtcn-nuclear-column`, and `--mtcn-mtcn-column`. If these are not supplied, the
script tries common aliases such as `mt_median_coverage`,
`nuclear_median_coverage`, and `mtcn_median`. Missing mtCN metrics are reported as
`NA` and logged in `collection_warnings.log`.

## QC metrics

* `Percent_100` is calculated as
  `100 * count(coverage > coverage_threshold) / total_positions` from the merged
  coverage file. The default `--coverage-threshold` is `100`.
* `MAD` is calculated as
  `median(abs(coverage - median(coverage))) / median(coverage)`, matching R's
  `mad(coverage, constant = 1) / median(coverage)`. If median coverage is zero,
  `MAD` is `NA`.
* VCF heteroplasmy/homoplasmy counts use allele frequencies parsed from the sample
  `FORMAT` `AF` field. Multiple AF values are counted independently.
  * heteroplasmy: `AF >= --low-hetero` and `AF < --low-homo`
  * homoplasmy: `AF >= --low-homo`
  * Only `PASS` (and `.`) variants are counted by default. Use
    `--include-filtered` to count filtered variants too.

## Usage

```bash
python qc_analysis/scripts/collect_variant_calling_results.py \
  --input-root /path/to/variant_calling_root \
  --outdir /path/to/collected_variant_calling_results \
  --metadata config/sample_ref_file.tsv \
  --low-hetero 0.05 \
  --low-homo 0.95 \
  --coverage-threshold 100
```

Config mode is also supported via the shared preprocessing-QC YAML file:

```bash
python qc_analysis/scripts/collect_variant_calling_results.py \
  --config config/qc_preprocessing.yaml
```

Set `collect_variant_calling.copy_files: true` in `config/qc_preprocessing.yaml`
or pass `--copy-files` to copy VCF and mtCN files instead of creating symlinks.
