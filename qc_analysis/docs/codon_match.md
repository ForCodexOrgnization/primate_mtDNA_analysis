# Codon match

`run_codon_match.py` validates its inputs, builds compact one-to-many indexes, resolves `sample -> reference_key`, annotates every lifted VCF record, and atomically publishes a VCF and per-sample summary. It never filters records. A separate merge operation atomically creates the all-samples summary.

## Inputs and preflight

Strict validation is the default (`codon_match.settings.strict_input_validation: true`). Required columns are:

| table | required columns |
|---|---|
| reference codon | `reference_key`, `pos`, `gene`, `strand`, `codon_pos_in_triplet`, `codon_seq`, `ref_base_genome` |
| historical sample codon | `sample`, `pos`, `gene`, `strand`, `codon_pos_in_triplet`, `codon_seq`, `ref_base_genome` |
| human codon | `pos`, `gene`, `strand`, `codon_pos_in_triplet`, `codon_seq` |
| sample-reference map | `sample`, `reference_key` |

Positions must be positive integers; strand is `+`/`-`; phase is 1–3; codons are exactly three valid IUPAC DNA symbols (`ACGTRYSWKMBDHVN`); genomic reference bases are A/C/G/T; and keys are nonempty. Codons are normalized to uppercase. Valid ambiguous codons are retained and counted during preflight rather than treated as fatal, while empty codons, invalid symbols, invalid lengths, conflicting sample mappings, and inconsistent genomic bases remain fatal. Preflight writes `codon_annotation_ambiguous_codons.tsv` when ambiguity is present. Validate without a VCF, optionally writing overlap diagnostics, with:

```bash
python qc_analysis/scripts/run_codon_match.py --config config/qc_preprocessing.yaml --validate-inputs
python qc_analysis/scripts/run_codon_match.py --config config/qc_preprocessing.yaml --validate-inputs --report-overlaps overlaps.tsv
bash qc_analysis/scripts/run_qc_preprocessing.sh codon_match_validate config/qc_preprocessing.yaml
```

## Matching and overlaps

Every CDS annotation at a position is retained, including ATP8/ATP6 and ND4L/ND4. Exact biological duplicate rows are removed before matching and reported through `MTCODON_DUPLICATE_ANNOTATIONS` and an optional duplicate diagnostics report. Overlap means **more than one unique nonempty gene**, not more than one raw row. Annotation counts and pair counts use deduplicated rows. Every source-human pair is evaluated and deterministic scoring/tie-breaking selects the representative.

Alternate codons are constructed only for single-base A/C/G/T `SRC_REF` and `SRC_ALT` and a fully resolved A/C/G/T source codon. Minus-strand ALT is complemented. Before construction, genomic-orientation `SRC_REF` is compared directly (never complemented) with `ref_base_genome`; a disagreement produces `SOURCE_REF_MISMATCH`, `MTCODON_MATCH=no`, and alternate codon `.`. An exact match is possible only when both source and human codons are fully resolved; ambiguity is not expanded into possible codons.

`strict_gene_phase_status` is the preferred setting. It only chooses whether failures are categorized as `GENE_MISMATCH`/`PHASE_MISMATCH` or collapsed to `MISMATCH`; phase-mismatched variants never pass. Legacy `strict_phase_match` remains supported. The preferred value wins, with a warning only when both values conflict.

Status precedence is: `MISSING_COORD`, `SKIPPED_NONCODING`, `NO_HUMAN_CODON`, `SOURCE_REF_MISMATCH`, `UNSUPPORTED_NON_SNV`, `PASS`, `AMBIGUOUS_CODON`, `GENE_MISMATCH`, `PHASE_MISMATCH`, `MISMATCH`. `AMBIGUOUS_CODON` means gene/phase-compatible annotations exist but none of those pairs has two resolved codons; it does not mask a gene or phase mismatch.

## INFO schema

All fields below have `Number=1,Type=String` except the explicit groups:

* String scalar: `MTCODON_STATUS`, `MTCODON_SUPPORTED_SNV`, `MTCODON_MATCH`, `MTCODON_STRICT_PHASE`, `MTCODON_GENE_MATCH`, `MTCODON_PHASE_MATCH`, `MTCODON_PRIMATE_GENE`, `MTCODON_PRIMATE_CODON`, `MTCODON_PRIMATE_ALT_CODON`, `MTCODON_PRIMATE_PHASE`, `MTCODON_HUMAN_GENE`, `MTCODON_HUMAN_CODON`, `MTCODON_HUMAN_PHASE`, `MTCODON_OVERLAPPING_CDS`, `MTCODON_AMBIGUOUS_BEST_MATCH`, `MTCODON_SOURCE_REF_MATCH`, `MTCODON_DUPLICATE_ANNOTATIONS`, `MTCODON_SOURCE_CODON_RESOLVED`, `MTCODON_HUMAN_CODON_RESOLVED`, `MTCODON_ANY_RESOLVED_PAIR`. The representative resolved fields use `yes`, `no`, or `NA`; the any-pair field uses `yes` or `no`.
* `Number=1,Type=Integer`: `MTCODON_N_PRIMATE_ANNOTATIONS`, `MTCODON_N_HUMAN_ANNOTATIONS`, `MTCODON_N_PAIR_CANDIDATES`.
* `Number=.,Type=String`: `MTCODON_PRIMATE_GENES`, `MTCODON_HUMAN_GENES`, `MTCODON_MATCHING_GENES`.

VCFs also receive version and source-table provenance metadata. Summaries record version, table paths, reference key, and effective strict settings.

## Sequential and parallel operation

Each annotation invocation writes only `<sample>.codon_match_summary.tsv`; it never touches the merged report, so independent sample jobs are safe:

```bash
python qc_analysis/scripts/run_codon_match.py --config config/qc_preprocessing.yaml --sample SAMPLE_A
python qc_analysis/scripts/run_codon_match.py --config config/qc_preprocessing.yaml --sample SAMPLE_B
python qc_analysis/scripts/run_codon_match.py --config config/qc_preprocessing.yaml --merge-summaries
# wrapper equivalent
bash qc_analysis/scripts/run_qc_preprocessing.sh codon_match_merge config/qc_preprocessing.yaml
```

Merge scans only per-sample summaries, checks identical schemas, rejects conflicting rows for a sample, sorts samples, and atomically publishes `all_samples.codon_match_summary.tsv`. The wrapper `all` workflow runs one multi-sample annotation owner and then merges; parallel sample tasks must not run the merge themselves. VCFs, summaries, merge output, and diagnostics use same-directory temporary files followed by `os.replace()`.
