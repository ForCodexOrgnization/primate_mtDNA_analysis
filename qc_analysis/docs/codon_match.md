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

Positions must be positive integers; strand is `+`/`-`; phase is 1–3; codons are exactly three valid IUPAC DNA symbols (`ACGTRYSWKMBDHVN`); required genomic reference bases are exactly one symbol from the same alphabet; and keys are nonempty. Codons and reference bases are normalized to uppercase. A/C/G/T bases are **resolved**; R/Y/S/W/K/M/B/D/H/V/N are valid but **ambiguous**. Valid ambiguity is retained and counted during preflight rather than treated as fatal, while empty required values, invalid symbols, invalid lengths, conflicting sample mappings, and positions containing conflicting resolved bases (for example A and G) remain fatal. A resolved and an ambiguous annotation at one position (for example G and N), or multiple ambiguous annotations, are not inconsistent and remain distinct. Preflight reports `source_ambiguous_ref_base_rows`, `source_ambiguous_ref_base_positions`, `source_resolved_ref_base_rows`, and `source_conflicting_resolved_ref_base_positions`; only the last is fatal. It writes `codon_annotation_ambiguous_codons.tsv` and `codon_annotation_ambiguous_reference_bases.tsv` diagnostics when applicable. Validate without a VCF, optionally writing overlap diagnostics, with:

```bash
python qc_analysis/scripts/run_codon_match.py --config config/qc_preprocessing.yaml --validate-inputs
python qc_analysis/scripts/run_codon_match.py --config config/qc_preprocessing.yaml --validate-inputs --report-overlaps overlaps.tsv
bash qc_analysis/scripts/run_qc_preprocessing.sh codon_match_validate config/qc_preprocessing.yaml
```

## Matching and overlaps

Every CDS annotation at a position is retained, including ATP8/ATP6 and ND4L/ND4. Exact biological duplicate rows are removed before matching and reported through `MTCODON_DUPLICATE_ANNOTATIONS` and an optional duplicate diagnostics report. Overlap means **more than one unique nonempty gene**, not more than one raw row. Annotation counts and pair counts use deduplicated rows. Every source-human pair is evaluated and deterministic scoring/tie-breaking selects the representative.

Alternate codons are constructed only for single-base A/C/G/T `SRC_REF` and `SRC_ALT`, a fully resolved A/C/G/T source codon, and a resolved matching `ref_base_genome` for that particular annotation. Minus-strand ALT is complemented. Genomic-orientation `SRC_REF` is compared directly (never complemented) with each annotation's `ref_base_genome`. A disagreement between two resolved bases is a confirmed `SOURCE_REF_MISMATCH`; an ambiguous genomic base makes the comparison unknown (`MTCODON_SOURCE_REF_MATCH=NA`) and the alternate codon `.` rather than inventing a base or reporting a mismatch. A candidate can pass only when its source reference base and both codons are resolved, its source reference agrees with `SRC_REF`, and gene, phase, and codon comparisons all match. Thus a resolved overlapping candidate can still pass when another annotation is ambiguous.

`strict_gene_phase_status` is the preferred setting. It only chooses whether failures are categorized as `GENE_MISMATCH`/`PHASE_MISMATCH` or collapsed to `MISMATCH`; phase-mismatched variants never pass. Legacy `strict_phase_match` remains supported. The preferred value wins, with a warning only when both values conflict.

Status precedence is: `MISSING_COORD`, `SKIPPED_NONCODING`, `NO_HUMAN_CODON`, `SOURCE_REF_MISMATCH`, `UNSUPPORTED_NON_SNV`, `PASS`, `AMBIGUOUS_SOURCE_REF`, `AMBIGUOUS_CODON`, `GENE_MISMATCH`, `PHASE_MISMATCH`, `MISMATCH`. `AMBIGUOUS_SOURCE_REF` means compatible candidates exist but all have unresolved genomic reference bases. `AMBIGUOUS_CODON` means gene/phase-compatible candidates with usable reference bases exist but none of those pairs has two resolved codons. Neither ambiguity status masks a gene or phase mismatch.

## INFO schema

All fields below have `Number=1,Type=String` except the explicit groups:

* String scalar: `MTCODON_STATUS`, `MTCODON_SUPPORTED_SNV`, `MTCODON_MATCH`, `MTCODON_STRICT_PHASE`, `MTCODON_GENE_MATCH`, `MTCODON_PHASE_MATCH`, `MTCODON_PRIMATE_GENE`, `MTCODON_PRIMATE_CODON`, `MTCODON_PRIMATE_ALT_CODON`, `MTCODON_PRIMATE_PHASE`, `MTCODON_HUMAN_GENE`, `MTCODON_HUMAN_CODON`, `MTCODON_HUMAN_PHASE`, `MTCODON_OVERLAPPING_CDS`, `MTCODON_AMBIGUOUS_BEST_MATCH`, `MTCODON_SOURCE_REF_MATCH`, `MTCODON_SOURCE_REF_RESOLVED`, `MTCODON_ANY_RESOLVED_SOURCE_REF`, `MTCODON_DUPLICATE_ANNOTATIONS`, `MTCODON_SOURCE_CODON_RESOLVED`, `MTCODON_HUMAN_CODON_RESOLVED`, `MTCODON_ANY_RESOLVED_PAIR`. Representative resolved fields use `yes`, `no`, or `NA`; any-candidate fields use `yes` or `no`. `MTCODON_SOURCE_REF_RESOLVED=no` identifies a selected valid ambiguous reference base, while `NA` means there is no selected annotation or value.
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
