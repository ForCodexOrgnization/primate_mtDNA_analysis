# Global liftover anchor discovery

Mitochondrial genomes are circular, but FASTA files are linear strings. Two homologous mtDNA references can therefore contain the same biology while starting at different FASTA coordinates. The global-anchor workflow discovers a reproducible circular cut point before coordinate liftover.

Core principle:

> Global MSA is used to identify a common homologous circular cut point. Pairwise species–human alignment is used to generate the final coordinate map.

## Why replace the first exact shared k-mer anchor?

The old per-sample helper found the first exact k-mer shared with human and used that immediately as the production rotation anchor. That is brittle because it depends on the pairwise search order, exact identity, FASTA starts, and per-sample repetition. The new workflow uses exact shared k-mers only as a coarse pre-rotation aid so the circular genomes can be linearized similarly for a multi-species alignment.

## Workflow

Run before coordinate liftover:

```bash
python3 qc_analysis/scripts/discover_global_liftover_anchor.py --config config/qc_preprocessing.yaml
```

The script:

1. reads `global_anchor_discovery.sample_ref_file`;
2. resolves each `species_fasta` and selected mitochondrial record;
3. hashes normalized sequence as SHA256;
4. deduplicates identical reference sequences;
5. filters references by length, N fraction, empty sequence, invalid bases, and ID collisions;
6. coarsely rotates eligible references using explicit/sample anchors when provided or pairwise shared k-mers;
7. builds `unique_references.coarse_rotated.fa` containing human plus eligible unique references;
8. runs MAFFT for the multi-species alignment;
9. scores every candidate window;
10. selects one deterministic homologous anchor column;
11. projects that column back to each original FASTA coordinate.

## Reference IDs and hashes

The primary sequence identity is:

```python
sha256(sequence.upper().encode()).hexdigest()
```

Whitespace is removed before hashing. If a manifest row supplies `reference_id`, that ID is used. Otherwise, the stable ID is derived from species, FASTA basename, and the first 12 characters of the sequence SHA256. This same helper is used by discovery and runtime liftover.

## Candidate windows and ranking

Windows are controlled by `candidate_window_size`. Each window reports occupancy, minimum occupancy, major-allele fraction, Shannon entropy, gap fraction, homopolymer run, human gaps, eligible-reference gaps, and distance from the MSA edges. A window is eligible when it satisfies configured thresholds such as `min_window_mean_occupancy`, `max_window_gap_fraction`, and `max_homopolymer_run`.

Eligible windows and columns are ranked deterministically by occupancy, conservation, gap fraction, entropy, edge distance, and finally lower alignment column number. Input filesystem order and manifest order are not used as hidden tie breakers.

## Fallbacks

`GLOBAL_MSA_ANCHOR` rows are preferred. `FAMILY_MSA_ANCHOR` rows can be used when present. `PAIRWISE_FALLBACK` is runtime-only and must be explicitly enabled. Position 1 is never silently assigned in production when validated anchors are required.

## Outputs

`results/qc/global_anchor/` contains:

- `unique_reference_manifest.tsv`
- `excluded_references.tsv`
- `unique_references.coarse_rotated.fa`
- `all_references.aligned.fa`
- `global_anchor_candidates.tsv`
- `global_anchor_selection.tsv`
- `reference_anchor_positions.tsv`
- `reference_anchor_failures.tsv`
- `family_anchor_candidates.tsv`
- `family_anchor_positions.tsv`
- `global_anchor_summary.tsv`

Coordinate liftover consumes `reference_anchor_positions.tsv`.

## Reviewing and overriding

Review `global_anchor_candidates.tsv` for nearby windows with similar conservation or excessive gaps. Manual sample overrides can still use `rotate_anchor` in the sample manifest; these are reported as `SAMPLE_OVERRIDE` in coordinate-liftover QC and take priority over the reference-level table.
