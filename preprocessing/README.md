# Preprocessing workflow for primate mtDNA QC

This module prepares species-level and sample-level inputs for mitochondrial variant calling and in-house NUMT score analysis.

## Step 0. Reference discovery

`preprocessing/scripts/run_reference_discovery.sh` runs `find_primate_wg_chrM_refs.py` on `data/metadata/preprint_species_list.tsv`, a local RefSeq mitochondrion FASTA, and a primate tree. The raw output is written to `results/preprocessing/reference_discovery/species_reference_chrM_summary.tsv`. After manual review, copy or symlink the stable manifest to `data/metadata/species_reference_chrM_summary.tsv`.

## Step 1. Reference materialization

`preprocessing/scripts/build_reference_materialization_manifest.R` reads the reviewed discovery manifest and writes `reference_materialization_manifest.tsv` to both `results/preprocessing/reference_materialization/` and `references/manifests/`. It classifies chrM as `embedded_in_wg_ref`, `independent_chrM_ref`, or `missing_chrM_ref` and constructs expected local FASTA paths.

`preprocessing/scripts/materialize_references.sh` downloads WG FASTA files, indexes them, extracts embedded chrM records from the matching WG FASTA, downloads or extracts independent chrM FASTA files, and writes `references/manifests/in_house_score_reference_inputs.tsv`.

## Step 2. In-house score and minimal NUMT mask selection

The in-house score scripts consume `references/manifests/in_house_score_reference_inputs.tsv`. Embedded chrM rows must use chrM FASTA extracted from the WG reference; independent chrM rows must use the independent chrM FASTA. The merged result should be `results/preprocessing/in_house_score/merged_in_house_score.tsv`.

## Step 3. Prepare variant-calling inputs

`preprocessing/scripts/prepare_variant_calling_inputs.R` joins `data/metadata/sample_metadata.tsv`, `references/manifests/in_house_score_reference_inputs.tsv`, and the merged in-house score table to produce `results/preprocessing/variant_calling_inputs/variant_calling_input_table.tsv`. Downstream mtDNA variant calling should use this table as its main input, but the downstream pipeline is not modified here.
