# Preprocessing workflow for primate mtDNA QC

This module prepares species-level and sample-level inputs for mitochondrial variant calling and in-house NUMT score analysis.

## Quick start

1. Edit `config/preprocessing_paths.yaml` so the species list, local RefSeq mitochondrion FASTA, primate tree, sample metadata, and output paths match your HPC workspace.
2. Run reference discovery:

   ```bash
   bash preprocessing/scripts/run_preprocessing.sh reference_discovery config/preprocessing_paths.yaml
   ```

3. Manually review `results/preprocessing/reference_discovery/species_reference_chrM_summary.tsv`. Confirm that the selected WG and chrM references are biologically appropriate before any downstream step.
4. Copy or symlink the reviewed manifest to `data/metadata/species_reference_chrM_summary.tsv`:

   ```bash
   cp results/preprocessing/reference_discovery/species_reference_chrM_summary.tsv data/metadata/species_reference_chrM_summary.tsv
   # or:
   ln -sf ../../results/preprocessing/reference_discovery/species_reference_chrM_summary.tsv data/metadata/species_reference_chrM_summary.tsv
   ```

5. Run reference materialization:

   ```bash
   bash preprocessing/scripts/run_preprocessing.sh reference_materialization config/preprocessing_paths.yaml
   ```

6. Run in-house score and minimal NUMT mask selection:

   ```bash
   bash preprocessing/scripts/run_preprocessing.sh in_house_score config/preprocessing_paths.yaml
   ```

7. Prepare variant-calling inputs:

   ```bash
   bash preprocessing/scripts/run_preprocessing.sh variant_inputs config/preprocessing_paths.yaml
   ```

8. Optionally render Quarto reports:

   ```bash
   bash preprocessing/scripts/run_preprocessing.sh reports config/preprocessing_paths.yaml
   ```

After manual reference review, you may run the downstream preprocessing stages together:

```bash
bash preprocessing/scripts/run_preprocessing.sh post_reference_review config/preprocessing_paths.yaml
```

The `all` command intentionally runs only reference discovery and then stops, because `species_reference_chrM_summary.tsv` requires manual review before reference materialization, in-house scoring, or variant input preparation.

Downloaded WG FASTA files, chrM FASTA files, FASTA indexes, BLAST outputs, and other generated reference artifacts are HPC-local outputs. Do not commit these large downloaded/generated reference files to GitHub; commit only small metadata/configuration files and reviewed manifests when appropriate.

## Step 0. Reference discovery

`preprocessing/scripts/run_reference_discovery.sh` runs `find_primate_wg_chrM_refs.py` on `data/metadata/all_species_list.txt`, a local RefSeq mitochondrion FASTA, and a primate tree. The raw output is written to `results/preprocessing/reference_discovery/species_reference_chrM_summary.tsv`. After manual review, copy or symlink the stable manifest to `data/metadata/species_reference_chrM_summary.tsv`.

## Step 1. Reference materialization

`preprocessing/scripts/build_reference_materialization_manifest.R` reads the reviewed discovery manifest and writes `reference_materialization_manifest.tsv` to both `results/preprocessing/reference_materialization/` and `references/manifests/`. It classifies chrM as `embedded_in_wg_ref`, `independent_chrM_ref`, or `missing_chrM_ref` and constructs expected local FASTA paths.

`preprocessing/scripts/materialize_references.sh` downloads WG FASTA files, indexes them, extracts embedded chrM records from the matching WG FASTA, downloads or extracts independent chrM FASTA files, and writes `references/manifests/in_house_score_reference_inputs.tsv`.

## Step 2. In-house score and minimal NUMT mask selection

The in-house score scripts consume `references/manifests/in_house_score_reference_inputs.tsv`. Embedded chrM rows must use chrM FASTA extracted from the WG reference; independent chrM rows must use the independent chrM FASTA. The merged result should be `results/preprocessing/in_house_score/merged_in_house_score.tsv`.

## Step 3. Prepare variant-calling inputs

`preprocessing/scripts/prepare_variant_calling_inputs.R` joins `data/metadata/sample_metadata.tsv`, `references/manifests/in_house_score_reference_inputs.tsv`, and the merged in-house score table to produce `results/preprocessing/variant_calling_inputs/variant_calling_input_table.tsv`. Downstream mtDNA variant calling should use this table as its main input, but the downstream pipeline is not modified here.
