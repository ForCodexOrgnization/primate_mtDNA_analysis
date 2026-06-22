# Preprocessing workflow for primate mtDNA QC

This module prepares species-level and sample-level inputs for mitochondrial variant calling and in-house NUMT score analysis.

## Quick start

1. Edit `config/preprocessing_paths.yaml` so the species list, local RefSeq mitochondrion FASTA, primate tree, sample metadata, and output paths match your HPC workspace.

2. If you want one command to run every preprocessing stage, use `all_steps`:

   ```bash
   bash preprocessing/scripts/run_preprocessing.sh all_steps config/preprocessing_paths.yaml
   ```

   `all_steps` copies the raw reference-discovery summary into the reviewed-manifest path and continues through reference materialization, in-house score, and variant input preparation. Use this shortcut only when unreviewed reference choices are acceptable, such as exploratory runs.

If you want to manually review reference choices before downstream preprocessing, use the safer staged workflow:

1. Run reference discovery:

   ```bash
   bash preprocessing/scripts/run_preprocessing.sh reference_discovery config/preprocessing_paths.yaml
   ```

2. Manually review `results/preprocessing/reference_discovery/species_reference_chrM_summary.tsv`. Confirm that the selected WG and chrM references are biologically appropriate before any downstream step.
3. Copy or symlink the reviewed manifest to `data/metadata/species_reference_chrM_summary.tsv`:

   ```bash
   cp results/preprocessing/reference_discovery/species_reference_chrM_summary.tsv data/metadata/species_reference_chrM_summary.tsv
   # or:
   ln -sf ../../results/preprocessing/reference_discovery/species_reference_chrM_summary.tsv data/metadata/species_reference_chrM_summary.tsv
   ```

4. Run reference materialization:

   ```bash
   bash preprocessing/scripts/run_preprocessing.sh reference_materialization config/preprocessing_paths.yaml
   ```

5. Run in-house score and minimal NUMT mask selection:

   ```bash
   bash preprocessing/scripts/run_preprocessing.sh in_house_score config/preprocessing_paths.yaml
   ```

6. Prepare variant-calling inputs:

   ```bash
   bash preprocessing/scripts/run_preprocessing.sh variant_inputs config/preprocessing_paths.yaml
   ```

7. Optionally render Quarto reports:

   ```bash
   bash preprocessing/scripts/run_preprocessing.sh reports config/preprocessing_paths.yaml
   ```

After manual reference review, you may run the downstream preprocessing stages together:

```bash
bash preprocessing/scripts/run_preprocessing.sh post_reference_review config/preprocessing_paths.yaml
```

The `all` command intentionally runs only reference discovery and then stops, because `species_reference_chrM_summary.tsv` requires manual review before reference materialization, in-house scoring, or variant input preparation.

Downloaded WG FASTA files, chrM FASTA files, FASTA indexes, BLAST outputs, and other generated reference artifacts are HPC-local outputs. Do not commit these large downloaded/generated reference files to GitHub; commit only small metadata/configuration files and reviewed manifests when appropriate.

## Input files and sources

Before running preprocessing, make sure these user-provided inputs are available and that `config/preprocessing_paths.yaml` points to them:

- `data/metadata/all_species_list.txt`: the project species list used by reference discovery. This repository now stores the former `preprint_species_list.tsv` under this more general name. It should include a species column accepted by `find_primate_wg_chrM_refs.py` such as `species`, `target_species`, `GENUS_SPECIES`, or `FINAL_PRIMATE_NAME`; optional columns such as `sample_count` and `preprint_REFERENCE_SPECIES` can guide summaries and fallback reference choices.
- `mitochondrion.1.1.genomic.fna.gz`: the local RefSeq mitochondrion release FASTA used to find same-species or fallback chrM references. Download it from the NCBI RefSeq mitochondrion release directory: <https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/>. Store this large FASTA on HPC storage and set `mito_fasta` in `config/preprocessing_paths.yaml` to the downloaded file path.
- Primate Newick tree: the phylogeny used for nearest-species fallback searches. Store it on HPC/project storage and set `tree_newick` in `config/preprocessing_paths.yaml`.
- `data/metadata/sample_metadata.tsv`: the sample table used when preparing final variant-calling inputs. It should include `sample_id`, `target_species` or `species`, `cram_path`, and `cram_index_path`.
- `data/metadata/species_reference_chrM_summary.tsv`: the reviewed reference-discovery manifest. Create it by manually reviewing `results/preprocessing/reference_discovery/species_reference_chrM_summary.tsv`, or let `all_steps` copy the raw discovery output into this path for exploratory unreviewed runs.

## Step 0. Reference discovery

`preprocessing/scripts/run_reference_discovery.sh` runs `find_primate_wg_chrM_refs.py` on `data/metadata/all_species_list.txt`, a local RefSeq mitochondrion FASTA, and a primate tree. The raw output is written to `results/preprocessing/reference_discovery/species_reference_chrM_summary.tsv`. After manual review, copy or symlink the stable manifest to `data/metadata/species_reference_chrM_summary.tsv`.

## Step 1. Reference materialization

`preprocessing/scripts/build_reference_materialization_manifest.R` reads the reviewed discovery manifest and writes `reference_materialization_manifest.tsv` to both `results/preprocessing/reference_materialization/` and `references/manifests/`. It classifies chrM as `embedded_in_wg_ref`, `independent_chrM_ref`, or `missing_chrM_ref` and constructs expected local FASTA paths.

`preprocessing/scripts/materialize_references.sh` downloads WG FASTA files, indexes them, extracts embedded chrM records from the matching WG FASTA, downloads or extracts independent chrM FASTA files, and writes `references/manifests/in_house_score_reference_inputs.tsv`.

## Step 2. In-house score and minimal NUMT mask selection

The in-house score scripts consume `references/manifests/in_house_score_reference_inputs.tsv`. Embedded chrM rows must use chrM FASTA extracted from the WG reference; independent chrM rows must use the independent chrM FASTA. The merged result should be `results/preprocessing/in_house_score/merged_in_house_score.tsv`.

## Step 3. Prepare variant-calling inputs

`preprocessing/scripts/prepare_variant_calling_inputs.R` joins `data/metadata/sample_metadata.tsv`, `references/manifests/in_house_score_reference_inputs.tsv`, and the merged in-house score table to produce `results/preprocessing/variant_calling_inputs/variant_calling_input_table.tsv`. Downstream mtDNA variant calling should use this table as its main input, but the downstream pipeline is not modified here.
