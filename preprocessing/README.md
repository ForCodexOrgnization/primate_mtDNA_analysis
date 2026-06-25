# Preprocessing workflow for primate mtDNA QC

This module prepares species-level reference inputs for mitochondrial variant calling and in-house NUMT score analysis.

## Quick start

1. Edit `config/preprocessing_paths.yaml` so the species list, local RefSeq mitochondrion FASTA, primate tree, tool commands, and output paths match your HPC workspace.

2. If you want one command to run every preprocessing stage, use `all_steps`:

   ```bash
   bash preprocessing/scripts/run_preprocessing.sh all_steps config/preprocessing_paths.yaml
   ```

   `all_steps` copies the raw reference-discovery summary into the reviewed-manifest path and continues through reference materialization, in-house score, and variant-calling reference package generation. Use this shortcut only when unreviewed reference choices are acceptable, such as exploratory runs.

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

6. Build variant-calling reference packages:

   ```bash
   bash preprocessing/scripts/run_preprocessing.sh variant_references config/preprocessing_paths.yaml
   ```

7. Optionally render Quarto reports:

   ```bash
   bash preprocessing/scripts/run_preprocessing.sh reports config/preprocessing_paths.yaml
   ```

After manual reference review, you may run the downstream preprocessing stages together:

```bash
bash preprocessing/scripts/run_preprocessing.sh post_reference_review config/preprocessing_paths.yaml
```

The `all` command intentionally runs only reference discovery and then stops, because `species_reference_chrM_summary.tsv` requires manual review before reference materialization, in-house scoring, or variant-reference package generation.

Downloaded WG FASTA files, chrM FASTA files, FASTA indexes, BLAST outputs, and other generated reference artifacts are HPC-local outputs. Do not commit these large downloaded/generated reference files to GitHub; commit only small metadata/configuration files and reviewed manifests when appropriate.


## HPC environment setup

Before each preprocessing step runs, `preprocessing/scripts/run_preprocessing.sh` can source an optional setup script configured in `config/preprocessing_paths.yaml`. Use this for HPC-specific commands such as `module load`, `conda activate`, or site-specific executable paths. The runner checks required commands before starting each step and fails early with a clear message if a required tool is still unavailable.

Example config:

```yaml
environment_setup_script: "config/preprocessing_hpc_env.sh"
rscript_command: "Rscript"
python_command: "python3"
wget_command: "wget"
samtools_command: "samtools"
bwa_command: "bwa"
gatk_command: "gatk"
curl_command: "curl"
efetch_command: "efetch"
```

A safe example setup file is included at `config/preprocessing_hpc_env.sh` and is referenced by the default config. Edit that file for your cluster before running steps that require modules or a conda environment. For example:

```bash
module load R samtools wget
# or activate a project environment:
# source /path/to/miniconda3/etc/profile.d/conda.sh
# conda activate primate-mtdna
```

If your cluster provides tools under non-standard names or absolute paths, set the matching `*_command` key instead of relying on `PATH`.

## Step 0. Reference discovery

`preprocessing/scripts/run_reference_discovery.sh` runs `find_primate_wg_chrM_refs.py` on `data/metadata/all_species_list.txt`, a local RefSeq mitochondrion FASTA, and a primate tree. The raw output is written to `results/preprocessing/reference_discovery/species_reference_chrM_summary.tsv`. After manual review, copy or symlink the stable manifest to `data/metadata/species_reference_chrM_summary.tsv`.

For faster discovery on HPC, set `reference_discovery_threads` in `config/preprocessing_paths.yaml` to analyze multiple target species concurrently. NCBI E-utility calls remain rate-limited by the configured `delay`, so increase thread count cautiously and keep a valid NCBI email/API key when using higher concurrency.

## Step 1. Reference materialization

`preprocessing/scripts/build_reference_materialization_manifest.R` reads the reviewed discovery manifest and writes `reference_materialization_manifest.tsv` to both `results/preprocessing/reference_materialization/` and `references/manifests/`. It classifies chrM as `embedded_in_wg_ref`, `independent_chrM_ref`, or `missing_chrM_ref` and constructs expected local FASTA paths.

Reference discovery checks NCBI assembly reports for the selected WG candidates and their paired GCA/GCF assemblies, so embedded chrM rows should already prefer the partner that actually carries chrM before materialization begins.

The discovery step also writes `assembly_chrM_diagnostics.tsv` in the discovery output directory, listing every WG assembly report checked (including paired GCA/GCF partners) and the observed chrM status/contig/length. Use it with `species_reference_chrM_summary.tsv` to see whether the paired GCA was checked, whether its assembly report was unavailable, or whether chrM was rejected by length.

`preprocessing/scripts/materialize_references.sh` downloads WG FASTA files, stores NCBI `genomic.fna.gz` downloads as decompressed `.genome.fa` files so `samtools faidx` can index them, extracts embedded chrM records from the matching WG FASTA, falls back from a non-chrM-bearing GCA/GCF assembly to its paired NCBI assembly before retrying embedded extraction, downloads or extracts independent chrM FASTA files, and writes `references/manifests/in_house_score_reference_inputs.tsv`.

## Step 2. In-house score and minimal NUMT mask selection

The in-house score scripts consume `references/manifests/in_house_score_reference_inputs.tsv`. This is a reference-level step: it does not inspect sample CRAM files for NUMT detection. Instead, NUMT candidates are generated by BLASTing each selected chrM FASTA against the selected WG FASTA in the same manifest row. Embedded chrM rows must use chrM FASTA extracted from the WG reference; independent chrM rows must use the independent chrM FASTA. The merged result should be `results/preprocessing/in_house_score/merged_in_house_score.tsv`.

## Step 3. Build variant-calling reference packages

After `in_house_score`, build the reference packages used by downstream mitochondrial variant calling:

```bash
bash preprocessing/scripts/run_preprocessing.sh variant_references config/preprocessing_paths.yaml
```

This is still a reference-level step and does **not** require `data/metadata/sample_metadata.tsv` or any sample CRAM/CRAI paths. It consumes `references/manifests/in_house_score_reference_inputs.tsv` and `results/preprocessing/in_house_score/merged_in_house_score.tsv`, then writes per-reference FASTA packages and a manifest under `references/variant_calling/` by default. Configure the destination with `variant_calling_reference_out_root` in `config/preprocessing_paths.yaml`.

The generated reference package includes the whole-genome FASTA with configured NUMT masking applied where eligible, the chrM FASTA, shifted chrM FASTA, FASTA indexes, sequence dictionaries, BWA indexes, interval files, shift-back chain files, and `variant_calling_reference_manifest.tsv`.
