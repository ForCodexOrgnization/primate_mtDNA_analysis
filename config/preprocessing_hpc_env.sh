#!/usr/bin/env bash
# Optional HPC environment setup for preprocessing.
#
# This file is sourced by preprocessing/scripts/run_preprocessing.sh when
# environment_setup_script points here in config/preprocessing_paths.yaml.
# It is intentionally safe by default: edit it for your cluster before running
# steps that require Rscript, samtools, wget/curl, or NCBI EDirect efetch.
#
# Examples -- uncomment and adapt to your HPC site:
#
module load R/4.4.2
module load SAMtools/1.21-GCC-13.3.0
module load BWA/0.7.18-GCCcore-13.3.0
module load GATK/4.6.1.0-GCCcore-13.3.0-Java-17
#
# Or activate a conda/mamba environment:
# source /path/to/miniconda3/etc/profile.d/conda.sh
# conda activate primate-mtdna
#
# If your executables have non-standard names or absolute paths, prefer setting
# the corresponding *_command key in config/preprocessing_paths.yaml.
