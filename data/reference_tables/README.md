# Reference tables

This directory contains small static reference inputs used by QC workflows.

- `human_chrM.fa` is the configured raw human mitochondrial reference FASTA for
  coordinate liftover (`config/qc_preprocessing.yaml`). It reserves the canonical
  rCRS/NC_012920.1 length of 16,569 bp and should be replaced with the exact
  NC_012920.1 base sequence before running reference-sensitive liftover on real
  VCF/COV inputs.
