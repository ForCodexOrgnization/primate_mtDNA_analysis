# QC results

Generated QC preprocessing outputs should be written under this directory.

The default `config/qc_preprocessing.yaml` now places collected variant-calling
inputs in `results/qc/collected_variant_calling_results/` and coordinate
liftover outputs in `results/qc/coordinate_liftover/`. These generated result
subdirectories are created by the QC scripts when needed and should not be
committed unless they are intentionally curated small artifacts.

`results/qc/intraspecies_contamination/` contains original-coordinate core
variant tables, intra-species mixture-evidence tables, plots, and run logs. It
is independent of coordinate liftover.

## Downstream annotation outputs

`codon_match`, `trna_match`, and `rrna_match` write annotation-only VCFs and reports under their corresponding directories. They retain every input VCF record and run after coordinate liftover.
