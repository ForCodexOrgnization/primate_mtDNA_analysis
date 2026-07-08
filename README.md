# Primate mtDNA Analysis

This repository keeps the active primate mtDNA preprocessing, QC, and variant-calling code together while avoiding unused placeholder files and duplicate output roots.

## Repository structure

```text
preprocessing/          # Reference discovery/materialization and variant-calling input preparation
variant_calling/        # Production WDL/Nextflow variant-calling workflows
qc_analysis/            # QC preprocessing utilities and coordinate liftover scripts
config/                 # Active workflow configuration files
results/                # Regenerated reports, tables, figures, and QC result products
data/                   # Curated small data files when present
```

Generated artifacts should be written under `results/`. The old `outputs/` root has been removed because it duplicated `results/`; QC preprocessing now writes under `results/qc/` by default.
