# Primate mtDNA Analysis

This repository is organized to keep the production variant-calling workflow and downstream reproducible analysis in the same project while maintaining a clear logical separation between them.

## Repository structure

```text
variant_calling/        # Production Nextflow workflow
downstream_analysis/    # Quarto/R-based reproducible reports
config/                 # Global parameters, thresholds, paths, and sample metadata
data/                   # Small metadata, frozen inputs, and toy data
results/                # Regenerated reports, tables, and figures
docs/                   # User manuals and file documentation
```

The directories are scaffolded first and can be populated incrementally as the workflow, reports, configuration files, data assets, outputs, and documentation are developed.
