#!/usr/bin/env Rscript

# Render the human contamination filtering QC report from qc_analysis/.

report <- file.path("analysis", "05_human_contamination_filtering.qmd")
status <- system2("quarto", c("render", report))
quit(save = "no", status = status)
