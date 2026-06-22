# Preprocessing to variant calling

The final preprocessing output is `results/preprocessing/variant_calling_inputs/variant_calling_input_table.tsv`, one row per sample. It contains CRAM paths, WG FASTA and index, chrM FASTA and index, reference pairing status, chrM context, final reference strategy, minimal NUMT mask status, sample QC status, and manual-review fields.

Downstream mtDNA variant calling should use this table as its main input rather than rebuilding reference choices.
