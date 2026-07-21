# MITOS2 annotation integration

MITOS2 is invoked only through its conda environment, not as a presumed system command:

```bash
module load miniconda/24.11.3
conda activate mitos2
```

`run_mitos2_annotation.py` activates that environment in a login shell, probes `mitos2`, `mitos`, then `runmitos.py`, and records the selected executable in `results/qc/mitos2_annotation/mitos2_annotation_summary.tsv`.

The workflow deduplicates final materialized chrM FASTAs from the resolved reference manifest. Those FASTAs are the coordinate truth for all emitted positions. MITOS2 supplies CDS, tRNA, and rRNA *intervals*; its tRNA/rRNA output does not provide secondary-structure stem/loop information and does not replace tRNAscan paired-site annotations or human-guided rRNA stem/loop annotation. The interval table may support future fallback region tables.

The first integration target is CDS/codon fallback. `build_primate_codon_table.py` selects a single annotation source for each sample: valid nonzero GenBank CDS rows first, otherwise MITOS2 rows for that final reference. It never combines sources within a sample. Raw/parsed MITOS2 tables remain separate, and the source comparison is written to `results/qc/codon_table_build/genbank_vs_mitos2_cds_comparison.tsv`. The chosen table used downstream is `data/reference_tables/all_primate_position_codon_table.tsv`.

```bash
bash qc_analysis/scripts/run_qc_preprocessing.sh mitos2_annotation config/qc_preprocessing.yaml
bash qc_analysis/scripts/run_qc_preprocessing.sh build_primate_codon_table config/qc_preprocessing.yaml
bash qc_analysis/scripts/run_qc_preprocessing.sh all config/qc_preprocessing.yaml
```
