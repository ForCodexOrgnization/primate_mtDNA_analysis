# Reference discovery inputs and outputs

Inputs: species table with `species` plus optional `sample_count`, `preprint_REFERENCE_SPECIES`, and `REFERENCE_SPECIES`; local `mitochondrion.1.1.genomic.fna.gz`; and `primate_tree.nwk`.

Outputs are written under `results/preprocessing/reference_discovery/`: `species_reference_chrM_summary.tsv`, status counts, candidate WG references, nuccore mitochondrial hits, and a cache directory. The reviewed stable manifest belongs at `data/metadata/species_reference_chrM_summary.tsv` and is the only discovery table consumed downstream.
