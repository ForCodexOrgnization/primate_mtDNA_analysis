#!/usr/bin/env bash
set -euo pipefail
SPECIES_TABLE=${SPECIES_TABLE:-data/metadata/all_species_list.txt}
MITO_FASTA=${MITO_FASTA:-/path/to/mitochondrion.1.1.genomic.fna.gz}
TREE_NEWICK=${TREE_NEWICK:-/path/to/primate_tree.nwk}
OUTDIR=${OUTDIR:-results/preprocessing/reference_discovery}
EMAIL=${EMAIL:-your_email@yale.edu}
MAX_NEAREST=${MAX_NEAREST:-200}
DELAY=${DELAY:-0.34}
mkdir -p "$OUTDIR"
python3 preprocessing/scripts/find_primate_wg_chrM_refs.py \
  --species "$SPECIES_TABLE" \
  --mito-fasta "$MITO_FASTA" \
  --tree "$TREE_NEWICK" \
  --outdir "$OUTDIR" \
  --email "$EMAIL" \
  --max-nearest "$MAX_NEAREST" \
  --delay "$DELAY"
