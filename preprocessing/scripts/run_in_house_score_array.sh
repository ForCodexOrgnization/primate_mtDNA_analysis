#!/usr/bin/env bash
set -euo pipefail
REF_INPUTS=${REF_INPUTS:-references/manifests/in_house_score_reference_inputs.tsv}
echo "Submit/iterate in-house score jobs using ${REF_INPUTS}; do not mix embedded and independent chrM FASTA paths." >&2
