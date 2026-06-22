#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 2 ]]; then echo "usage: $0 ACCESSION OUT_FASTA [LOCAL_MITO_FASTA]" >&2; exit 2; fi
ACC=$1; OUT=$2; LOCAL=${3:-}
MIN_LEN=${MIN_MITO_LEN:-14000}; MAX_LEN=${MAX_MITO_LEN:-25000}
mkdir -p "$(dirname "$OUT")"
if [[ -n "$LOCAL" && -f "$LOCAL" ]]; then
  samtools faidx "$LOCAL" || true
  if [[ -f "$LOCAL.fai" ]] && awk -v a="$ACC" '$1==a{found=1} END{exit !found}' "$LOCAL.fai"; then samtools faidx "$LOCAL" "$ACC" > "$OUT"; fi
fi
if [[ ! -s "$OUT" ]]; then
  if command -v efetch >/dev/null 2>&1; then efetch -db nucleotide -id "$ACC" -format fasta > "$OUT"; else curl -L "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ACC}&rettype=fasta&retmode=text" > "$OUT"; fi
fi
samtools faidx "$OUT"
len=$(awk 'NR==1{print $2}' "$OUT.fai")
[[ "$len" -ge "$MIN_LEN" && "$len" -le "$MAX_LEN" ]] || { echo -e "failure\tlength_outside_${MIN_LEN}_${MAX_LEN}\t$len" >&2; exit 1; }
echo -e "success\t$ACC\t$len"
