#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 3 ]]; then echo "usage: $0 WG_FASTA OUT_CHRM_FASTA CANDIDATE..." >&2; exit 2; fi
WG_FASTA=$1; OUT=$2; shift 2
MIN_LEN=${MIN_MITO_LEN:-14000}; MAX_LEN=${MAX_MITO_LEN:-25000}
mkdir -p "$(dirname "$OUT")"
samtools faidx "$WG_FASTA"
match=""
for c in "$@"; do
  [[ -z "${c:-}" ]] && continue
  if awk -v c="$c" 'BEGIN{ok=1} $1==c{ok=0} END{exit ok}' "$WG_FASTA.fai"; then match="$c"; break; fi
  m=$(awk -v c="$c" '$1 ~ "^" c "(\\.|$)" {print $1; exit}' "$WG_FASTA.fai")
  [[ -n "$m" ]] && { match="$m"; break; }
done
[[ -n "$match" ]] || { echo -e "failure\tno_candidate_found" >&2; exit 1; }
samtools faidx "$WG_FASTA" "$match" > "$OUT"
samtools faidx "$OUT"
len=$(awk 'NR==1{print $2}' "$OUT.fai")
[[ "$len" -ge "$MIN_LEN" && "$len" -le "$MAX_LEN" ]] || { echo -e "failure\tlength_outside_${MIN_LEN}_${MAX_LEN}\t$len" >&2; exit 1; }
echo -e "success\t$match\t$len"
