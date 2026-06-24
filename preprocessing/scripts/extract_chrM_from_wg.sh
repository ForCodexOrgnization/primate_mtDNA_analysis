#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 3 ]]; then echo "usage: $0 WG_FASTA OUT_CHRM_FASTA CANDIDATE..." >&2; exit 2; fi
WG_FASTA=$1; OUT=$2; shift 2
MIN_LEN=${MIN_MITO_LEN:-14000}; MAX_LEN=${MAX_MITO_LEN:-25000}
SAMTOOLS_COMMAND=${SAMTOOLS_COMMAND:-samtools}
mkdir -p "$(dirname "$OUT")"
"$SAMTOOLS_COMMAND" faidx "$WG_FASTA"
match=""
for c in "$@"; do
  [[ -z "${c:-}" ]] && continue
  if awk -v c="$c" 'BEGIN{ok=1} $1==c{ok=0} END{exit ok}' "$WG_FASTA.fai"; then match="$c"; break; fi
  m=$(awk -v c="$c" '$1 ~ "^" c "(\\.|$)" {print $1; exit}' "$WG_FASTA.fai")
  [[ -n "$m" ]] && { match="$m"; break; }
done
if [[ -z "$match" ]]; then
  # Do not rely only on manifest/assembly-report accessions.  Some NCBI GCA/GCF
  # partners use mitochondrial names such as chrMT/MT/M in the local FASTA even
  # when the expected accession fields differ, so scan the actual local FASTA
  # headers and keep only complete-mitogenome-sized records.
  match=$(
    awk -v min="$MIN_LEN" -v max="$MAX_LEN" '
      function emit() {
        if (name != "" && len >= min && len <= max && is_mito(header, name)) {
          score = mito_score(header, name)
          if (score > best_score) {
            best_score = score
            best_name = name
          }
        }
      }
      function is_mito(h, n, lower_h, lower_n) {
        lower_h = tolower(h)
        lower_n = tolower(n)
        return (lower_n ~ /^(chrm|chrmt|mt|m)$/ || lower_n ~ /(^|[_|.:-])(chrm|chrmt|mt|mitochondrion|mitochondrial)([_|.:-]|$)/ || lower_h ~ /(^|[[:space:]_.,;:|()-])(chrm|chrmt|mt|mitochondrion|mitochondrial|mitochondria)([[:space:]_.,;:|()-]|$)/)
      }
      function mito_score(h, n, lower_h, lower_n, s) {
        lower_h = tolower(h)
        lower_n = tolower(n)
        s = 0
        if (lower_n == "chrm" || lower_n == "chrmt" || lower_n == "mt" || lower_n == "m") s += 100
        if (lower_h ~ /complete[[:space:]_-]+mitochond/) s += 50
        if (lower_h ~ /mitochond/) s += 20
        return s
      }
      /^>/ {
        emit()
        header = substr($0, 2)
        split(header, parts, /[[:space:]]+/)
        name = parts[1]
        len = 0
        next
      }
      {
        gsub(/[[:space:]]/, "", $0)
        len += length($0)
      }
      END {
        emit()
        if (best_name != "") print best_name
      }
    ' "$WG_FASTA"
  )
fi
[[ -n "$match" ]] || { echo -e "failure\tno_candidate_found" >&2; exit 1; }
"$SAMTOOLS_COMMAND" faidx "$WG_FASTA" "$match" > "$OUT"
"$SAMTOOLS_COMMAND" faidx "$OUT"
len=$(awk 'NR==1{print $2}' "$OUT.fai")
[[ "$len" -ge "$MIN_LEN" && "$len" -le "$MAX_LEN" ]] || { echo -e "failure\tlength_outside_${MIN_LEN}_${MAX_LEN}\t$len" >&2; exit 1; }
echo -e "success\t$match\t$len"
