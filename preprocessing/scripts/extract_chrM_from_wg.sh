#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 3 ]]; then echo "usage: $0 WG_FASTA OUT_CHRM_FASTA CANDIDATE..." >&2; exit 2; fi
WG_FASTA=$1; OUT=$2; shift 2
MIN_LEN=${MIN_MITO_LEN:-14000}; MAX_LEN=${MAX_MITO_LEN:-25000}
SAMTOOLS_COMMAND=${SAMTOOLS_COMMAND:-samtools}
mkdir -p "$(dirname "$OUT")"
"$SAMTOOLS_COMMAND" faidx "$WG_FASTA"

try_extract() {
  local seq_name=$1 len
  [[ -z "${seq_name:-}" ]] && return 1
  if ! "$SAMTOOLS_COMMAND" faidx "$WG_FASTA" "$seq_name" > "$OUT"; then
    rm -f "$OUT" "$OUT.fai"
    echo -e "warning\tcandidate_extract_failed\t${seq_name}" >&2
    return 1
  fi
  if ! "$SAMTOOLS_COMMAND" faidx "$OUT"; then
    rm -f "$OUT" "$OUT.fai"
    echo -e "warning\tcandidate_output_index_failed\t${seq_name}" >&2
    return 1
  fi
  len=$(awk 'NR==1{print $2}' "$OUT.fai")
  if [[ "$len" =~ ^[0-9]+$ && "$len" -ge "$MIN_LEN" && "$len" -le "$MAX_LEN" ]]; then
    echo -e "success\t$seq_name\t$len"
    return 0
  fi
  rm -f "$OUT" "$OUT.fai"
  echo -e "warning\tcandidate_length_outside_${MIN_LEN}_${MAX_LEN}\t${seq_name}\t${len:-missing}" >&2
  return 1
}

tried_file=$(mktemp)
cleanup() { rm -f "$tried_file"; }
trap cleanup EXIT
try_once() {
  local seq_name=$1
  [[ -z "${seq_name:-}" ]] && return 1
  if awk -v s="$seq_name" '$0 == s {found=1} END{exit found ? 0 : 1}' "$tried_file"; then
    return 1
  fi
  printf '%s\n' "$seq_name" >> "$tried_file"
  try_extract "$seq_name"
}

for c in "$@"; do
  [[ -z "${c:-}" ]] && continue
  if awk -v c="$c" 'BEGIN{ok=1} $1==c{ok=0} END{exit ok}' "$WG_FASTA.fai"; then
    try_once "$c" && exit 0
  fi
  while IFS= read -r m; do
    try_once "$m" && exit 0
  done < <(awk -v c="$c" '$1 ~ "^" c "(\\.|$)" {print $1}' "$WG_FASTA.fai")
done

# Do not rely only on manifest/assembly-report accessions.  Some NCBI GCA/GCF
# partners use mitochondrial names such as chrMT/MT/M in the local FASTA even
# when the expected accession fields differ, so scan the actual local FASTA
# headers and keep only complete-mitogenome-sized records.
while IFS=$'\t' read -r _score seq_name; do
  try_once "$seq_name" && exit 0
done < <(
  awk -v min="$MIN_LEN" -v max="$MAX_LEN" '
    function emit() {
      if (name != "" && len >= min && len <= max && is_mito(header, name)) {
        print mito_score(header, name) "\t" name
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
    END { emit() }
  ' "$WG_FASTA" | sort -k1,1nr
)

echo -e "failure\tno_extractable_complete_mito_candidate_found" >&2
exit 1
