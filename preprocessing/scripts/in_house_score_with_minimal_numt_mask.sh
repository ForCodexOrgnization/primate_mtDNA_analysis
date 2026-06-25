#!/bin/bash
#SBATCH --job-name=merge_ref
#SBATCH --output=log/preprocessing/numt_%A_%a.out
#SBATCH --error=log/preprocessing/numt_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00

set -euo pipefail

# =============================================================================
# In-house reference classifier + reference-level NUMT mask candidate generator
#
# Main fixes in this version:
#   1. Species list is de-duplicated with sort -u, avoiding repeated rows in the
#      merged summary.
#   2. Minimal chrM-covering set is calculated using the UNION of true query
#      intervals per merged subject fragment, not min(qstart)~max(qend).
#      This prevents one short fragment from falsely covering the whole chrM.
#   3. Summary header uses Bash $'...' so tabs are real tabs, not literal \t.
#   4. Mask BED/fragment TSV files are written for #C-Ambiguous and #C-likely_comp references.
#
# Recommended use:
#   SUBMIT_ARRAY=1 bash preprocessing/scripts/in_house_score_with_minimal_numt_mask.sh
#
# Or submit manually after checking the manifest row count:
#   N=$(python3 - <<'PY_N'
# import csv
# seen = set()
# with open("references/manifests/in_house_score_reference_inputs.tsv", newline="") as h:
#     for r in csv.DictReader(h, delimiter="\t"):
#         key = tuple((r.get(c) or "").strip()
#                     for c in ("target_species", "wg_fasta_path", "chrM_fasta_path"))
#         if all(key):
#             seen.add(key)
# print(len(seen))
# PY_N
#   )
#   jid=$(sbatch --array=1-${N}%50 preprocessing/scripts/in_house_score_with_minimal_numt_mask.sh | awk '{print $4}')
#   sbatch --dependency=afterok:${jid} --export=ALL,MERGE_ONLY=1 preprocessing/scripts/in_house_score_with_minimal_numt_mask.sh
#
# Or after array finishes:
#   MERGE_ONLY=1 bash preprocessing/scripts/in_house_score_with_minimal_numt_mask.sh
# =============================================================================

# -------------------- Config --------------------
# Reference-level inputs produced by materialize_references.sh. This step does
# not read CRAM/sample files; NUMT candidates are detected by BLASTing each chrM
# FASTA against its selected whole-genome reference FASTA.
REF_INPUTS="${REF_INPUTS:-references/manifests/in_house_score_reference_inputs.tsv}"
OUTDIR="${OUTDIR:-results/preprocessing/in_house_score}"
MERGED_IN_HOUSE_SCORE="${MERGED_IN_HOUSE_SCORE:-${OUTDIR}/merged_in_house_score.tsv}"
PYTHON_COMMAND="${PYTHON_COMMAND:-python3}"
MAKEBLASTDB_COMMAND="${MAKEBLASTDB_COMMAND:-makeblastdb}"
BLASTN_COMMAND="${BLASTN_COMMAND:-blastn}"
ARRAY_CONCURRENCY="${ARRAY_CONCURRENCY:-50}"
IN_HOUSE_SCORE_LOG_DIR="${IN_HOUSE_SCORE_LOG_DIR:-log/preprocessing}"

mkdir -p "$OUTDIR" "$IN_HOUSE_SCORE_LOG_DIR" "${OUTDIR}/numt_candidates" "${OUTDIR}/numt_beds"
if [[ -n "${BLAST_MODULE:-}" ]] && command -v module >/dev/null 2>&1; then
  module load "$BLAST_MODULE"
fi

# Use Slurm value if available, otherwise safe local default.
THREADS="${SLURM_CPUS_PER_TASK:-4}"

# -------------------- Tunable classification thresholds --------------------
CHRMSIZE_RATIO_MIN="0.80"
CHRMSIZE_RATIO_MAX="1.20"
A_TOP_MERGED_RATIO_MIN="0.95"
CLEAN_TOP_MERGED_RATIO_MIN="0.90"
MAJOR_TOP_MERGED_RATIO_MIN="0.60"
MEGA_MERGED_LEN_MIN="10000"
FRAG_CUM_MERGED_RATIO_MIN="0.90"
FRAG_CONTIGS_MIN="3"
FRAG_INTERVALS_MIN="3"

# -------------------- Tunable NUMT mask-candidate thresholds --------------------
# These thresholds are intentionally moderate for first-pass summary generation.
# After reviewing the summary, tighten or loosen them and rerun if needed.
NUMT_MIN_PIDENT="90"
NUMT_MIN_ALN_LEN="100"
NUMT_MAX_EVALUE="1e-3"
NUMT_MIN_BITSCORE="0"
NUMT_TARGET_CHRM_COV="0.95"
NUMT_PAD_BP="50"
MASK_REF_TYPES="${MASK_REF_TYPES:-#C-likely_comp,#C-Ambiguous}"
A_MASK_MODE="${A_MASK_MODE:-diagnostic_only}"

# -------------------- Final mask selection --------------------
# #C-Ambiguous and #C-likely_comp references receive a minimal chrM-covering
# FINAL NUMT mask. Other reference classes keep header-only final mask files.

# -------------------- Summary header --------------------
BASE_SUMMARY_HEADER=$'Species\tREF_TYPE\tMTLIKE_PATTERN\tValidAnnotatedMitoContig\tValidAnnotatedMitoLength\tHasValidAnnotatedMito\tTopContig\tTopLength\tMitoRefLength\tMaxHitRatio\tTotalHitRatio\tContigRatio\tScore\tM\t_dCplus\t_dCminus\t_dTplus\t_dTminus\tTopContigMergedLen\tTopContigMergedRatio\tLongestMergedHitLen\tCumulativeMergedHitLength\tMergedHitContigsN\tMergedHitIntervalsN\tNonTopLongestMergedHitLen\tNonTopCumulativeMergedHitLength\tNonTopMergedHitContigsN\tNonTopMergedHitIntervalsN'
NUMT_SUMMARY_HEADER=$'AllHighConfFragmentsN\tAllHighConfContigsN\tAllHighConfSubjectLen\tAllHighConfLongestFragment\tAllHighConfChrMCoveredBp\tAllHighConfChrMCoverageRatio\tMinimalFragmentsN\tMinimalContigsN\tMinimalSubjectLen\tMinimalLongestFragment\tMinimalChrMCoveredBp\tMinimalChrMCoverageRatio\tMinimalTargetReached\tMinimalVsAllFragmentRatio\tMinimalVsAllSubjectLenRatio\tMaskPriority\tMinimalMaskBED\tFullMaskBED\tCandidateTSV'
A_NUMT_SUMMARY_HEADER=$'A_NonChrM_HighConfFragmentsN\tA_NonChrM_HighConfContigsN\tA_NonChrM_SubjectLen\tA_NonChrM_LongestFragment\tA_NonChrM_ChrMCoveredBp\tA_NonChrM_ChrMCoverageRatio\tA_NonChrM_MinimalFragmentsN\tA_NonChrM_MinimalContigsN\tA_NonChrM_MinimalSubjectLen\tA_NonChrM_MinimalChrMCoverageRatio\tA_NonChrM_MinimalTargetReached'
SUMMARY_HEADER="${BASE_SUMMARY_HEADER}"$'\t'"${NUMT_SUMMARY_HEADER}"$'\t'"${A_NUMT_SUMMARY_HEADER}"

# -------------------- Helpers --------------------
ts() { date "+%F %T"; }
log(){ echo "[$(ts)] [INFO] $*" >&2; }
warn(){ echo "[$(ts)] [WARN] $*" >&2; }
err(){ echo "[$(ts)] [ERROR] $*" >&2; }
safe_id() {
  printf "%s" "$1" | tr ' /' '__' | tr -cd '[:alnum:]_.-'
}

acquire_merge_lock() {
  local lockdir="$1"
  local waited=0
  while ! mkdir "$lockdir" 2>/dev/null; do
    if (( waited >= 3600 )); then
      err "Timed out waiting for merge lock: ${lockdir}"
      exit 1
    fi
    sleep 5
    waited=$((waited + 5))
  done
}

release_merge_lock() {
  local lockdir="$1"
  rmdir "$lockdir" 2>/dev/null || true
}


write_failure_summary() {
  local summary_file="$1" species="$2" ref_type="$3" pattern="$4" annotated_contig="$5" annotated_len="$6" has_annotated="$7" mito_len="$8" priority="$9"
  "$PYTHON_COMMAND" - "$SUMMARY_HEADER" "$species" "$ref_type" "$pattern" "$annotated_contig" "$annotated_len" "$has_annotated" "$mito_len" "$priority" <<'PY_FAIL_SUMMARY'
import sys
header, species, ref_type, pattern, annot, annot_len, has_annot, mito_len, priority = sys.argv[1:]
cols = header.split("\t")
row = {c: "0" for c in cols}
row.update({
    "Species": species, "REF_TYPE": ref_type, "MTLIKE_PATTERN": pattern,
    "ValidAnnotatedMitoContig": annot, "ValidAnnotatedMitoLength": annot_len,
    "HasValidAnnotatedMito": has_annot, "TopContig": "NA", "TopLength": "0",
    "MitoRefLength": mito_len, "MinimalTargetReached": "no", "MaskPriority": priority,
    "MinimalMaskBED": "NA", "FullMaskBED": "NA", "CandidateTSV": "NA",
    "A_NonChrM_MinimalTargetReached": "no",
})
print("\t".join(row.get(c, "0") for c in cols))
PY_FAIL_SUMMARY
}

append_valid_summary_body() {
  local summary_file="$1" merged_file="$2" expected_cols="$3"
  local cleaned_body
  cleaned_body="$(mktemp "${OUTDIR}/.summary_body.XXXXXX")"
  if LC_ALL=C tr -d '\000' < "$summary_file" | awk -F'\t' -v n="$expected_cols" '
      NR == 1 { next }
      NF == n { print; next }
      { bad++ }
      END { if (bad) exit 2 }
    ' > "$cleaned_body"; then
    cat "$cleaned_body" >> "$merged_file"
  else
    warn "Skipping malformed summary ${summary_file}; expected every data row to have ${expected_cols} tab-delimited columns. Re-run its array task before trusting the merged table."
  fi
  rm -f "$cleaned_body"
}

fasta_header_full() {
  local sp="$1" contig="$2"
  local norm="${sp// /_}"
  printf ">%s isolate=JT067 genome=assembly, contig: %s, whole genome shotgun sequence\n" "$norm" "$contig"
}

format_sequence() { fold -w 80 | sed 's/[[:space:]]*$//'; }

# ---- score calculator ----
# Keeps the current in-house score formula unchanged.
calc_score_py() {
  local max_hit_ratio="$1"
  local total_hit_ratio="$2"
  local contig_ratio="$3"

  "$PYTHON_COMMAND" - "$max_hit_ratio" "$total_hit_ratio" "$contig_ratio" <<'PY'
import math, sys
MHR = float(sys.argv[1]); THR = float(sys.argv[2]); CR = float(sys.argv[3])
M = max(0.0, min(1.0, MHR))
def safe_log2(x):
    if x <= 0: return None
    return math.log(x, 2)
lc = safe_log2(CR)
if lc is None:
    dC_plus = 0.0; dC_minus = 50.0
else:
    dC_plus = max(0.0, lc); dC_minus = max(0.0, -lc)
lt = safe_log2(THR)
if lt is None:
    dT_plus = 0.0; dT_minus = 50.0
else:
    dT_plus = max(0.0, lt); dT_minus = max(0.0, -lt)
Kc_plus = 1.5; Kc_minus = 1.5; Kt_plus = 0.8; Kt_minus = 0.3
score = M * math.exp(-(Kc_plus*dC_plus + Kc_minus*dC_minus + Kt_plus*dT_plus + Kt_minus*dT_minus))
print(f"{score:.6g}\t{M:.6g}\t{dC_plus:.6g}\t{dC_minus:.6g}\t{dT_plus:.6g}\t{dT_minus:.6g}")
PY
}

# ---- select a valid annotated mitochondrial contig ----
select_annotated_mito_contig_py() {
  local genome="$1"
  local mito_len="$2"

  "$PYTHON_COMMAND" - "$genome" "$mito_len" "$CHRMSIZE_RATIO_MIN" "$CHRMSIZE_RATIO_MAX" <<'PY'
import sys

genome = sys.argv[1]
mito_len = int(float(sys.argv[2])) if sys.argv[2] not in ("", "NA") else 0
ratio_min = float(sys.argv[3]); ratio_max = float(sys.argv[4])
candidates = []
name = None; header = None; seq_len = 0

def flush():
    global name, header, seq_len
    if name is None:
        return
    h = header.lower()
    tokens = set(h.replace('|', ' ').replace(',', ' ').replace(';', ' ').split())
    if ("mitochondr" in h) or ("mitochondrial" in h) or ("chrm" in h) or ("mtdna" in h) or ("mt" in tokens):
        candidates.append((name, seq_len, header))

with open(genome) as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith(">"):
            flush()
            header = line[1:]
            name = header.split()[0]
            seq_len = 0
        else:
            seq_len += len(line.strip())
flush()

if mito_len <= 0 or not candidates:
    print("NA\t0\t0"); sys.exit(0)
valid = []
for c, l, h in candidates:
    ratio = l / mito_len if mito_len > 0 else 0
    if ratio_min <= ratio <= ratio_max:
        valid.append((abs(l - mito_len), c, l))
if not valid:
    print("NA\t0\t0")
else:
    _, c, l = sorted(valid)[0]
    print(f"{c}\t{l}\t1")
PY
}

# ---- merge BLAST subject intervals and summarize pre-NUMT-like signal ----
merge_blast_stats_py() {
  local blast_file="$1"
  local annotated_contig="$2"
  local mito_len="$3"

  "$PYTHON_COMMAND" - "$blast_file" "$annotated_contig" "$mito_len" <<'PY'
import sys
from collections import defaultdict
blast_file = sys.argv[1]; annotated = sys.argv[2]
mito_len = int(float(sys.argv[3])) if sys.argv[3] not in ("", "NA") else 0
intervals = defaultdict(list)
with open(blast_file) as f:
    for line in f:
        if not line.strip(): continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 7: fields = line.split()
        if len(fields) < 7: continue
        contig = fields[1]
        sstart = int(float(fields[4])); send = int(float(fields[5]))
        a, b = sorted((sstart, send))
        intervals[contig].append((a, b))

def merge_ints(ints):
    if not ints: return []
    ints = sorted(ints)
    merged = [list(ints[0])]
    for a, b in ints[1:]:
        if a <= merged[-1][1] + 1:
            if b > merged[-1][1]: merged[-1][1] = b
        else:
            merged.append([a, b])
    return [(a, b) for a, b in merged]

merged_by_contig = {c: merge_ints(v) for c, v in intervals.items()}
merged_len_by_contig = {c: sum(b-a+1 for a,b in ints) for c, ints in merged_by_contig.items()}
longest_interval_by_contig = {c: max((b-a+1 for a,b in ints), default=0) for c, ints in merged_by_contig.items()}
if not merged_len_by_contig:
    print("NoHit\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0"); sys.exit(0)
if annotated and annotated != "NA" and annotated in merged_len_by_contig:
    top_contig = annotated
else:
    top_contig = sorted(merged_len_by_contig, key=lambda c: (merged_len_by_contig[c], longest_interval_by_contig[c]), reverse=True)[0]
top_merged_len = merged_len_by_contig.get(top_contig, 0)
top_ratio = top_merged_len / mito_len if mito_len > 0 else 0.0
longest_merged_hit = max(longest_interval_by_contig.values())
cum_merged_len = sum(merged_len_by_contig.values())
contigs_n = sum(1 for v in merged_len_by_contig.values() if v > 0)
intervals_n = sum(len(v) for v in merged_by_contig.values())
non_top_contigs = [c for c in merged_len_by_contig if c != top_contig]
non_top_longest = max((longest_interval_by_contig[c] for c in non_top_contigs), default=0)
non_top_cum = sum(merged_len_by_contig[c] for c in non_top_contigs)
non_top_contigs_n = sum(1 for c in non_top_contigs if merged_len_by_contig[c] > 0)
non_top_intervals_n = sum(len(merged_by_contig[c]) for c in non_top_contigs)
print(f"{top_contig}\t{top_merged_len}\t{top_ratio:.6g}\t{longest_merged_hit}\t{cum_merged_len}\t{contigs_n}\t{intervals_n}\t{non_top_longest}\t{non_top_cum}\t{non_top_contigs_n}\t{non_top_intervals_n}")
PY
}

# ---- REF_TYPE classifier ----
classify_ref_py() {
  local has_valid_annotated_mito="$1"
  local top_contig_len="$2"
  local mito_len="$3"
  local top_merged_len="$4"
  local top_merged_ratio="$5"
  local longest_merged_hit="$6"
  local cumulative_merged_len="$7"
  local merged_contigs_n="$8"
  local merged_intervals_n="$9"
  local score="${10}"
  local mhr="${11}"
  local thr="${12}"
  local cr="${13}"

  "$PYTHON_COMMAND" - "$has_valid_annotated_mito" "$top_contig_len" "$mito_len" \
    "$top_merged_len" "$top_merged_ratio" "$longest_merged_hit" \
    "$cumulative_merged_len" "$merged_contigs_n" "$merged_intervals_n" \
    "$score" "$mhr" "$thr" "$cr" \
    "$CHRMSIZE_RATIO_MIN" "$CHRMSIZE_RATIO_MAX" "$A_TOP_MERGED_RATIO_MIN" \
    "$CLEAN_TOP_MERGED_RATIO_MIN" "$MAJOR_TOP_MERGED_RATIO_MIN" "$MEGA_MERGED_LEN_MIN" \
    "$FRAG_CUM_MERGED_RATIO_MIN" "$FRAG_CONTIGS_MIN" "$FRAG_INTERVALS_MIN" <<'PY'
import sys
has_valid_annot = int(float(sys.argv[1])); top_contig_len = float(sys.argv[2]); mito_len = float(sys.argv[3])
top_merged_len = float(sys.argv[4]); top_merged_ratio = float(sys.argv[5]); longest_merged_hit = float(sys.argv[6])
cum_merged_len = float(sys.argv[7]); merged_contigs_n = int(float(sys.argv[8])); merged_intervals_n = int(float(sys.argv[9]))
score = float(sys.argv[10]); mhr = float(sys.argv[11]); thr = float(sys.argv[12]); cr = float(sys.argv[13])
chrmsize_min = float(sys.argv[14]); chrmsize_max = float(sys.argv[15]); a_top_merged_min = float(sys.argv[16])
clean_top_merged_min = float(sys.argv[17]); major_top_merged_min = float(sys.argv[18]); mega_len_min = float(sys.argv[19])
frag_cum_ratio_min = float(sys.argv[20]); frag_contigs_min = int(float(sys.argv[21])); frag_intervals_min = int(float(sys.argv[22]))
top_len_ratio = top_contig_len / mito_len if mito_len > 0 else 0.0
cum_ratio = cum_merged_len / mito_len if mito_len > 0 else 0.0
valid_chrM_sized_contig = chrmsize_min <= top_len_ratio <= chrmsize_max
is_A = has_valid_annot == 1 and valid_chrM_sized_contig and top_merged_ratio >= a_top_merged_min
has_clean_complete_top_contig = valid_chrM_sized_contig and top_merged_ratio >= clean_top_merged_min
has_major_top_contig = top_merged_ratio >= major_top_merged_min or top_merged_len >= mega_len_min or longest_merged_hit >= mega_len_min
has_fragmented_cumulative_burden = (not has_major_top_contig) and cum_ratio >= frag_cum_ratio_min and merged_contigs_n >= frag_contigs_min and merged_intervals_n >= frag_intervals_min
if is_A:
    label = "#A"; reason = "annotated_complete_chrM"
elif has_clean_complete_top_contig:
    label = "#C-likely_comp"; reason = "clean_complete_top_contig"
elif has_major_top_contig:
    label = "#C-likely_comp"; reason = "mega_or_major_top_contig"
elif has_fragmented_cumulative_burden:
    label = "#C-Ambiguous"; reason = "fragmented_cumulative_burden"
else:
    label = "#C-likely_incomp"; reason = "weak_or_no_mtlike_signal"
print(f"{label}\t{reason}")
PY
}

# ---- reference-level NUMT candidate generation and minimal set cover ----
# Output one TSV summary line plus writes:
#   candidate TSV: all high-confidence merged reference-level mt-like fragments
#   full TSV/BED: all high-confidence fragments for mask-eligible REF_TYPE values
#   minimal TSV/BED: greedy minimal chrM-covering selected fragments for mask-eligible REF_TYPE values
numt_mask_summary_py() {
  local species="$1"
  local blast_file="$2"
  local mito_len="$3"
  local ref_type="$4"
  local species_id="${5:-$(safe_id "$species")}"
  local candidate_tsv="${OUTDIR}/numt_candidates/${species_id}.highconf_reference_mtlike_candidates.tsv"
  local full_tsv="${OUTDIR}/numt_beds/${species_id}.full_highconf.tsv"
  local minimal_tsv="${OUTDIR}/numt_beds/${species_id}.minimal_chrMcover.tsv"
  local full_bed="${OUTDIR}/numt_beds/${species_id}.full_highconf.bed"
  local minimal_bed="${OUTDIR}/numt_beds/${species_id}.minimal_chrMcover.bed"
  local final_tsv="${OUTDIR}/numt_beds/${species_id}.FINAL_numt_mask.tsv"
  local final_bed="${OUTDIR}/numt_beds/${species_id}.FINAL_numt_mask.bed"

  "$PYTHON_COMMAND" - "$species" "$blast_file" "$mito_len" "$ref_type" \
    "$candidate_tsv" "$full_tsv" "$minimal_tsv" "$full_bed" "$minimal_bed" "$final_tsv" "$final_bed" \
    "$NUMT_MIN_PIDENT" "$NUMT_MIN_ALN_LEN" "$NUMT_MAX_EVALUE" "$NUMT_MIN_BITSCORE" \
    "$NUMT_TARGET_CHRM_COV" "$NUMT_PAD_BP" "$ANNOTATED_CONTIG_FOR_NUMT" "$MASK_REF_TYPES" "$A_MASK_MODE" <<'PY'
import sys, math
from collections import defaultdict

species, blast_file, mito_len_s, ref_type = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
candidate_tsv, full_tsv, minimal_tsv, full_bed, minimal_bed, final_tsv, final_bed = sys.argv[5:12]
min_pident = float(sys.argv[12]); min_aln_len = int(float(sys.argv[13]))
max_evalue = float(sys.argv[14]); min_bitscore = float(sys.argv[15])
target_cov = float(sys.argv[16]); pad_bp = int(float(sys.argv[17]))
annotated_mito_contig = sys.argv[18]
mask_ref_types = {x.strip() for x in sys.argv[19].split(",") if x.strip()}
a_mask_mode = sys.argv[20]
if a_mask_mode not in ("diagnostic_only", "mask_if_requested"):
    raise SystemExit(f"Unsupported A_MASK_MODE={a_mask_mode!r}; expected diagnostic_only or mask_if_requested")
mito_len = int(float(mito_len_s)) if mito_len_s not in ("", "NA") else 0
mask_requested = ref_type in mask_ref_types
should_final_mask = mask_requested and (ref_type != "#A" or a_mask_mode == "mask_if_requested")

cand_header = "Species\tFragmentID\tSubjectContig\tSubjectStart1\tSubjectEnd1\tSubjectStart0\tSubjectEnd0\tSubjectLen\tQueryIntervals\tQueryCoveredBp\tPidentMax\tPidentMean\tBitscoreMax\tEvalueMin\tHSPsN"
bed_header = "#chrom\tstart0\tend0\tfragment_id\tspecies\tq_covered_bp\tsubject_len\tpident_max\tbitscore_max"

# Always initialize output files, so stale old results cannot be reused.
for path, header in [(candidate_tsv, cand_header), (full_tsv, cand_header), (minimal_tsv, cand_header), (final_tsv, cand_header)]:
    with open(path, "w") as out:
        out.write(header + "\n")
for path in [full_bed, minimal_bed, final_bed]:
    with open(path, "w") as out:
        out.write(bed_header + "\n")

if mito_len <= 0:
    empty = [0,0,0,0,0,"0",0,0,0,0,0,"0","no","NA","NA","not_evaluable",minimal_bed,full_bed,candidate_tsv]
    print("\t".join(map(str, empty)))
    sys.exit(0)

# Read filtered HSPs grouped by subject contig.
# BLAST outfmt 6 columns:
# qseqid sseqid qstart qend sstart send length pident mismatch evalue bitscore
hsps_by_contig = defaultdict(list)
try:
    fh = open(blast_file)
except FileNotFoundError:
    fh = []

for line in fh:
    if not line.strip():
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 11:
        f = line.split()
    if len(f) < 11:
        continue
    qstart, qend = int(float(f[2])), int(float(f[3]))
    sstart, send = int(float(f[4])), int(float(f[5]))
    aln_len = int(float(f[6])); pident = float(f[7]); evalue = float(f[9]); bitscore = float(f[10])
    if aln_len < min_aln_len or pident < min_pident or evalue > max_evalue or bitscore < min_bitscore:
        continue
    contig = f[1]
    if ref_type == "#A" and annotated_mito_contig not in ("", "NA") and contig == annotated_mito_contig:
        continue
    sa, sb = sorted((sstart, send))
    qa, qb = sorted((qstart, qend))
    # Clamp query coordinates to chrM length.
    qa = max(1, min(qa, mito_len)); qb = max(1, min(qb, mito_len))
    hsps_by_contig[contig].append({
        "sa": sa, "sb": sb, "qa": qa, "qb": qb,
        "pident": pident, "evalue": evalue, "bitscore": bitscore,
        "aln_len": aln_len,
    })
if hasattr(fh, "close"):
    fh.close()

# Merge HSPs into subject-side fragments. Crucial fix:
# each merged fragment retains a LIST of true query intervals; query coverage is the union of those intervals.
def merge_hsps_by_subject(hsps):
    if not hsps:
        return []
    hsps = sorted(hsps, key=lambda x: (x["sa"], x["sb"]))
    groups = []
    cur = [hsps[0]]
    cur_end = hsps[0]["sb"]
    for h in hsps[1:]:
        if h["sa"] <= cur_end + 1:
            cur.append(h)
            if h["sb"] > cur_end:
                cur_end = h["sb"]
        else:
            groups.append(cur)
            cur = [h]
            cur_end = h["sb"]
    groups.append(cur)
    return groups

def q_bases_from_intervals(qints):
    bases = set()
    for a, b in qints:
        a = max(1, min(a, mito_len)); b = max(1, min(b, mito_len))
        if a > b: a, b = b, a
        # 0-based half-open base set.
        bases.update(range(a-1, b))
    return bases

def merge_query_intervals(qints):
    if not qints:
        return []
    xs = sorted((min(a,b), max(a,b)) for a,b in qints)
    merged = [list(xs[0])]
    for a,b in xs[1:]:
        if a <= merged[-1][1] + 1:
            if b > merged[-1][1]: merged[-1][1] = b
        else:
            merged.append([a,b])
    return [(a,b) for a,b in merged]

def qints_to_str(qints):
    if not qints:
        return "NA"
    return ";".join(f"{a}-{b}" for a,b in qints)

def write_candidate(out, c):
    out.write("\t".join(map(str, [
        species, c["id"], c["contig"], c["s1"], c["s2"], c["start0"], c["end0"], c["subject_len"],
        qints_to_str(c["q_intervals_merged"]), c["q_covered_bp"],
        f'{c["pident_max"]:.6g}', f'{c["pident_mean"]:.6g}', f'{c["bitscore_max"]:.6g}', f'{c["evalue_min"]:.6g}', c["hsps_n"]
    ])) + "\n")

def write_bed(out, c):
    out.write("\t".join(map(str, [
        c["contig"], c["start0"], c["end0"], c["id"], species, c["q_covered_bp"], c["subject_len"],
        f'{c["pident_max"]:.6g}', f'{c["bitscore_max"]:.6g}'
    ])) + "\n")

candidates = []
frag_i = 0
for contig, hsps in hsps_by_contig.items():
    for group in merge_hsps_by_subject(hsps):
        frag_i += 1
        s1 = min(h["sa"] for h in group)
        s2 = max(h["sb"] for h in group)
        start0 = max(0, s1 - 1 - pad_bp)
        end0 = s2 + pad_bp
        subject_len = end0 - start0
        qints = [(h["qa"], h["qb"]) for h in group]
        qints_merged = merge_query_intervals(qints)
        qbases = q_bases_from_intervals(qints_merged)
        pidents = [h["pident"] for h in group]
        bitscores = [h["bitscore"] for h in group]
        evalues = [h["evalue"] for h in group]
        candidates.append({
            "id": f"{species}_frag{frag_i:06d}",
            "contig": contig,
            "s1": s1, "s2": s2,
            "start0": start0, "end0": end0,
            "subject_len": subject_len,
            "q_intervals_merged": qints_merged,
            "q_bases": qbases,
            "q_covered_bp": len(qbases),
            "pident_max": max(pidents),
            "pident_mean": sum(pidents)/len(pidents),
            "bitscore_max": max(bitscores),
            "evalue_min": min(evalues),
            "hsps_n": len(group),
        })

# Write all candidate fragments for all species.
with open(candidate_tsv, "a") as out:
    for c in sorted(candidates, key=lambda x: (x["contig"], x["start0"], x["end0"])):
        write_candidate(out, c)

all_frag_n = len(candidates)
all_contigs_n = len(set(c["contig"] for c in candidates)) if candidates else 0
all_subject_len = sum(c["subject_len"] for c in candidates)
all_longest = max((c["subject_len"] for c in candidates), default=0)
all_q_bases = set()
for c in candidates:
    all_q_bases |= c["q_bases"]
all_q_cov_bp = len(all_q_bases)
all_q_cov_ratio = all_q_cov_bp / mito_len if mito_len > 0 else 0.0

# Greedy set cover over true union-of-query bases.
target_bases_n = int(math.ceil(target_cov * mito_len))
covered = set()
selected = []
remaining = candidates[:]
while len(covered) < target_bases_n and remaining:
    best_idx = None
    best_key = None
    for i, c in enumerate(remaining):
        n_new = len(c["q_bases"] - covered)
        key = (n_new, c["subject_len"], c["pident_max"], c["bitscore_max"], -c["evalue_min"])
        if best_key is None or key > best_key:
            best_key = key
            best_idx = i
    if best_idx is None or best_key[0] <= 0:
        break
    best = remaining.pop(best_idx)
    selected.append(best)
    covered |= best["q_bases"]

min_frag_n = len(selected)
min_contigs_n = len(set(c["contig"] for c in selected)) if selected else 0
min_subject_len = sum(c["subject_len"] for c in selected)
min_longest = max((c["subject_len"] for c in selected), default=0)
min_q_cov_bp = len(covered)
min_q_cov_ratio = min_q_cov_bp / mito_len if mito_len > 0 else 0.0
min_target_reached = "yes" if min_q_cov_ratio >= target_cov else "no"
frag_ratio = min_frag_n / all_frag_n if all_frag_n > 0 else "NA"
len_ratio = min_subject_len / all_subject_len if all_subject_len > 0 else "NA"
frag_ratio_s = f"{frag_ratio:.6g}" if frag_ratio != "NA" else "NA"
len_ratio_s = f"{len_ratio:.6g}" if len_ratio != "NA" else "NA"

# Write full/minimal diagnostics for every REF_TYPE. FINAL is populated only when policy allows it.
full_sorted = sorted(candidates, key=lambda x: (x["contig"], x["start0"], x["end0"]))
minimal_sorted = selected[:]
with open(full_tsv, "a") as out_tsv, open(full_bed, "a") as out_bed:
    for c in full_sorted:
        write_candidate(out_tsv, c)
        write_bed(out_bed, c)
with open(minimal_tsv, "a") as out_tsv, open(minimal_bed, "a") as out_bed:
    for c in minimal_sorted:
        write_candidate(out_tsv, c)
        write_bed(out_bed, c)

final_selected = []
if ref_type == "#A":
    if should_final_mask and all_frag_n > 0:
        final_selected = minimal_sorted
        final_strategy = "A_FINAL_minimal_non_chrM_mask" if min_target_reached == "yes" else "A_FINAL_minimal_non_chrM_mask_target_not_reached"
    else:
        final_strategy = "A_diagnostic_non_chrM_numt_candidates_only"
elif ref_type == "#C-likely_incomp":
    final_strategy = "no_mask_likely_incomplete_reference"
elif ref_type == "#C-likely_comp":
    if all_frag_n == 0:
        final_strategy = "C_likely_comp_no_highconf_candidates_check_thresholds"
    elif should_final_mask:
        final_selected = minimal_sorted
        final_strategy = "C_likely_comp_FINAL_minimal_mask" if min_target_reached == "yes" else "C_likely_comp_FINAL_minimal_mask_target_not_reached"
    else:
        final_strategy = "C_likely_comp_diagnostic_only_mask_not_requested"
elif ref_type == "#C-Ambiguous":
    if all_frag_n == 0:
        final_strategy = "C_Ambiguous_no_highconf_candidates_check_thresholds"
    elif should_final_mask:
        final_selected = minimal_sorted
        final_strategy = "C_Ambiguous_FINAL_minimal_mask" if min_target_reached == "yes" else "C_Ambiguous_FINAL_minimal_mask_target_not_reached"
    else:
        final_strategy = "C_Ambiguous_diagnostic_only_mask_not_requested"
else:
    final_strategy = "no_mask_noneligible_reference"

with open(final_tsv, "a") as out_tsv, open(final_bed, "a") as out_bed:
    for c in final_selected:
        write_candidate(out_tsv, c)
        write_bed(out_bed, c)

priority = final_strategy
if ref_type == "#A":
    a_diag = [all_frag_n, all_contigs_n, all_subject_len, all_longest, all_q_cov_bp, f"{all_q_cov_ratio:.6g}",
              min_frag_n, min_contigs_n, min_subject_len, f"{min_q_cov_ratio:.6g}", min_target_reached]
else:
    a_diag = [0,0,0,0,0,"0",0,0,0,"0","no"]

print("	".join(map(str, [
    all_frag_n, all_contigs_n, all_subject_len, all_longest, all_q_cov_bp, f"{all_q_cov_ratio:.6g}",
    min_frag_n, min_contigs_n, min_subject_len, min_longest, min_q_cov_bp, f"{min_q_cov_ratio:.6g}",
    min_target_reached, frag_ratio_s, len_ratio_s, priority, minimal_bed, full_bed, candidate_tsv,
    *a_diag
])))
PY
}

# ---- merge all per-species summary files ----
merge_all_summaries() {
  local merged="${OUTDIR}/all_species.in_house_summary.with_numt_mask.tsv"
  local missing="${OUTDIR}/all_species.missing_summary.txt"
  local masked_min="${OUTDIR}/masked_refs.minimal_chrMcover.fragments.tsv"
  local masked_full="${OUTDIR}/masked_refs.full_highconf.fragments.tsv"
  local masked_min_bed="${OUTDIR}/masked_refs.minimal_chrMcover.bed"
  local masked_full_bed="${OUTDIR}/masked_refs.full_highconf.bed"
  local masked_final="${OUTDIR}/masked_refs.FINAL_numt_mask.fragments.tsv"
  local masked_final_bed="${OUTDIR}/masked_refs.FINAL_numt_mask.bed"
  local a_candidates="${OUTDIR}/A.non_chrM_numt_candidates.fragments.tsv"
  local a_candidates_bed="${OUTDIR}/A.non_chrM_numt_candidates.bed"
  local a_min="${OUTDIR}/A.non_chrM_minimal_chrMcover.fragments.tsv"
  local a_min_bed="${OUTDIR}/A.non_chrM_minimal_chrMcover.bed"
  local lockdir="${OUTDIR}/.merge_all_summaries.lock"
  local tmp_prefix
  tmp_prefix="$(mktemp -d "${OUTDIR}/.merge_tmp.XXXXXX")"
  local expected_cols
  expected_cols=$(awk -F'\t' '{print NF; exit}' <<< "$SUMMARY_HEADER")

  acquire_merge_lock "$lockdir"
  trap 'release_merge_lock "$lockdir"; rm -rf "$tmp_prefix"' RETURN

  log "Merging per-reference summaries into ${merged}"
  printf "%s\n" "$SUMMARY_HEADER" > "${tmp_prefix}/merged"
  : > "${tmp_prefix}/missing"

  local idx sp summary_id summary_file
  for idx in "${!ALL_SPECIES[@]}"; do
    sp="${ALL_SPECIES[$idx]}"
    summary_id="${SUMMARY_IDS[$idx]}"
    summary_file="${OUTDIR}/${summary_id}.summary.tsv"
    if [[ -s "$summary_file" ]]; then
      append_valid_summary_body "$summary_file" "${tmp_prefix}/merged" "$expected_cols"
    else
      printf "%s\t%s\n" "$sp" "$summary_id" >> "${tmp_prefix}/missing"
    fi
  done

  mv -f "${tmp_prefix}/merged" "$merged"
  mv -f "${tmp_prefix}/missing" "$missing"

  # Merge mask-eligible (#C-Ambiguous and #C-likely_comp) fragment lists. Keep one header.
  printf "%s\n" $'Species\tFragmentID\tSubjectContig\tSubjectStart1\tSubjectEnd1\tSubjectStart0\tSubjectEnd0\tSubjectLen\tQueryIntervals\tQueryCoveredBp\tPidentMax\tPidentMean\tBitscoreMax\tEvalueMin\tHSPsN' > "${tmp_prefix}/masked_min"
  printf "%s\n" $'Species\tFragmentID\tSubjectContig\tSubjectStart1\tSubjectEnd1\tSubjectStart0\tSubjectEnd0\tSubjectLen\tQueryIntervals\tQueryCoveredBp\tPidentMax\tPidentMean\tBitscoreMax\tEvalueMin\tHSPsN' > "${tmp_prefix}/masked_full"
  printf "%s\n" $'#chrom\tstart0\tend0\tfragment_id\tspecies\tq_covered_bp\tsubject_len\tpident_max\tbitscore_max' > "${tmp_prefix}/masked_min_bed"
  printf "%s\n" $'#chrom\tstart0\tend0\tfragment_id\tspecies\tq_covered_bp\tsubject_len\tpident_max\tbitscore_max' > "${tmp_prefix}/masked_full_bed"
  printf "%s\n" $'Species\tFragmentID\tSubjectContig\tSubjectStart1\tSubjectEnd1\tSubjectStart0\tSubjectEnd0\tSubjectLen\tQueryIntervals\tQueryCoveredBp\tPidentMax\tPidentMean\tBitscoreMax\tEvalueMin\tHSPsN' > "${tmp_prefix}/masked_final"
  printf "%s\n" $'#chrom\tstart0\tend0\tfragment_id\tspecies\tq_covered_bp\tsubject_len\tpident_max\tbitscore_max' > "${tmp_prefix}/masked_final_bed"

  for idx in "${!ALL_SPECIES[@]}"; do
    sp="${ALL_SPECIES[$idx]}"
    sp_id="${SUMMARY_IDS[$idx]}"
    [[ -s "${OUTDIR}/numt_beds/${sp_id}.minimal_chrMcover.tsv" ]] && tail -n +2 "${OUTDIR}/numt_beds/${sp_id}.minimal_chrMcover.tsv" >> "${tmp_prefix}/masked_min"
    [[ -s "${OUTDIR}/numt_beds/${sp_id}.full_highconf.tsv" ]] && tail -n +2 "${OUTDIR}/numt_beds/${sp_id}.full_highconf.tsv" >> "${tmp_prefix}/masked_full"
    [[ -s "${OUTDIR}/numt_beds/${sp_id}.minimal_chrMcover.bed" ]] && tail -n +2 "${OUTDIR}/numt_beds/${sp_id}.minimal_chrMcover.bed" >> "${tmp_prefix}/masked_min_bed"
    [[ -s "${OUTDIR}/numt_beds/${sp_id}.full_highconf.bed" ]] && tail -n +2 "${OUTDIR}/numt_beds/${sp_id}.full_highconf.bed" >> "${tmp_prefix}/masked_full_bed"
    [[ -s "${OUTDIR}/numt_beds/${sp_id}.FINAL_numt_mask.tsv" ]] && tail -n +2 "${OUTDIR}/numt_beds/${sp_id}.FINAL_numt_mask.tsv" >> "${tmp_prefix}/masked_final"
    [[ -s "${OUTDIR}/numt_beds/${sp_id}.FINAL_numt_mask.bed" ]] && tail -n +2 "${OUTDIR}/numt_beds/${sp_id}.FINAL_numt_mask.bed" >> "${tmp_prefix}/masked_final_bed"
    if [[ -s "${OUTDIR}/${sp_id}.summary.tsv" ]] && awk -F'	' 'NR==2 && $2=="#A" {found=1} END{exit !found}' "${OUTDIR}/${sp_id}.summary.tsv"; then
      [[ -s "${OUTDIR}/numt_candidates/${sp_id}.highconf_reference_mtlike_candidates.tsv" ]] && tail -n +2 "${OUTDIR}/numt_candidates/${sp_id}.highconf_reference_mtlike_candidates.tsv" >> "${tmp_prefix}/a_candidates"
      [[ -s "${OUTDIR}/numt_beds/${sp_id}.full_highconf.bed" ]] && tail -n +2 "${OUTDIR}/numt_beds/${sp_id}.full_highconf.bed" >> "${tmp_prefix}/a_candidates_bed"
      [[ -s "${OUTDIR}/numt_beds/${sp_id}.minimal_chrMcover.tsv" ]] && tail -n +2 "${OUTDIR}/numt_beds/${sp_id}.minimal_chrMcover.tsv" >> "${tmp_prefix}/a_min"
      [[ -s "${OUTDIR}/numt_beds/${sp_id}.minimal_chrMcover.bed" ]] && tail -n +2 "${OUTDIR}/numt_beds/${sp_id}.minimal_chrMcover.bed" >> "${tmp_prefix}/a_min_bed"
    fi
  done


  mv -f "${tmp_prefix}/masked_min" "$masked_min"
  mv -f "${tmp_prefix}/masked_full" "$masked_full"
  mv -f "${tmp_prefix}/masked_min_bed" "$masked_min_bed"
  mv -f "${tmp_prefix}/masked_full_bed" "$masked_full_bed"
  mv -f "${tmp_prefix}/masked_final" "$masked_final"
  mv -f "${tmp_prefix}/masked_final_bed" "$masked_final_bed"
  mv -f "${tmp_prefix}/a_candidates" "$a_candidates"
  mv -f "${tmp_prefix}/a_candidates_bed" "$a_candidates_bed"
  mv -f "${tmp_prefix}/a_min" "$a_min"
  mv -f "${tmp_prefix}/a_min_bed" "$a_min_bed"

  local n_total n_merged n_missing
  n_total="${#ALL_SPECIES[@]}"
  n_merged=$(($(wc -l < "$merged") - 1))
  n_missing=$(wc -l < "$missing")
  log "Merge complete: ${n_merged}/${n_total} summaries merged; ${n_missing} missing."
  log "Merged mask-eligible minimal fragments: ${masked_min}"
  log "Merged mask-eligible full fragments: ${masked_full}"
  log "Merged mask-eligible FINAL mask fragments: ${masked_final}"
  log "Merged #A non-chrM diagnostic candidates: ${a_candidates}"
  log "REF_TYPE vs MaskPriority summary:"
  awk -F'	' 'NR==1{for(i=1;i<=NF;i++){if($i=="REF_TYPE") r=i; if($i=="MaskPriority") m=i} next} {k=$r"	"$m; c[k]++} END{for(k in c) print "  " k "	" c[k]}' "$merged" >&2
  awk -F'	' 'NR==1{for(i=1;i<=NF;i++){if($i=="REF_TYPE") r=i; if($i=="MaskPriority") m=i} next} $r=="#C-likely_comp" || $r=="#C-Ambiguous" {eligible++} $m ~ /FINAL_.*mask/ {applied++} END{print "[INFO] final mask applied count: " applied+0 "; default eligible rows: " eligible+0 > "/dev/stderr"}' "$merged"
  if [[ "$n_missing" -gt 0 ]]; then
    warn "Missing summary list: ${missing}"
  fi
}

process_species() {

  local index="$1"
  species="${SPECIES_NAMES[$index]}"
  GENOME="${WG_FASTA_PATHS[$index]}"
  MITO_REF="${CHRM_FASTA_PATHS[$index]}"
  SPECIES_ID="${SUMMARY_IDS[$index]}"
  BLAST_DB="${OUTDIR}/${SPECIES_ID}_db"
  STATUS_FILE="${OUTDIR}/${SPECIES_ID}.status"
  SUMMARY_FILE="${OUTDIR}/${SPECIES_ID}.summary.tsv"

  echo "Running" > "$STATUS_FILE"
  log "Processing ${species}"

  if [[ ! -f "$GENOME" || ! -f "$MITO_REF" ]]; then
    warn "Missing input files for ${species}"
    printf "%s\n" "$SUMMARY_HEADER" > "$SUMMARY_FILE"
    printf "%s\tMissingFile\tmissing_input\tNA\t0\t0\tNA\t0\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tno\tNA\tNA\tmissing_input\tNA\tNA\tNA\n" "$species" >> "$SUMMARY_FILE"
    echo "Failed" > "$STATUS_FILE"
    return 0
  fi

  TEMP_MITO="${OUTDIR}/${SPECIES_ID}_temp_mito.fa"
  BLAST_TEMP="${OUTDIR}/${SPECIES_ID}_blast_all.tsv"
  TOP_BLAST_OUT="${OUTDIR}/${SPECIES_ID}_blast_top.tsv"

  fasta_header_full "$species" "chrM" > "$TEMP_MITO"
  grep -v '^>' "$MITO_REF" | format_sequence >> "$TEMP_MITO"

  MITO_LEN=$(awk '/^>/ {if (seqlen){print seqlen}; seqlen=0; next} {gsub(/[[:space:]]/, ""); seqlen += length($0)} END {print seqlen+0}' "$TEMP_MITO")
  [[ -z "${MITO_LEN:-}" ]] && MITO_LEN=0

  read -r ANNOTATED_CONTIG ANNOTATED_CONTIG_LEN HAS_VALID_ANNOTATED_MITO < <(
    select_annotated_mito_contig_py "$GENOME" "$MITO_LEN"
  )

  if ! "$MAKEBLASTDB_COMMAND" -in "$GENOME" -dbtype nucl -out "$BLAST_DB" &> "${IN_HOUSE_SCORE_LOG_DIR}/${SPECIES_ID}_makeblastdb.log"; then
    err "makeblastdb failed for ${species}"
    printf "%s\n" "$SUMMARY_HEADER" > "$SUMMARY_FILE"
    printf "%s\tMakeDBFail\tmakeblastdb_failed\t%s\t%s\t%s\tNA\t0\t%s\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tno\tNA\tNA\tmakeblastdb_failed\tNA\tNA\tNA\n" \
      "$species" "$ANNOTATED_CONTIG" "$ANNOTATED_CONTIG_LEN" "$HAS_VALID_ANNOTATED_MITO" "$MITO_LEN" >> "$SUMMARY_FILE"
    echo "Failed" > "$STATUS_FILE"
    return 0
  fi

  "$BLASTN_COMMAND" -query "$TEMP_MITO" \
    -db "$BLAST_DB" \
    -out "$BLAST_TEMP" \
    -outfmt "6 qseqid sseqid qstart qend sstart send length pident mismatch evalue bitscore" \
    -evalue 1e-3 \
    -perc_identity 90 \
    -num_threads "$THREADS" \
    -max_target_seqs 10000 \
    -task blastn &> "${IN_HOUSE_SCORE_LOG_DIR}/${SPECIES_ID}_blastn.log" || warn "blastn had non-zero exit for ${species}"

  SCORE="0"; M_CLIP="0"; DCPLUS="0"; DCMINUS="0"; DTPLUS="0"; DTMINUS="0"
  REF_TYPE="#C-likely_incomp"; MTLIKE_PATTERN="weak_or_no_mtlike_signal"
  TOP_CONTIG="NoHit"; TOP_CONTIG_LEN=0; MAX_HIT_RATIO=0; TOTAL_HIT_RATIO=0; CONTIG_RATIO=0
  TOP_MERGED_LEN=0; TOP_MERGED_RATIO=0; LONGEST_MERGED_HIT=0; CUM_MERGED_LEN=0; MERGED_CONTIGS_N=0; MERGED_INTERVALS_N=0
  NONTOP_LONGEST_MERGED=0; NONTOP_CUM_MERGED=0; NONTOP_CONTIGS_N=0; NONTOP_INTERVALS_N=0

  if [[ ! -s "$BLAST_TEMP" ]]; then
    warn "No BLAST hits for ${species}"
    : > "$TOP_BLAST_OUT"
    read -r SCORE M_CLIP DCPLUS DCMINUS DTPLUS DTMINUS < <(calc_score_py "$MAX_HIT_RATIO" "$TOTAL_HIT_RATIO" "$CONTIG_RATIO")
  else
    read -r TOP_CONTIG TOP_MERGED_LEN TOP_MERGED_RATIO LONGEST_MERGED_HIT CUM_MERGED_LEN MERGED_CONTIGS_N MERGED_INTERVALS_N NONTOP_LONGEST_MERGED NONTOP_CUM_MERGED NONTOP_CONTIGS_N NONTOP_INTERVALS_N < <(
      merge_blast_stats_py "$BLAST_TEMP" "$ANNOTATED_CONTIG" "$MITO_LEN"
    )

    awk -v c="$TOP_CONTIG" '$2==c' "$BLAST_TEMP" > "$TOP_BLAST_OUT"
    MAX_HIT_LENGTH=$(awk '{print $7}' "$TOP_BLAST_OUT" | sort -nr | head -1)
    [[ -z "${MAX_HIT_LENGTH:-}" ]] && MAX_HIT_LENGTH=0

    TOP_CONTIG_LEN=$(awk -v c="$TOP_CONTIG" '
      BEGIN{keep=0; len=0}
      /^>/ {header=$0; sub(/^>/,"",header); split(header,a," "); keep=(a[1]==c); next}
      keep {gsub(/[[:space:]]/, "", $0); len += length($0)}
      END{print len+0}
    ' "$GENOME")
    [[ -z "${TOP_CONTIG_LEN:-}" ]] && TOP_CONTIG_LEN=0

    TOTAL_HIT_LENGTH=$(awk '{sum+=$7} END {print sum+0}' "$TOP_BLAST_OUT")
    [[ -z "${TOTAL_HIT_LENGTH:-}" ]] && TOTAL_HIT_LENGTH=0

    if [[ "$MITO_LEN" -le 0 ]]; then
      MAX_HIT_RATIO=0; TOTAL_HIT_RATIO=0; CONTIG_RATIO=0
    else
      read -r MAX_HIT_RATIO TOTAL_HIT_RATIO CONTIG_RATIO < <(
        "$PYTHON_COMMAND" - "$MAX_HIT_LENGTH" "$TOTAL_HIT_LENGTH" "$TOP_CONTIG_LEN" "$MITO_LEN" <<'PY'
import sys
max_hit, total_hit, top_len, mito_len = map(float, sys.argv[1:])
if mito_len <= 0:
    print("0 0 0")
else:
    print(f"{max_hit/mito_len:.6g} {total_hit/mito_len:.6g} {top_len/mito_len:.6g}")
PY
      )
    fi

    read -r SCORE M_CLIP DCPLUS DCMINUS DTPLUS DTMINUS < <(calc_score_py "$MAX_HIT_RATIO" "$TOTAL_HIT_RATIO" "$CONTIG_RATIO")

    read -r REF_TYPE MTLIKE_PATTERN < <(classify_ref_py \
      "$HAS_VALID_ANNOTATED_MITO" "$TOP_CONTIG_LEN" "$MITO_LEN" \
      "$TOP_MERGED_LEN" "$TOP_MERGED_RATIO" "$LONGEST_MERGED_HIT" \
      "$CUM_MERGED_LEN" "$MERGED_CONTIGS_N" "$MERGED_INTERVALS_N" \
      "$SCORE" "$MAX_HIT_RATIO" "$TOTAL_HIT_RATIO" "$CONTIG_RATIO")
  fi

  ANNOTATED_CONTIG_FOR_NUMT="$ANNOTATED_CONTIG"
  read -r ALL_FRAG_N ALL_CONTIG_N ALL_SUBJECT_LEN ALL_LONGEST ALL_Q_BP ALL_Q_RATIO MIN_FRAG_N MIN_CONTIG_N MIN_SUBJECT_LEN MIN_LONGEST MIN_Q_BP MIN_Q_RATIO MIN_TARGET FRAG_RATIO LEN_RATIO MASK_PRIORITY MIN_BED FULL_BED CAND_TSV A_NONCHRM_FRAG_N A_NONCHRM_CONTIG_N A_NONCHRM_SUBJECT_LEN A_NONCHRM_LONGEST A_NONCHRM_Q_BP A_NONCHRM_Q_RATIO A_NONCHRM_MIN_FRAG_N A_NONCHRM_MIN_CONTIG_N A_NONCHRM_MIN_SUBJECT_LEN A_NONCHRM_MIN_Q_RATIO A_NONCHRM_MIN_TARGET < <(
    numt_mask_summary_py "$species" "$BLAST_TEMP" "$MITO_LEN" "$REF_TYPE" "$SPECIES_ID"
  )

  {
    printf "%s\n" "$SUMMARY_HEADER"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$species" "$REF_TYPE" "$MTLIKE_PATTERN" "$ANNOTATED_CONTIG" "$ANNOTATED_CONTIG_LEN" "$HAS_VALID_ANNOTATED_MITO" \
      "${TOP_CONTIG:-NA}" "${TOP_CONTIG_LEN:-0}" "$MITO_LEN" \
      "${MAX_HIT_RATIO:-0}" "${TOTAL_HIT_RATIO:-0}" "${CONTIG_RATIO:-0}" \
      "${SCORE:-0}" "${M_CLIP:-0}" "${DCPLUS:-0}" "${DCMINUS:-0}" "${DTPLUS:-0}" "${DTMINUS:-0}" \
      "${TOP_MERGED_LEN:-0}" "${TOP_MERGED_RATIO:-0}" "${LONGEST_MERGED_HIT:-0}" "${CUM_MERGED_LEN:-0}" \
      "${MERGED_CONTIGS_N:-0}" "${MERGED_INTERVALS_N:-0}" \
      "${NONTOP_LONGEST_MERGED:-0}" "${NONTOP_CUM_MERGED:-0}" "${NONTOP_CONTIGS_N:-0}" "${NONTOP_INTERVALS_N:-0}" \
      "$ALL_FRAG_N" "$ALL_CONTIG_N" "$ALL_SUBJECT_LEN" "$ALL_LONGEST" "$ALL_Q_BP" "$ALL_Q_RATIO" \
      "$MIN_FRAG_N" "$MIN_CONTIG_N" "$MIN_SUBJECT_LEN" "$MIN_LONGEST" "$MIN_Q_BP" "$MIN_Q_RATIO" "$MIN_TARGET" \
      "$FRAG_RATIO" "$LEN_RATIO" "$MASK_PRIORITY" "$MIN_BED" "$FULL_BED" "$CAND_TSV" \
      "$A_NONCHRM_FRAG_N" "$A_NONCHRM_CONTIG_N" "$A_NONCHRM_SUBJECT_LEN" "$A_NONCHRM_LONGEST" "$A_NONCHRM_Q_BP" "$A_NONCHRM_Q_RATIO" \
      "$A_NONCHRM_MIN_FRAG_N" "$A_NONCHRM_MIN_CONTIG_N" "$A_NONCHRM_MIN_SUBJECT_LEN" "$A_NONCHRM_MIN_Q_RATIO" "$A_NONCHRM_MIN_TARGET"
  } > "$SUMMARY_FILE"

  rm -f "$TEMP_MITO" "$TOP_BLAST_OUT" 2>/dev/null || true
  rm -f "${BLAST_DB}".n* 2>/dev/null || true
  # Keep BLAST_TEMP because it is useful for debugging. Uncomment below if you prefer cleanup.
  # rm -f "$BLAST_TEMP" 2>/dev/null || true

  echo "Done" > "$STATUS_FILE"
  log "Successfully processed ${species}"
}

# -------------------- Execution --------------------
if [[ ! -s "$REF_INPUTS" ]]; then
  err "Missing or empty in-house score reference input manifest: ${REF_INPUTS}"
  exit 1
fi

mapfile -t REFERENCE_ROWS < <(
  "$PYTHON_COMMAND" - "$REF_INPUTS" <<'PY'
import csv, sys
manifest = sys.argv[1]
seen = set()
with open(manifest, newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        species = (row.get("target_species") or "").strip()
        wg = (row.get("wg_fasta_path") or "").strip()
        chrm = (row.get("chrM_fasta_path") or "").strip()
        if not species or not wg or not chrm:
            continue
        key = (species, wg, chrm)
        if key in seen:
            continue
        seen.add(key)
        ordinal = len(seen)
        safe = ''.join(ch if ch.isalnum() or ch in '_.-' else '_' for ch in species.replace(' ', '_').replace('/', '_'))
        print("\t".join([species, wg, chrm, f"{safe}_ref{ordinal:04d}"]))
PY
)

if [[ "${#REFERENCE_ROWS[@]}" -eq 0 ]]; then
  err "No usable reference-level rows found in ${REF_INPUTS}; required columns are target_species, wg_fasta_path, chrM_fasta_path."
  exit 1
fi

ALL_SPECIES=()
SPECIES_NAMES=()
WG_FASTA_PATHS=()
CHRM_FASTA_PATHS=()
SUMMARY_IDS=()
for row in "${REFERENCE_ROWS[@]}"; do
  IFS=$'\t' read -r species wg_fasta chrM_fasta summary_id <<< "$row"
  ALL_SPECIES+=("$species")
  SPECIES_NAMES+=("$species")
  WG_FASTA_PATHS+=("$wg_fasta")
  CHRM_FASTA_PATHS+=("$chrM_fasta")
  SUMMARY_IDS+=("$summary_id")
done

N_SPECIES="${#ALL_SPECIES[@]}"

if [[ "${MERGE_ONLY:-0}" == "1" ]]; then
  merge_all_summaries
  tmp_merged_score="$(mktemp "${OUTDIR}/.merged_in_house_score.XXXXXX")"
  cp "${OUTDIR}/all_species.in_house_summary.with_numt_mask.tsv" "$tmp_merged_score"
  mv -f "$tmp_merged_score" "$MERGED_IN_HOUSE_SCORE"
  exit 0
fi

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  if [[ "${SUBMIT_ARRAY:-0}" == "1" ]]; then
    log "Submitting ${N_SPECIES} in-house score array tasks with concurrency ${ARRAY_CONCURRENCY}."
    jid=$(sbatch --array="1-${N_SPECIES}%${ARRAY_CONCURRENCY}" "$0" | awk '{print $4}')
    log "Submitted array job ${jid}; submitting dependent MERGE_ONLY job."
    sbatch --dependency="afterok:${jid}" --export=ALL,MERGE_ONLY=1 "$0"
    exit 0
  fi
  err "SLURM_ARRAY_TASK_ID is not set. Submit array with: sbatch --array=1-${N_SPECIES}%${ARRAY_CONCURRENCY} $0, use SUBMIT_ARRAY=1 bash $0, or merge with: MERGE_ONLY=1 bash $0"
  exit 1
fi

if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N_SPECIES )); then
  err "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} is outside species list range 1-${N_SPECIES}"
  exit 1
fi

SPECIES_INDEX=$((SLURM_ARRAY_TASK_ID-1))
SPECIES="${ALL_SPECIES[$SPECIES_INDEX]}"
log "Selected species index ${SLURM_ARRAY_TASK_ID}/${N_SPECIES}: ${SPECIES}"

{
  process_species "$SPECIES_INDEX"
} || {
  err "Processing failed for $SPECIES"
  echo "Failed" > "${OUTDIR}/${SUMMARY_IDS[$SPECIES_INDEX]}.status"
  exit 1
}

exit 0
