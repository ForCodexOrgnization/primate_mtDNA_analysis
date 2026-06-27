#!/usr/bin/env python3
"""Check integrity of variant-calling reference packages.

This script verifies NUMT-mask application and chrM append/rename behavior for
successful rows in the variant-calling reference manifest.

Example:
    python preprocessing/scripts/check_variant_reference_integrity.py \
      --repo-root . \
      --manifest references/variant_calling/variant_calling_reference_manifest.tsv \
      --mask-ref-types "#C-likely_comp,#C-Ambiguous,#A" \
      --out results/preprocessing/reports/variant_reference_integrity_check.tsv

If only manifest shards exist, first merge the shards into a temporary manifest
and pass that merged TSV to --manifest.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

FINAL_MASK_PRIORITIES = {
    "C_likely_comp_FINAL_minimal_mask",
    "C_likely_comp_FINAL_minimal_mask_target_not_reached",
    "C_Ambiguous_FINAL_minimal_mask",
    "C_Ambiguous_FINAL_minimal_mask_target_not_reached",
    "A_FINAL_minimal_non_chrM_mask",
    "A_FINAL_minimal_non_chrM_mask_target_not_reached",
}

MISSING_VALUES = {"", "NA", "N/A", "NULL", "None", "null", "nan"}
PREFERRED_COLUMNS = [
    "target_species",
    "safe_species_id",
    "REF_TYPE",
    "MaskPriority",
    "expected_chrM_mode",
    "whole_chrM_count",
    "whole_chrM_equals_Ref_chrM",
    "whole_chrM_len",
    "Ref_chrM_len",
    "manifest_mask_applied",
    "expected_mask_applied",
    "bed_exists",
    "bed_n_intervals",
    "bed_total_bp",
    "mask_interval_bp_present_in_whole",
    "mask_interval_N_count",
    "mask_interval_nonN_count",
    "missing_bed_contig_count",
    "whole_non_chrM_total_N",
    "check_status",
    "check_message",
]
EXTRA_COLUMNS = [
    "ValidAnnotatedMitoContig",
    "HasValidAnnotatedMito",
    "numt_mask_bed",
    "whole_fasta",
    "chrM_fasta",
    "whole_chrM_md5",
    "Ref_chrM_md5",
    "whole_chrM_N_count",
    "Ref_chrM_N_count",
    "valid_contig_still_present_in_whole",
    "missing_bed_contigs_first10",
]
OUTPUT_COLUMNS = PREFERRED_COLUMNS + EXTRA_COLUMNS


def is_missing(value: object) -> bool:
    return str(value or "").strip() in MISSING_VALUES


def yes_no(value: bool) -> str:
    return "yes" if value else "no"


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [{key: (value or "") for key, value in row.items()} for row in reader]


def fasta_iter(path: Path) -> Iterable[Tuple[str, str]]:
    name: Optional[str] = None
    seq: List[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                header = line[1:].strip()
                name = header.split()[0] if header else ""
                seq = []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def read_fasta_records(path: Path) -> List[Tuple[str, str]]:
    return list(fasta_iter(path))


def read_fasta_dict(path: Path) -> Dict[str, str]:
    records: Dict[str, str] = {}
    for name, seq in fasta_iter(path):
        records[name] = seq
    return records


def md5_seq(seq: Optional[str]) -> str:
    if seq is None:
        return ""
    return hashlib.md5(seq.encode("utf-8")).hexdigest()


def read_bed(path: Optional[Path]) -> List[Tuple[str, int, int]]:
    if path is None or not path.exists():
        return []
    intervals: List[Tuple[str, int, int]] = []
    with path.open() as handle:
        for raw in handle:
            if not raw.strip() or raw.lstrip().startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 3:
                fields = raw.split()
            if len(fields) < 3:
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                continue
            if end > start:
                intervals.append((fields[0], start, end))
    return intervals


def resolve(repo_root: Path, value: object) -> Optional[Path]:
    if is_missing(value):
        return None
    path = Path(str(value).strip())
    return path if path.is_absolute() else repo_root / path


def expected_mask(row: Dict[str, str], mask_ref_types: Sequence[str], bed: Optional[Path], intervals: Sequence[Tuple[str, int, int]]) -> bool:
    return (
        row.get("REF_TYPE", "") in set(mask_ref_types)
        and row.get("MaskPriority", "") in FINAL_MASK_PRIORITIES
        and bed is not None
        and bed.exists()
        and len(intervals) > 0
    )


def count_n_in_interval(seq: str, start: int, end: int) -> int:
    segment = seq[max(0, start):min(len(seq), end)]
    return sum(1 for base in segment if base.upper() == "N")


def count_non_n_in_interval(seq: str, start: int, end: int) -> int:
    segment = seq[max(0, start):min(len(seq), end)]
    return sum(1 for base in segment if base.upper() != "N")


def expected_chrm_mode(row: Dict[str, str]) -> str:
    if row.get("HasValidAnnotatedMito", "").strip() == "1" and not is_missing(row.get("ValidAnnotatedMitoContig", "")):
        return "rename_embedded_mito_contig_to_chrM"
    return "append_chrM_to_whole_ref"


def base_output_row(row: Dict[str, str]) -> Dict[str, str]:
    return {col: "" for col in OUTPUT_COLUMNS} | {
        "target_species": row.get("target_species", ""),
        "safe_species_id": row.get("safe_species_id", ""),
        "REF_TYPE": row.get("REF_TYPE", ""),
        "MaskPriority": row.get("MaskPriority", ""),
        "ValidAnnotatedMitoContig": row.get("ValidAnnotatedMitoContig", ""),
        "HasValidAnnotatedMito": row.get("HasValidAnnotatedMito", ""),
        "numt_mask_bed": row.get("numt_mask_bed", ""),
        "whole_fasta": row.get("whole_fasta", ""),
        "chrM_fasta": row.get("chrM_fasta", ""),
        "expected_chrM_mode": expected_chrm_mode(row),
        "manifest_mask_applied": yes_no(row.get("numt_mask_applied_to_whole_ref", "").strip().lower() == "yes"),
    }


def check_row(row: Dict[str, str], repo_root: Path, mask_ref_types: Sequence[str]) -> Dict[str, str]:
    out = base_output_row(row)
    messages: List[str] = []
    try:
        whole_path = resolve(repo_root, row.get("whole_fasta"))
        chrm_path = resolve(repo_root, row.get("chrM_fasta"))
        bed_path = resolve(repo_root, row.get("numt_mask_bed"))
        bed_intervals = read_bed(bed_path)
        out["bed_exists"] = yes_no(bed_path is not None and bed_path.exists())
        out["bed_n_intervals"] = str(len(bed_intervals))
        out["bed_total_bp"] = str(sum(end - start for _, start, end in bed_intervals))
        exp_mask = expected_mask(row, mask_ref_types, bed_path, bed_intervals)
        out["expected_mask_applied"] = yes_no(exp_mask)
        manifest_mask = row.get("numt_mask_applied_to_whole_ref", "").strip().lower() == "yes"
        if exp_mask and not manifest_mask:
            messages.append("expected_mask_yes_but_manifest_no")
        if not exp_mask and manifest_mask:
            messages.append("manifest_mask_yes_but_expected_no")
        if whole_path is None or not whole_path.exists():
            messages.append("missing_whole_fasta")
            return finish(out, messages)
        if chrm_path is None or not chrm_path.exists():
            messages.append("missing_chrM_fasta")
            return finish(out, messages)

        whole_records = read_fasta_records(whole_path)
        chrm_records = read_fasta_records(chrm_path)
        whole = {name: seq for name, seq in whole_records}
        chrm = {name: seq for name, seq in chrm_records}
        whole_chrM_records = [(name, seq) for name, seq in whole_records if name == "chrM"]
        ref_chrM_records = [(name, seq) for name, seq in chrm_records if name == "chrM"]
        whole_chrm_seq = whole_chrM_records[0][1] if len(whole_chrM_records) == 1 else whole.get("chrM")
        ref_chrm_seq = ref_chrM_records[0][1] if len(ref_chrM_records) == 1 else chrm.get("chrM")
        whole_chrm_count = len(whole_chrM_records)
        out["whole_chrM_count"] = str(whole_chrm_count)
        out["whole_chrM_len"] = str(len(whole_chrm_seq or ""))
        out["Ref_chrM_len"] = str(len(ref_chrm_seq or ""))
        out["whole_chrM_equals_Ref_chrM"] = yes_no(whole_chrm_seq is not None and whole_chrm_seq == ref_chrm_seq)
        out["whole_chrM_md5"] = md5_seq(whole_chrm_seq)
        out["Ref_chrM_md5"] = md5_seq(ref_chrm_seq)
        out["whole_chrM_N_count"] = str(count_n_in_interval(whole_chrm_seq or "", 0, len(whole_chrm_seq or "")))
        out["Ref_chrM_N_count"] = str(count_n_in_interval(ref_chrm_seq or "", 0, len(ref_chrm_seq or "")))
        out["whole_non_chrM_total_N"] = str(sum(count_n_in_interval(seq, 0, len(seq)) for name, seq in whole.items() if name != "chrM"))
        if whole_chrm_count != 1:
            messages.append(f"whole_chrM_count_not_one:{whole_chrm_count}")
        if [name for name, _seq in chrm_records] != ["chrM"]:
            messages.append("Ref_chrM_record_not_exactly_chrM")
        if whole_chrm_seq != ref_chrm_seq:
            messages.append("whole_chrM_sequence_differs_from_Ref_chrM")
        if out["whole_chrM_N_count"] != out["Ref_chrM_N_count"]:
            messages.append("chrM_N_count_changed_between_whole_and_Ref_chrM")
        valid_contig = row.get("ValidAnnotatedMitoContig", "").strip()
        still_present = out["expected_chrM_mode"] == "rename_embedded_mito_contig_to_chrM" and valid_contig != "chrM" and valid_contig in whole
        out["valid_contig_still_present_in_whole"] = yes_no(still_present)
        if still_present:
            messages.append("valid_mito_contig_still_present_after_chrM_rename")

        missing = []
        bp_present = n_count = non_n_count = 0
        if any(contig == "chrM" for contig, _, _ in bed_intervals):
            messages.append("mask_BED_contains_chrM_contig")
        for contig, start, end in bed_intervals:
            if contig == "chrM":
                continue
            seq = whole.get(contig)
            if seq is None:
                missing.append(contig)
                continue
            clipped_start, clipped_end = max(0, start), min(len(seq), end)
            if clipped_end <= clipped_start:
                continue
            bp_present += clipped_end - clipped_start
            n_count += count_n_in_interval(seq, clipped_start, clipped_end)
            non_n_count += count_non_n_in_interval(seq, clipped_start, clipped_end)
        missing_unique = sorted(set(missing))
        out["mask_interval_bp_present_in_whole"] = str(bp_present)
        out["mask_interval_N_count"] = str(n_count)
        out["mask_interval_nonN_count"] = str(non_n_count)
        out["missing_bed_contig_count"] = str(len(missing_unique))
        out["missing_bed_contigs_first10"] = ",".join(missing_unique[:10])
        if exp_mask and non_n_count > 0:
            messages.append(f"BED_intervals_not_fully_masked_nonN_bp:{non_n_count}")
        if exp_mask and n_count == 0:
            messages.append("expected_mask_but_no_N_inside_BED_intervals")
        if exp_mask and bp_present == 0:
            messages.append("expected_mask_but_no_BED_intervals_found_in_whole")
    except Exception as exc:  # keep checking subsequent manifest rows
        messages.append(f"row_check_error:{type(exc).__name__}:{exc}")
    return finish(out, messages)


def finish(out: Dict[str, str], messages: Sequence[str]) -> Dict[str, str]:
    out["check_status"] = "fail" if messages else "pass"
    out["check_message"] = ";".join(messages)
    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--manifest", default="references/variant_calling/variant_calling_reference_manifest.tsv")
    parser.add_argument("--mask-ref-types", default="#C-likely_comp,#C-Ambiguous,#A")
    parser.add_argument("--out", default="results/preprocessing/reports/variant_reference_integrity_check.tsv")
    parser.add_argument("--max-rows", type=int, default=0, help="Maximum success rows to check; 0 means all success rows.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    manifest = resolve(repo_root, args.manifest)
    if manifest is None or not manifest.exists():
        sys.exit(f"ERROR: manifest not found: {manifest or args.manifest}")
    rows = [row for row in read_tsv(manifest) if row.get("build_status", "") == "success"]
    if args.max_rows and args.max_rows > 0:
        rows = rows[:args.max_rows]
    mask_ref_types = [item.strip() for item in args.mask_ref_types.split(",") if item.strip()]
    results = [check_row(row, repo_root, mask_ref_types) for row in rows]
    out_path = resolve(repo_root, args.out)
    if out_path is None:
        sys.exit("ERROR: --out resolved to an empty path")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(results)

    counts = Counter(row["check_status"] for row in results)
    print(f"wrote: {out_path}")
    print(f"pass: {counts.get('pass', 0)}")
    print(f"fail: {counts.get('fail', 0)}")
    failed = [row for row in results if row["check_status"] == "fail"][:20]
    if failed:
        print("first_failed_rows:")
        for row in failed:
            print("\t".join([row.get("target_species", ""), row.get("REF_TYPE", ""), row.get("MaskPriority", ""), row.get("check_message", "")]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
