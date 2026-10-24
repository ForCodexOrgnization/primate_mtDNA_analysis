#!/usr/bin/env python3
"""Collect original full-depth VCFs for the 10-sample HOM downsampling comparison.

Safe matching: identify sample from basename, never select arbitrarily between
multiple VCFs. Copies to unique baseline/ names; writes a source inventory.
"""
import argparse
import csv
import gzip
import os
import shutil
from pathlib import Path

DEFAULT_BASE = Path("/vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis_coverage_only/results/qc_coverage_only/downsampling_pilot")
DEFAULT_SOURCE = Path("/vast/palmer/scratch/lake_nicole/lt692/primate_final_analysis_combine/vcf")
EXTS = (".vcf.gz", ".vcf.bgz", ".vcf")
def is_vcf(path):
    return path.is_file() and any(path.name.endswith(e) for e in EXTS)
def sample_match(name, sample):
    return name.startswith(sample + ".") or name.startswith(sample + "_") or name.startswith(sample + "-")
def copy_checked(src, dst, dry_run):
    if dry_run:
        return "would_copy"
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        if dst.stat().st_size == src.stat().st_size:
            return "exists_same_size"
        raise RuntimeError(f"Refusing to overwrite different-sized file: {dst}")
    temp = dst.with_name("." + dst.name + ".tmp." + str(os.getpid()))
    try:
        shutil.copy2(src, temp)
        os.replace(temp, dst)
    finally:
        temp.unlink(missing_ok=True)
    return "copied"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", type=Path, default=DEFAULT_BASE)
    ap.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    with (args.base / "candidate_downsampling_jobs.tsv").open(newline="") as f:
        entries = list(csv.DictReader(f, delimiter="\t"))
    samples = sorted({r["sample"] for r in entries})
    out = args.output or args.base / "collected_hom_results" / "baseline_vcf"
    if not args.source.is_dir():
        raise SystemExit(f"ERROR: missing source directory: {args.source}")
    all_vcfs = [p for p in args.source.rglob("*") if is_vcf(p)]
    report = []
    for sample in samples:
        candidates = sorted(p for p in all_vcfs if sample_match(p.name, sample))
        if len(candidates) == 0:
            status, src, dst = "MISSING", "", ""
        elif len(candidates) > 1:
            status, src, dst = "AMBIGUOUS", ";".join(map(str,candidates)), ""
        else:
            source = candidates[0]
            suffix = next(e for e in EXTS if source.name.endswith(e))
            target = out / f"{sample}.original_full_depth{suffix}"
            status = copy_checked(source, target, args.dry_run)
            src, dst = str(source), str(target)
        report.append(dict(sample=sample,status=status,source=src,collected_path=dst,
                           candidate_count=len(candidates)))
        print(f"{sample}: {status} ({len(candidates)} matching VCFs)")
        if len(candidates) > 1:
            for candidate in candidates:
                print(f"  CANDIDATE {candidate}")
    if not args.dry_run:
        out.mkdir(parents=True, exist_ok=True)
        with (out / "baseline_manifest.tsv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(report[0]),delimiter="\t")
            writer.writeheader()
            writer.writerows(report)
    print(f"Baseline directory: {out}")
    print("Review ambiguous or missing rows before HOM recall analysis.")

if __name__ == "__main__":
    main()
