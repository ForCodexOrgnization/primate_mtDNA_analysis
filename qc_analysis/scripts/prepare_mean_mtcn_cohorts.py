#!/usr/bin/env python3
"""Prepare two isolated QC result trees while reusing immutable inputs.

The original qc_coverage_only and previous DP50 outputs are never overwritten.
Only qc-specific reports are generated in the new cohort directories.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ORIGINAL = ROOT / "results/qc_coverage_only"
DP50 = ROOT / "results/qc_coverage_only_dp50"
COHORTS = {
    "hom10": {
        "root": ROOT / "results/qc_cohort_hom10_mean_mtcn",
        "config": ROOT / "config/qc_cohort_hom10_mean_mtcn.yaml",
    },
    "homhet40": {
        "root": ROOT / "results/qc_cohort_homhet40_mean_mtcn",
        "config": ROOT / "config/qc_cohort_homhet40_mean_mtcn.yaml",
    },
}
REBUILD = {
    "sample_variant_filtering", "local_heteroplasmy_qc",
    "final_filter_hom_only", "homoplasmy_qc_only",
}
RRNA_FOLDER = "vcf_rrna"


def symlink(src: Path, dst: Path):
    if not src.exists():
        raise RuntimeError("Missing required input: %s" % src)
    if dst.is_symlink():
        if dst.resolve() != src.resolve():
            raise RuntimeError("Conflicting link: %s" % dst)
        return
    if dst.exists():
        raise RuntimeError("Refusing to replace existing path: %s" % dst)
    dst.symlink_to(src.resolve(), target_is_directory=src.is_dir())


def check_collection_summary():
    summary = ORIGINAL / "updated_collection_summary.tsv"
    if not summary.is_file():
        raise RuntimeError("Missing sample summary: %s" % summary)
    with summary.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or ())
        required = {"sample", "mt_mean_coverage", "nuclear_median_coverage"}
        missing = required - fields
        if missing:
            raise RuntimeError(
                "Collection summary lacks required mean-mtCN inputs: %s"
                % ", ".join(sorted(missing))
            )
        seen = total = 0
        for row in reader:
            total += 1
            if row.get("mt_mean_coverage") not in (None, "", "NA") and row.get("nuclear_median_coverage") not in (None, "", "NA"):
                seen += 1
    print("Collection summary: %d samples, %d have both mtDNA mean and nuclear median depth" % (total, seen), flush=True)
    if seen == 0:
        raise RuntimeError("No samples have both inputs for mean-based mtCN")


def check_repaired_vcfs():
    old = ORIGINAL / "rrna_match" / RRNA_FOLDER
    repaired = DP50 / "rrna_match" / RRNA_FOLDER
    srcs = sorted(x for x in old.glob("*.lifted.codon.trna.rrna.vcf") if x.is_file())
    if not srcs:
        raise RuntimeError("No original rRNA VCFs")
    absent = [p.name for p in srcs if not (repaired / p.name).is_file()]
    still_links = [p.name for p in srcs if (repaired / p.name).is_symlink()]
    if absent or still_links:
        raise RuntimeError(
            "Repaired DP50 rRNA VCFs incomplete: missing=%d symlinks=%d (examples %s %s)"
            % (len(absent), len(still_links), absent[:3], still_links[:3])
        )
    print("Repaired rRNA VCFs available: %d" % len(srcs), flush=True)


def prepare(cohort):
    spec = COHORTS[cohort]
    new = spec["root"]
    if not spec["config"].is_file():
        raise RuntimeError("Missing cohort config: %s" % spec["config"])
    new.mkdir(parents=True, exist_ok=True)
    for src in sorted(ORIGINAL.iterdir()):
        dst = new / src.name
        if src.name in REBUILD:
            if dst.is_symlink():
                raise RuntimeError("Output directory is a symlink: %s" % dst)
            dst.mkdir(exist_ok=True)
        elif src.name == "pre_liftover_variant_qc":
            if cohort == "hom10":
                symlink(DP50 / src.name, dst)
            else:
                if dst.is_symlink():
                    raise RuntimeError("HOM+HET DP100 prelift must not be a symlink: %s" % dst)
                dst.mkdir(exist_ok=True)
                (dst / "reports").mkdir(exist_ok=True)
                (dst / "vcf_source_qc").mkdir(exist_ok=True)
        elif src.name == "rrna_match":
            if dst.is_symlink():
                raise RuntimeError("Refusing to write into linked rRNA annotation directory: %s" % dst)
            dst.mkdir(exist_ok=True)
            for sub in sorted(src.iterdir()):
                if sub.name == RRNA_FOLDER:
                    symlink(DP50 / "rrna_match" / RRNA_FOLDER, dst / RRNA_FOLDER)
                else:
                    symlink(sub, dst / sub.name)
        else:
            symlink(src, dst)

    for name in REBUILD:
        target = new / name
        if not target.exists():
            target.mkdir()
        if target.is_symlink():
            raise RuntimeError("Unsafe symlink for output %s" % target)
    (new / "slurm_logs").mkdir(exist_ok=True)
    print("READY %s: %s" % (cohort, new), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", choices=["hom10", "homhet40", "both"], default="both")
    args = parser.parse_args()
    check_collection_summary()
    check_repaired_vcfs()
    for cohort in (COHORTS if args.cohort == "both" else (args.cohort,)):
        prepare(cohort)


if __name__ == "__main__":
    main()
