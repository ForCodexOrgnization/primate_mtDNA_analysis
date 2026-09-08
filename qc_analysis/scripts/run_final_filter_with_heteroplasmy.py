#!/usr/bin/env python3
"""Run terminal QC filtering and apply native-coordinate heteroplasmy/NUMT removals.

The base final_filter remains the owner of sample-level contamination/QC decisions.
This wrapper then consumes local_heteroplasmy_qc/reports/numt_variants_to_remove.tsv
and removes exactly those variants by immutable SOURCE_* identity.  NUMT-associated
samples are annotated, not failed wholesale by default.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from qc_analysis.lib.simple_yaml import read_simple_yaml


def resolve(value) -> Path:
    path = Path(str(value)).expanduser()
    return path if path.is_absolute() else ROOT / path


def read_tsv(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv_atomic(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix=path.name + ".", dir=path.parent)
    os.close(fd)
    temp = Path(temp_name)
    try:
        with temp.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        temp.replace(path)
    finally:
        if temp.exists():
            temp.unlink()


def source_key(row: dict) -> tuple[str, str, str, str, str]:
    return (
        str(row.get("sample", "")),
        str(row.get("source_chrom", "")),
        str(row.get("source_pos", "")),
        str(row.get("source_ref", "")),
        str(row.get("source_alt", "")),
    )


def parse_info(value: str) -> dict[str, str]:
    out = {}
    for item in value.split(";"):
        if not item or item == ".":
            continue
        if "=" in item:
            key, val = item.split("=", 1)
            out[key] = val
    return out


def append_reason(existing: str, reason: str) -> str:
    existing = "" if existing in {"", ".", "PASS", "NOT_AVAILABLE"} else existing
    parts = [item for item in existing.split(";") if item]
    if reason not in parts:
        parts.append(reason)
    return ";".join(parts)


def filter_vcf(path: Path, sample: str, removal: dict[tuple, dict]) -> int:
    """Remove SOURCE-keyed variants and regenerate tabix index when possible."""
    opener = gzip.open if path.suffix == ".gz" else open
    removed = 0
    with tempfile.TemporaryDirectory(prefix="heteroplasmy_filter_", dir=path.parent) as tmpdir:
        plain = Path(tmpdir) / (path.name[:-3] if path.name.endswith(".gz") else path.name)
        with opener(path, "rt") as src, plain.open("w", encoding="utf-8") as dst:
            for line in src:
                if line.startswith("#"):
                    dst.write(line)
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    dst.write(line)
                    continue
                info = parse_info(fields[7])
                key = (
                    sample,
                    info.get("SOURCE_CHROM", ""),
                    info.get("SOURCE_POS", ""),
                    info.get("SOURCE_REF", ""),
                    info.get("SOURCE_ALT", ""),
                )
                if key in removal:
                    removed += 1
                    continue
                dst.write(line)

        if path.suffix == ".gz":
            try:
                import pysam  # type: ignore

                bgz = Path(tmpdir) / path.name
                pysam.tabix_compress(str(plain), str(bgz), force=True)
                pysam.tabix_index(str(bgz), preset="vcf", force=True)
                shutil.copy2(bgz, path)
                tbi_src = Path(str(bgz) + ".tbi")
                if tbi_src.is_file():
                    shutil.copy2(tbi_src, Path(str(path) + ".tbi"))
            except ImportError:
                bgzip = shutil.which("bgzip")
                tabix = shutil.which("tabix")
                if not bgzip or not tabix:
                    raise RuntimeError("pysam or bgzip+tabix is required to rewrite final .vcf.gz outputs")
                compressed = subprocess.check_output([bgzip, "-c", str(plain)])
                path.write_bytes(compressed)
                subprocess.run([tabix, "-f", "-p", "vcf", str(path)], check=True)
        else:
            shutil.copy2(plain, path)
    return removed


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--skip-base-final-filter", action="store_true")
    args = parser.parse_args()

    cfg = read_simple_yaml(args.config)
    final_sec = cfg.get("final_filter", {})
    het_sec = cfg.get("local_heteroplasmy_qc", {})

    if not args.skip_base_final_filter:
        subprocess.run(
            [sys.executable, str(ROOT / "qc_analysis/scripts/run_final_filter_streaming.py"), "--config", str(args.config)],
            check=True,
        )

    final_dir = resolve(final_sec.get("output_dir", "results/qc/final_filter"))
    het_dir = resolve(het_sec.get("output_dir", "results/qc/local_heteroplasmy_qc")) / "reports"
    removal_path = resolve(final_sec.get("heteroplasmy_variant_removal_report", het_dir / "numt_variants_to_remove.tsv"))
    sample_het_path = resolve(final_sec.get("heteroplasmy_sample_report", het_dir / "local_heteroplasmy_sample_summary.tsv"))

    removal_rows = [row for row in read_tsv(removal_path) if row.get("filter_action", "").upper() == "REMOVE"]
    removal = {source_key(row): row for row in removal_rows}
    sample_het = {row.get("sample", ""): row for row in read_tsv(sample_het_path)}

    variant_report = final_dir / "reports" / "final_variant_qc.tsv"
    variants = read_tsv(variant_report)
    removed_report_rows = 0
    if variants:
        base_fields = list(variants[0])
        extra = ["heteroplasmy_filter_status", "heteroplasmy_filter_reason", "heteroplasmy_cluster_id", "heteroplasmy_numt_scope", "heteroplasmy_numt_tier"]
        fields = base_fields + [field for field in extra if field not in base_fields]
        for row in variants:
            hit = removal.get(source_key(row))
            if hit:
                row["heteroplasmy_filter_status"] = "REMOVE"
                row["heteroplasmy_filter_reason"] = hit.get("filter_reason", "NUMT_ASSOCIATED_HETEROPLASMY")
                row["heteroplasmy_cluster_id"] = hit.get("cluster_id", "")
                row["heteroplasmy_numt_scope"] = hit.get("numt_scope", "")
                row["heteroplasmy_numt_tier"] = hit.get("numt_tier", "")
                row["final_variant_status"] = "FAIL"
                row["final_variant_fail_reasons"] = append_reason(
                    row.get("final_variant_fail_reasons", ""),
                    "heteroplasmy_numt:" + hit.get("filter_reason", "REMOVE"),
                )
                removed_report_rows += 1
            else:
                row["heteroplasmy_filter_status"] = "KEEP"
                row["heteroplasmy_filter_reason"] = ""
                row["heteroplasmy_cluster_id"] = ""
                row["heteroplasmy_numt_scope"] = ""
                row["heteroplasmy_numt_tier"] = ""
        write_tsv_atomic(variant_report, variants, fields)

    sample_report = final_dir / "reports" / "final_sample_qc.tsv"
    samples = read_tsv(sample_report)
    if samples:
        base_fields = list(samples[0])
        extra = ["numt_sample", "n_numt_clusters", "n_numt_variants_removed", "fraction_het_removed"]
        fields = base_fields + [field for field in extra if field not in base_fields]
        for row in samples:
            het = sample_het.get(row.get("sample", ""), {})
            row["numt_sample"] = het.get("numt_sample", "NO")
            row["n_numt_clusters"] = het.get("n_numt_clusters", "0")
            row["n_numt_variants_removed"] = het.get("n_numt_variants_to_remove", "0")
            row["fraction_het_removed"] = het.get("fraction_het_removed", "0")
        write_tsv_atomic(sample_report, samples, fields)

    removed_vcf_records = 0
    final_vcf_dir = final_dir / "final_vcf"
    for path in sorted(final_vcf_dir.glob("*.vcf*")):
        if path.name.endswith((".tbi", ".csi")):
            continue
        sample = path.name.split(".vcf", 1)[0]
        removed_vcf_records += filter_vcf(path, sample, removal)

    summary = [{
        "blacklisted_source_variants": len(removal),
        "final_variant_report_rows_failed_by_heteroplasmy": removed_report_rows,
        "final_vcf_records_removed": removed_vcf_records,
        "numt_samples": sum(1 for row in sample_het.values() if row.get("numt_sample") == "YES"),
    }]
    write_tsv_atomic(
        final_dir / "reports" / "heteroplasmy_final_filter_summary.tsv",
        summary,
        ["blacklisted_source_variants", "final_variant_report_rows_failed_by_heteroplasmy", "final_vcf_records_removed", "numt_samples"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
