#!/usr/bin/env python3
"""Apply homoplasmic-only sample gates to PREVIOUS final QC results.

Does not modify original final_filter outputs or heteroplasmic statuses.
Uses the staged collection summary + previously merged per-base coverage.
Run after final_filter; rerunning is deterministic.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import math
import shutil
import subprocess
import tempfile
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from qc_analysis.lib.simple_yaml import read_simple_yaml


def path_of(x):
    p = Path(str(x)).expanduser()
    return p if p.is_absolute() else ROOT / p


def rows_from(path):
    with path.open(newline="", encoding="utf-8") as f:
        yield from csv.DictReader(f, delimiter="\t")


def number(raw):
    try:
        v = float(raw)
        return v if math.isfinite(v) else None
    except (ValueError, TypeError):
        return None


def percent_50_from_merged_coverage(path, min_positions):
    """Use same merged per-base denominator as original Percent_100."""
    if not path.is_file():
        return None, 0, "missing_merged_coverage"
    positions = set()
    covered = set()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"chrom", "pos", "target", "coverage"}.issubset(reader.fieldnames):
            raise ValueError(f"Unexpected merged coverage layout: {path}")
        for row in reader:
            key = (row["chrom"], row["pos"], row["target"])
            depth = number(row["coverage"])
            if depth is None:
                continue
            positions.add(key)
            if depth >= 50:
                covered.add(key)
    n = len(positions)
    if n < min_positions:
        return None, n, "too_few_coverage_positions"
    return 100 * len(covered) / n, n, "ok"


def write_report(path, columns, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="") as h:
        writer = csv.DictWriter(h, delimiter="\t", fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(path)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", required=True)
    ap.add_argument("--write-vcf", action="store_true", help="Also emit new VCFs, without changing original final_vcf")
    args = ap.parse_args()
    cfg = read_simple_yaml(args.config)
    sec = cfg.get("homoplasmy_qc") or {}
    if sec.get("enabled", True) is False:
        print("[homoplasmy_qc] disabled")
        return 0
    sample_sec = cfg.get("sample_variant_filtering") or {}
    collected = path_of((cfg.get("final_filter") or {}).get("collected_dir",
                         "results/qc/collected_variant_calling_results"))
    summary = path_of(sec.get("input_summary", sample_sec.get("input_summary",
                         str(collected / "reports/variant_calling_collection_summary.tsv"))))
    final_dir = path_of((cfg.get("final_filter") or {}).get("output_dir", "results/qc/final_filter"))
    final_report = final_dir / "reports/final_variant_qc.tsv"
    out = path_of(sec.get("output_dir", str(final_dir.parent / "homoplasmy_qc")))
    threshold_mtcn = float(sec.get("min_mtcn", 20))
    threshold_percent = float(sec.get("min_percent_50", 50))
    min_positions = int(sec.get("min_coverage_positions", 14000))
    if not summary.is_file() or not final_report.is_file():
        raise FileNotFoundError(f"Need existing summary and final QC: {summary}, {final_report}")

    evaluations = {}
    sample_rows = []
    for sample in rows_from(summary):
        name = sample.get("sample") or sample.get("Sample")
        if not name or name in evaluations:
            if name in evaluations:
                raise ValueError(f"Duplicate sample in summary: {name}")
            continue
        mt = number(sample.get("mt_median_coverage"))
        nuclear = number(sample.get("nuclear_median_coverage"))
        original = number(sample.get("mtcn_median"))
        nuclear_adjusted = 1.0 if nuclear == 0 else nuclear
        recalculated = (2.0 * mt / nuclear_adjusted
                        if mt is not None and nuclear_adjusted is not None and nuclear_adjusted > 0 else None)
        cov_path = collected / "collected_cov" / f"{name}.merged.max_depth.per_base_coverage.tsv"
        pct, n_positions, coverage_status = percent_50_from_merged_coverage(cov_path, min_positions)
        failures = []
        if recalculated is None or recalculated < threshold_mtcn:
            failures.append("HOM_mtCN_below_20_or_missing")
        if pct is None or pct < threshold_percent:
            failures.append("HOM_Percent_50_below_50_or_missing")
        result = {
            "sample": name, "species": sample.get("species", ""),
            "mt_median_coverage": mt, "nuclear_median_coverage_original": nuclear,
            "nuclear_median_coverage_adjusted": nuclear_adjusted,
            "mtcn_median_original": original, "mtcn_median_recalculated": recalculated,
            "Percent_50": pct, "n_coverage_positions": n_positions,
            "coverage_status": coverage_status,
            "homoplasmy_qc_status": "FAIL" if failures else "PASS",
            "homoplasmy_qc_fail_reasons": ";".join(failures)
        }
        evaluations[name] = result
        sample_rows.append(result)

    columns = list(sample_rows[0]) if sample_rows else []
    if not columns:
        raise ValueError("Empty sample collection summary")
    write_report(out / "reports/homoplasmy_sample_qc.tsv", columns, sample_rows)

    # Keep full previous final QC table, append new gate fields.
    # Only previous PASS calls with authoritative SOURCE HOM are affected.
    with final_report.open(newline="", encoding="utf-8") as h:
        reader = csv.DictReader(h, delimiter="\t")
        base_fields = reader.fieldnames or []
        for required in ("sample", "source_call_class", "final_variant_status"):
            if required not in base_fields:
                raise ValueError(f"Missing column in final QC: {required}")
        extra = ["homoplasmy_qc_status", "homoplasmy_qc_fail_reasons",
                 "mtcn_median_recalculated", "Percent_50",
                 "homoplasmy_adjusted_final_status", "homoplasmy_adjusted_fail_reasons"]
        path = out / "reports/final_variant_qc.homoplasmy_adjusted.tsv"
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(path.name + ".tmp")
        counts = {"hom_pass": 0, "hom_new_fail": 0, "other_unchanged": 0}
        with tmp.open("w", encoding="utf-8", newline="") as dest:
            writer = csv.DictWriter(dest, fieldnames=base_fields + [x for x in extra if x not in base_fields],
                                    delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            for row in reader:
                sample_result = evaluations.get(row["sample"])
                if not sample_result:
                    raise ValueError(f"Sample missing from source summary: {row['sample']}")
                is_hom = row["source_call_class"].upper() == "HOM"
                prior = row["final_variant_status"]
                new_fail = is_hom and prior == "PASS" and sample_result["homoplasmy_qc_status"] != "PASS"
                row["homoplasmy_qc_status"] = sample_result["homoplasmy_qc_status"]
                row["homoplasmy_qc_fail_reasons"] = sample_result["homoplasmy_qc_fail_reasons"]
                row["mtcn_median_recalculated"] = sample_result["mtcn_median_recalculated"]
                row["Percent_50"] = sample_result["Percent_50"]
                row["homoplasmy_adjusted_final_status"] = "FAIL" if new_fail else prior
                previous = row.get("final_variant_fail_reasons", "")
                reason = sample_result["homoplasmy_qc_fail_reasons"] if new_fail else ""
                row["homoplasmy_adjusted_fail_reasons"] = ";".join(x for x in (previous, reason) if x and x not in ("PASS", "NOT_AVAILABLE"))
                counts["hom_new_fail" if new_fail else "hom_pass" if is_hom and prior == "PASS" else "other_unchanged"] += 1
                writer.writerow(row)
        tmp.replace(path)

    if args.write_vcf:
        destdir = out / "final_vcf"
        destdir.mkdir(parents=True, exist_ok=True)
        try:
            import pysam
        except ImportError:
            pysam = None
        if not pysam and not (shutil.which("bgzip") and shutil.which("tabix")):
            raise RuntimeError("VCF export needs pysam or both bgzip and tabix")
        for source in sorted((final_dir / "final_vcf").glob("*.final.vcf.gz")):
            name = source.name.split(".final.vcf.gz")[0]
            assessment = evaluations.get(name)
            if assessment is None:
                raise ValueError(f"VCF sample missing summary: {name}")
            with tempfile.TemporaryDirectory(dir=destdir) as tmp_dir:
                uncompressed = Path(tmp_dir) / "filtered.vcf"
                with gzip.open(source, "rt") as inp, uncompressed.open("w") as out_vcf:
                    for line in inp:
                        if line.startswith("#"):
                            out_vcf.write(line)
                            continue
                        if assessment["homoplasmy_qc_status"] == "PASS":
                            out_vcf.write(line)
                        else:
                            fields = line.rstrip("\n").split("\t")
                            info = dict(i.split("=", 1) for i in fields[7].split(";") if "=" in i)
                            if info.get("SOURCE_CALL_CLASS", "").upper() != "HOM":
                                out_vcf.write(line)
                dest = destdir / source.name
                if pysam:
                    pysam.tabix_compress(str(uncompressed), str(dest), force=True)
                    pysam.tabix_index(str(dest), preset="vcf", force=True)
                else:
                    with dest.open("wb") as h:
                        subprocess.run(["bgzip", "-c", str(uncompressed)], stdout=h, check=True)
                    subprocess.run(["tabix", "-f", "-p", "vcf", str(dest)], check=True)

    print(f"[homoplasmy_qc] sample_report={out / 'reports/homoplasmy_sample_qc.tsv'}")
    print(f"[homoplasmy_qc] variant_report={out / 'reports/final_variant_qc.homoplasmy_adjusted.tsv'}")
    print(f"[homoplasmy_qc] sample_pass={sum(x['homoplasmy_qc_status']=='PASS' for x in sample_rows)}/{len(sample_rows)} counts={counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
