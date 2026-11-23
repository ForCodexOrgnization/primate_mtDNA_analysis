#!/usr/bin/env python3
"""Summarize mean-mtCN cohort QC and plot mtDNA mean coverage distributions."""
from __future__ import annotations

import csv
from pathlib import Path
from statistics import median
import sys

ROOT = Path(__file__).resolve().parents[2]
HOM = ROOT / "results/qc_cohort_hom10_mean_mtcn"
HET = ROOT / "results/qc_cohort_homhet40_mean_mtcn"


def rows(path):
    if not path.is_file():
        raise RuntimeError("Missing report: " + str(path))
    with path.open(newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def num(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return None


def sample_counts(cohort):
    path = cohort / "sample_variant_filtering/reports/sample_qc.tsv"
    records = list(rows(path))
    pass_col = "homoplasmy_qc_status" if cohort == HOM else "qc_status"
    passed = [r for r in records if r.get(pass_col) == "PASS"]
    samples = {r["sample"] for r in passed}
    species = {r["species"] for r in passed if r.get("species") not in (None, "", "NA")}
    return records, passed, samples, species


def final_counts(cohort):
    path = cohort / "final_filter_hom_only/reports/final_variant_qc.tsv"
    total = {"HOM": 0, "HET": 0, "OTHER": 0}
    samples, species = set(), set()
    for row in rows(path):
        if row.get("final_variant_status") != "PASS":
            continue
        cls = row.get("source_call_class", "OTHER")
        total[cls if cls in total else "OTHER"] += 1
        samples.add(row["sample"])
        if row.get("species") not in (None, "", "NA"):
            species.add(row["species"])
    return total, len(samples), len(species)


def main():
    summary = []
    cohort_data = []
    for name, cohort in (("HOM10", HOM), ("HOMHET40", HET)):
        records, passed, samples, species = sample_counts(cohort)
        variants, variant_samples, variant_species = final_counts(cohort)
        cohort_data.append((name, cohort, records, passed))
        summary.append({
            "cohort": name,
            "all_samples": len(records),
            "pass_samples": len(samples),
            "pass_species": len(species),
            "pass_variant_samples": variant_samples,
            "pass_variant_species": variant_species,
            "pass_hom_records": variants["HOM"],
            "pass_het_records": variants["HET"],
            "pass_other_records": variants["OTHER"],
        })

    output = HOM / "reports"
    output.mkdir(parents=True, exist_ok=True)
    summary_path = output / "mean_mtcn_cohorts_comparison.tsv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)
    print("COHORT SUMMARY")
    for row in summary:
        print(row)
    print("Wrote", summary_path)

    _, _, hom_records, hom_pass = cohort_data[0]
    original20 = []
    rescued10 = []
    for row in hom_pass:
        m = num(row.get("mtcn_qc_value"))
        cov = num(row.get("mt_mean_coverage"))
        if m is None or cov is None:
            continue
        if 10 <= m < 20:
            rescued10.append(cov)
        elif m >= 20:
            original20.append(cov)
    print("HOM10 mean mtDNA coverage, non-overlapping strata:")
    for title, vals in (("10 <= mean-mtCN < 20", rescued10),
                        ("mean-mtCN >= 20", original20)):
        print("%s: n=%d median_mt_mean_coverage=%s"
              % (title, len(vals), median(vals) if vals else "NA"))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib unavailable; skipped coverage figure")
        return
    fig, ax = plt.subplots(figsize=(8, 5))
    data = [v for v in (rescued10, original20) if v]
    labels = [title for title, vals in (
        ("10 <= mtCN < 20", rescued10), ("mtCN >= 20", original20)) if vals]
    if data:
        ax.boxplot(data, tick_labels=labels, showfliers=False)
        ax.set_ylabel("Mean mtDNA coverage (x)")
        ax.set_title("Mean mtDNA coverage by mean-based mtCN threshold")
        fig.tight_layout()
        dest = output / "mean_mtDNA_coverage_mtCN10_vs20.png"
        fig.savefig(dest, dpi=250)
        print("Wrote", dest)
    plt.close(fig)


if __name__ == "__main__":
    main()
