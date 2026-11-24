#!/usr/bin/env python3
"""NUMT-to-mtDNA BASE coverage only (no variants involved).

Use upstream per-sample NUMT best-hit TSVs for chrM coordinates, optionally
restricted to nuclear best-hit regions overlapping highconf NUMT BEDs.
Calculate unique mtDNA base coverage for each QC-passed sample and for
species/reference unions. High-confidence sensitivity tiers (>=10/50/90%
nuclear best-hit overlap) complement the production any-overlap rule.

The BED coordinate convention in the existing qc_analysis.lib.numt_annotation
is retained for production-comparable classification (inclusive endpoints).
Highconf BED coordinates are NOT directly interpreted as mtDNA positions.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from qc_analysis.lib.numt_annotation import read_besthit, read_highconf_bed
from qc_analysis.scripts.run_local_heteroplasmy_qc import read_metadata, patterns, choose_numt_dir
from qc_analysis.lib.simple_yaml import read_simple_yaml

COHORTS = {
    "hom10": ("config/qc_cohort_hom10_mean_mtcn.yaml", "results/qc_cohort_hom10_mean_mtcn"),
    "homhet40": ("config/qc_cohort_homhet40_mean_mtcn.yaml", "results/qc_cohort_homhet40_mean_mtcn"),
}
TIERS = ("any_besthit", "highconf_any_overlap", "highconf_nuclear_10pct",
         "highconf_nuclear_50pct", "highconf_nuclear_90pct")


def path_from_root(x):
    p = Path(x).expanduser()
    return p if p.is_absolute() else ROOT / p


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as f:
        yield from csv.DictReader(f, delimiter="\t")


def write_tsv(path, records, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def find_sample_file(directory, pats, sample):
    for pat in pats:
        p = directory / pat.format(sample=sample)
        if p.is_file():
            return p
    return None


def unique_bases(intervals, size):
    spans = sorted((max(1, min(a, b)), min(size, max(a, b)))
                   for a, b in intervals)
    covered = 0
    lo = hi = None
    for a, b in spans:
        if a > b:
            continue
        if lo is None:
            lo, hi = a, b
        elif a <= hi + 1:
            hi = max(hi, b)
        else:
            covered += hi - lo + 1
            lo, hi = a, b
    if lo is not None:
        covered += hi - lo + 1
    return covered


def nuclear_fraction(chrom, start, end, nuclear_bed):
    """Best-hit nuclear bases overlapped by highconf BED (union, capped at 1).

    Replicate original qc_analysis inclusive-overlap interpretation for
    the production 'any overlap' tier. For precise nucleotide accounting,
    validate the upstream BED and best-hit coordinate conventions.
    """
    pieces = [(max(start, bstart), min(end, bend))
              for bchrom, bstart, bend in nuclear_bed
              if bchrom == chrom and max(start, bstart) <= min(end, bend)]
    return unique_bases(pieces, end) / (end - start + 1) if pieces else 0.0


def finite_int(raw):
    try:
        return int(float(raw))
    except (TypeError, ValueError):
        return None


def length_map(filepath):
    lengths = {}
    if not filepath:
        return lengths
    for row in read_tsv(path_from_root(filepath)):
        k = (row["species"], row["reference_key"])
        n = finite_int(row["mt_length"])
        if n is None or n < 1000:
            raise ValueError(f"Invalid mt_length for {k}: {row.get('mt_length')}")
        lengths[k] = n
    return lengths


def quantile(values, q):
    values = sorted(values)
    if not values:
        return ""
    x = (len(values)-1)*q
    low = int(x)
    upper = min(low+1, len(values)-1)
    return values[low]*(1-(x-low))+values[upper]*(x-low)


def plot_results(outdir, sample_rows, species_rows):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib unavailable: TSV files are available, plots skipped.")
        return
    for scope, values, field in (
        ("sample_specific", sample_rows, "coverage_pct"),
        ("species_union", species_rows, "coverage_pct"),
    ):
        fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
        for i, cohort in enumerate(COHORTS):
            for j, tier in enumerate(("any_besthit", "highconf_any_overlap")):
                ax = axes[i][j]
                data = [float(r[field]) for r in values
                        if r["cohort"] == cohort and r["tier"] == tier]
                ax.hist(data, bins=range(0, 101, 5), edgecolor="white")
                if data:
                    ax.axvline(median(data), ls="--", lw=1.2)
                ax.set_title(f"{cohort} / {tier} (n={len(data)})")
                if i == 1:
                    ax.set_xlabel("mtDNA bases covered by NUMT regions (%)")
                if j == 0:
                    ax.set_ylabel("Samples" if scope == "sample_specific" else
                                  "Species/reference catalogues")
                ax.set_xlim(0, 100)
        fig.suptitle(f"NUMT base coverage: {scope.replace('_', ' ')}")
        fig.tight_layout()
        fig.savefig(outdir / f"{scope}_numt_coverage.png", dpi=240)
        plt.close(fig)

    # Sensitivity of sample-specific highconf definition; avoids treating all
    # nuclear best-hit intervals as equally high confidence.
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.7), sharey=True)
    chosen = TIERS[1:]
    for ax, cohort in zip(axes, COHORTS):
        distributions = [[float(r["coverage_pct"]) for r in sample_rows
                          if r["cohort"] == cohort and r["tier"] == tier]
                         for tier in chosen]
        if all(distributions):
            ax.boxplot(distributions, showfliers=False)
        ax.set_xticks(range(1, 5), ["Any overlap", ">=10%", ">=50%", ">=90%"],
                      rotation=20)
        ax.set_ylim(0, 100)
        ax.set_title(cohort)
        ax.set_ylabel("mtDNA coverage (%)")
    fig.suptitle("High-confidence nuclear-overlap threshold sensitivity")
    fig.tight_layout()
    fig.savefig(outdir / "highconf_threshold_sensitivity.png", dpi=240)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", choices=("both", "hom10", "homhet40"), default="both")
    ap.add_argument("--besthit-dir", type=Path)
    ap.add_argument("--highconf-dir", type=Path)
    ap.add_argument("--lengths-tsv", type=Path,
                    help="species/reference-specific lengths: species,reference_key,mt_length")
    ap.add_argument("--mt-length", type=int, default=16569,
                    help="PROXY denominator only; not valid for exact cross-species analysis")
    ap.add_argument("--output", type=Path,
                    default=ROOT / "results/numt_region_base_coverage")
    args = ap.parse_args()
    out = path_from_root(args.output)
    out.mkdir(parents=True, exist_ok=True)
    lengths = length_map(args.lengths_tsv)
    sample_rows, species_rows, qc_rows, summary_rows = [], [], [], []
    selected = COHORTS if args.cohort == "both" else (args.cohort,)

    for cohort in selected:
        cfgname, result_root = COHORTS[cohort]
        cfg = read_simple_yaml(path_from_root(cfgname))
        sec = cfg.get("local_heteroplasmy_qc", {})
        metadata_file = path_from_root(sec.get("sample_ref_file", "config/sample_ref_file.tsv"))
        species_map, ref_map = read_metadata(metadata_file)
        sample_qc = path_from_root(result_root) / "sample_variant_filtering/reports/sample_qc.tsv"
        if not sample_qc.exists():
            raise FileNotFoundError(sample_qc)
        qc_pass = sorted({r["sample"] for r in read_tsv(sample_qc)
                          if r.get("qc_status") == "PASS"})
        samples = [s for s in qc_pass if species_map.get(s)]
        best_pats = patterns(sec.get("numt_besthit_pattern"), [
            "{sample}.numt_vs_chrM.besthit.tsv", "{sample}.numt_besthit.tsv",
            "{sample}.besthit.tsv", "{sample}.tsv"])
        high_pats = patterns(sec.get("numt_highconf_pattern"),
                             ["{sample}.highconf_numt.bed"])
        fallback = ROOT / "data/numt_besthit"
        bestdir = (path_from_root(args.besthit_dir) if args.besthit_dir
                   else choose_numt_dir(sec.get("numt_besthit_dir"), best_pats,
                                        samples, [fallback]))
        highdir = (path_from_root(args.highconf_dir) if args.highconf_dir
                   else choose_numt_dir(sec.get("numt_highconf_dir"), high_pats,
                                        samples, [bestdir, fallback]))
        union_by_ref = defaultdict(lambda: defaultdict(list))
        counters = defaultdict(int)

        for sample in samples:
            species, ref = species_map[sample], ref_map.get(sample, species)
            key = (species, ref)
            bestfile = find_sample_file(bestdir, best_pats, sample)
            highfile = find_sample_file(highdir, high_pats, sample)
            if bestfile is None:
                qc_rows.append(dict(cohort=cohort, sample=sample, species=species,
                                    reference_key=ref, status="MISSING_BESTHIT",
                                    besthit_file="", highconf_file=""))
                counters["missing_besthit"] += 1
                continue
            if highfile is None:
                qc_rows.append(dict(cohort=cohort, sample=sample, species=species,
                                    reference_key=ref, status="MISSING_HIGHCONF_BED",
                                    besthit_file=str(bestfile), highconf_file=""))
                counters["missing_highconf"] += 1
            highbed = read_highconf_bed(highfile) if highfile else None
            by_tier = defaultdict(list)
            n_records = 0
            for row in read_besthit(bestfile):
                try:
                    coords = [finite_int(row[k]) for k in (
                        "nuclear_start", "nuclear_end", "chrM_start", "chrM_end")]
                except KeyError:
                    continue
                if any(x is None for x in coords):
                    continue
                ns, ne = sorted(coords[:2])
                ms, me = sorted(coords[2:])
                n_records += 1
                by_tier["any_besthit"].append((ms, me))
                if highbed is not None:
                    frac = nuclear_fraction(row.get("nuclear_chrom", ""), ns, ne, highbed)
                    if frac > 0:
                        by_tier["highconf_any_overlap"].append((ms, me))
                    if frac >= .10:
                        by_tier["highconf_nuclear_10pct"].append((ms, me))
                    if frac >= .50:
                        by_tier["highconf_nuclear_50pct"].append((ms, me))
                    if frac >= .90:
                        by_tier["highconf_nuclear_90pct"].append((ms, me))
            if n_records == 0:
                qc_rows.append(dict(cohort=cohort, sample=sample, species=species,
                                    reference_key=ref, status="NO_VALID_BESTHIT_ROWS",
                                    besthit_file=str(bestfile),
                                    highconf_file=str(highfile or "")))
                counters["no_valid_besthit"] += 1
                continue
            counters["samples_with_besthit"] += 1
            size = lengths.get(key, args.mt_length)
            if size <= 0:
                raise ValueError("Nonpositive mtDNA length")
            for tier in TIERS:
                if tier != "any_besthit" and highbed is None:
                    continue
                ivals = by_tier[tier]
                bp = unique_bases(ivals, size)
                union_by_ref[key][tier].extend(ivals)
                sample_rows.append(dict(cohort=cohort, sample=sample,
                    species=species, reference_key=ref, tier=tier,
                    mt_length=size, length_source="reference_specific" if key in lengths else "fixed_proxy",
                    covered_bp=bp, coverage_pct=round(100*bp/size, 4),
                    n_besthit_rows=n_records))
        for key, tiers in sorted(union_by_ref.items()):
            size = lengths.get(key, args.mt_length)
            for tier in TIERS:
                if tier not in tiers:
                    continue
                bp = unique_bases(tiers[tier], size)
                species_rows.append(dict(cohort=cohort, species=key[0],
                    reference_key=key[1], tier=tier, mt_length=size,
                    length_source="reference_specific" if key in lengths else "fixed_proxy",
                    covered_bp=bp, coverage_pct=round(100*bp/size, 4)))

        print(cohort, "QC-passed samples:", len(qc_pass),
              "with metadata:", len(samples), "with besthit:", counters["samples_with_besthit"],
              "missing besthit:", counters["missing_besthit"],
              "missing highconf:", counters["missing_highconf"])
        if not counters["samples_with_besthit"]:
            raise RuntimeError(f"{cohort}: no sample-specific besthit files matched in {bestdir}. "
                               "Use --besthit-dir and --highconf-dir")

    fields_s = ["cohort","sample","species","reference_key","tier","mt_length",
                "length_source","covered_bp","coverage_pct","n_besthit_rows"]
    fields_u = ["cohort","species","reference_key","tier","mt_length",
                "length_source","covered_bp","coverage_pct"]
    write_tsv(out/"sample_numt_region_base_coverage.tsv",sample_rows,fields_s)
    write_tsv(out/"species_union_numt_region_base_coverage.tsv",species_rows,fields_u)
    write_tsv(out/"numt_input_diagnostics.tsv",qc_rows,
              ["cohort","sample","species","reference_key","status","besthit_file","highconf_file"])
    for scope, tab in (("sample",sample_rows),("species_union",species_rows)):
        for cohort in selected:
            for tier in TIERS:
                vals = [r["coverage_pct"] for r in tab
                        if r["cohort"] == cohort and r["tier"] == tier]
                if not vals:
                    continue
                vals = sorted(vals)
                summary_rows.append(dict(cohort=cohort,scope=scope,tier=tier,n=len(vals),
                   median_pct=round(median(vals),4),q1_pct=round(quantile(vals,.25),4),
                   q3_pct=round(quantile(vals,.75),4),
                   n_coverage_ge90=sum(v>=90 for v in vals),
                   n_coverage_ge99=sum(v>=99 for v in vals)))
    write_tsv(out/"numt_region_base_coverage_summary.tsv",summary_rows,
              ["cohort","scope","tier","n","median_pct","q1_pct","q3_pct",
               "n_coverage_ge90","n_coverage_ge99"])
    for r in summary_rows:
        print("{cohort:8} {scope:13} {tier:24} n={n:5} median={median_pct:7.2f}% "
              "IQR={q1_pct:5.1f}-{q3_pct:5.1f}%".format(**r))
    if not lengths:
        print("WARNING: all fractions use proxy mtDNA length of", args.mt_length,
              "bp; provide --lengths-tsv for reference-specific denominators.",
              file=sys.stderr)
    plot_results(out,sample_rows,species_rows)
    print("OUTPUT:",out)

if __name__ == "__main__":
    main()
