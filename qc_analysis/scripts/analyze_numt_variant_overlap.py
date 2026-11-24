#!/usr/bin/env python3
"""Compare pre-filter source SNVs to final PASS SNVs at native mtDNA NUMT regions.

Restrict BOTH stages to the same final sample-QC PASS population.
Use source_pos and (species, reference_key) for interval overlap;
never compare lifted human_pos to native NUMT coordinates.
"""
from __future__ import annotations
import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
COHORTS = {
    "hom10": ROOT / "results/qc_cohort_hom10_mean_mtcn",
    "homhet40": ROOT / "results/qc_cohort_homhet40_mean_mtcn",
}

def tsv(path):
    with path.open(newline="", encoding="utf-8") as f:
        yield from csv.DictReader(f, delimiter="\t")

def metadata(path):
    rows = []
    with path.open(encoding="utf-8") as f:
        rows = [x.strip() for x in f if x.strip() and not x.startswith("#")]
    if not rows:
        raise ValueError("Empty metadata " + str(path))
    header = rows[0].split()
    is_header = bool(set(x.lower() for x in header) & {
        "sample", "species", "reference_key", "mt_reference", "chrm_reference"
    })
    result = {}
    if is_header:
        for row in csv.DictReader(rows, delimiter="\t" if "\t" in rows[0] else " "):
            sample = row.get("sample", row.get("Sample", ""))
            species = row.get("species", row.get("Species", ""))
            ref = next((str(row[k]) for k in (
                "reference_key", "mt_reference", "chrM_reference", "chrm_reference",
                "species_fasta", "reference_fasta"
            ) if row.get(k)), species)
            if sample and species:
                result[sample] = (species, ref)
    else:
        for line in rows:
            parts = line.split()
            if len(parts) >= 2:
                result[parts[0]] = (parts[1], parts[1])
    return result

def interval_index(path):
    data = defaultdict(lambda: defaultdict(list))
    for row in tsv(path):
        species, ref, tier = row["species"], row["reference_key"], row["tier"]
        try:
            a, b = int(row["chrm_start"]), int(row["chrm_end"])
        except ValueError:
            continue
        data[(species, ref)][tier].append((min(a, b), max(a, b)))
    if not data:
        raise ValueError("Empty NUMT catalogue: " + str(path))
    return data

def overlap(pos, data):
    high = any(a <= pos <= b for a, b in data.get("HIGH_CONF_NUMT", ()))
    any_hit = high or any(a <= pos <= b for tier, intervals in data.items()
                          for a, b in intervals)
    return int(any_hit), int(high)

def snv_key(row, sample, pos_col, ref_col, alt_col):
    try:
        pos = int(row[pos_col])
    except (KeyError, ValueError, TypeError):
        return None
    ref, alt = row.get(ref_col, "").upper(), row.get(alt_col, "").upper()
    if (pos < 1 or len(ref) != 1 or len(alt) != 1 or
        ref not in "ACGT" or alt not in "ACGT" or ref == alt):
        return None
    return sample, row.get("source_chrom", ""), pos, ref, alt

def load_before(path, eligible):
    output = {}
    for row in tsv(path):
        sample = row.get("sample", "")
        cls = row.get("source_call_class", "").upper()
        if sample not in eligible or cls not in ("HOM", "HET"):
            continue
        key = snv_key(row, sample, "source_pos", "source_ref", "source_alt")
        if key is not None:
            output[key] = cls
    return output

def load_after(path, eligible):
    output = {}
    for row in tsv(path):
        sample = row.get("sample", "")
        if sample not in eligible or row.get("final_variant_status") != "PASS":
            continue
        cls = row.get("source_call_class", "").upper()
        if cls not in ("HOM", "HET"):
            continue
        key = snv_key(row, sample, "source_pos", "source_ref", "source_alt")
        if key is not None:
            output[key] = cls
    return output

def summarize(cohort, root, meta, output_dir, intervals_override=None):
    sample_path = root / "final_filter_hom_only/reports/final_sample_qc.tsv"
    source_path = root / "pre_liftover_variant_qc/reports/source_variant_qc.tsv"
    final_path = root / "final_filter_hom_only/reports/final_variant_qc.tsv"
    intervals_path = intervals_override or root / "local_heteroplasmy_qc/reports/species_numt_intervals.tsv"
    for path in (sample_path, source_path, final_path, intervals_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    eligible = {r["sample"] for r in tsv(sample_path)
                if r.get("final_sample_status") == "PASS"}
    if not eligible:
        raise ValueError(f"{cohort}: no samples marked final_sample_status=PASS")
    catalog = interval_index(intervals_path)
    before = load_before(source_path, eligible)
    after = load_after(final_path, eligible)
    if not before:
        raise ValueError(f"{cohort}: no pre-filter source variants")
    unexpected = set(after) - set(before)
    if unexpected:
        raise ValueError(f"{cohort}: {len(unexpected)} final PASS keys missing from source QC; inspect provenance")

    grouped = defaultdict(Counter)
    unknown = Counter()
    detail = output_dir / f"{cohort}_numt_overlap_per_variant.tsv"
    with detail.open("w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(("cohort", "sample", "species", "reference_key", "source_chrom",
                         "source_pos", "source_ref", "source_alt", "call_class",
                         "final_status", "numt_besthit_overlap", "numt_highconf_overlap"))
        for key, cls in before.items():
            sample, chrom, pos, ref, alt = key
            if sample not in meta:
                unknown["missing_metadata"] += 1
                continue
            species, reference = meta[sample]
            if (species, reference) not in catalog:
                unknown["missing_interval_catalog"] += 1
                continue
            hit, high = overlap(pos, catalog[(species, reference)])
            stage = "PASS" if key in after else "REMOVED"
            group = "HOM" if cls == "HOM" else "HET"
            for label in ("ALL", group):
                for status in ("BEFORE", "AFTER" if stage == "PASS" else "REMOVED"):
                    if status == "BEFORE" or (status == "AFTER" and stage == "PASS") or (status == "REMOVED" and stage == "REMOVED"):
                        count = grouped[(label, status)]
                        count["n"] += 1
                        count["any"] += hit
                        count["high"] += high
            writer.writerow((cohort, sample, species, reference, chrom, pos, ref, alt,
                             group, stage, hit, high))
    summary_path = output_dir / f"{cohort}_numt_overlap_summary.tsv"
    with summary_path.open("w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(("cohort", "call_class", "stage", "variants", "any_numt_n",
                         "any_numt_pct", "highconf_numt_n", "highconf_numt_pct"))
        for cls in ("ALL", "HOM", "HET"):
            for stage in ("BEFORE", "AFTER", "REMOVED"):
                d = grouped[(cls, stage)]
                n = d["n"]
                writer.writerow((cohort, cls, stage, n, d["any"],
                                 f"{100*d['any']/n:.4f}" if n else "NA",
                                 d["high"], f"{100*d['high']/n:.4f}" if n else "NA"))
    print(f"{cohort}: final-PASS samples={len(eligible)}; raw unique calls={len(before)}, final PASS={len(after)}")
    print(f"  Interval-unassessable before calls: {dict(unknown)}")
    print(f"  Summary: {summary_path}")
    print(f"  Per-variant table: {detail}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", choices=("hom10", "homhet40", "both"), default="both")
    ap.add_argument("--metadata", type=Path, default=ROOT / "config/sample_ref_file.tsv")
    ap.add_argument("--output-dir", type=Path, default=ROOT / "results/numt_overlap_before_after")
    ap.add_argument("--intervals", type=Path, help="Optional NUMT species catalogue TSV (same reference system)")
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    meta = metadata(args.metadata)
    for name in (COHORTS if args.cohort == "both" else (args.cohort,)):
        summarize(name, COHORTS[name], meta, args.output_dir, args.intervals)

if __name__ == "__main__":
    main()
