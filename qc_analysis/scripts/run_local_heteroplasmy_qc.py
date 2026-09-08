#!/usr/bin/env python3
"""Native-coordinate heteroplasmy clustering, NUMT annotation and recurrence QC."""
from __future__ import annotations

import argparse
import csv
import math
import sys
import zlib
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from qc_analysis.lib.heteroplasmy_clusters import (
    ClusterConfig,
    critical_count,
    empirical_p,
    expand_clusters,
    independent_seeds,
    max_af_coherent_count,
    permutation_null,
    summarize_cluster,
)
from qc_analysis.lib.heteroplasmy_recurrence import assess_recurrence, build_species_sample_index
from qc_analysis.lib.numt_annotation import best_cluster_overlap, best_species_overlap, load_sample_numts, merge_species_intervals
from qc_analysis.lib.simple_yaml import read_simple_yaml


def resolve(value: str | Path) -> Path:
    p = Path(str(value)).expanduser()
    return p if p.is_absolute() else ROOT / p


def read_tsv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as h:
        return list(csv.DictReader(h, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def number(value):
    try:
        x = float(value)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def read_metadata(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    """Support headered metadata and the production headerless SAMPLE SPECIES file."""
    if not path.is_file():
        raise FileNotFoundError(f"sample metadata not found: {path}")
    lines = []
    with path.open(encoding="utf-8") as h:
        for line in h:
            line = line.strip()
            if line and not line.startswith("#"):
                lines.append(line)
    if not lines:
        raise RuntimeError(f"sample metadata is empty: {path}")

    first = lines[0].split()
    header = {x.lower() for x in first}
    has_header = bool(header.intersection({"sample", "species", "reference_key", "mt_reference", "chrm_reference"}))
    species: dict[str, str] = {}
    reference: dict[str, str] = {}

    if not has_header:
        for line_no, line in enumerate(lines, 1):
            f = line.split()
            if len(f) < 2:
                raise RuntimeError(f"invalid sample_ref_file line {line_no}: {line}")
            sample, sp = f[0], f[1]
            species[sample] = sp
            reference[sample] = sp
        return species, reference

    if "\t" in lines[0]:
        rows = list(csv.DictReader(lines, delimiter="\t"))
    else:
        names = first
        rows = [dict(zip(names, line.split())) for line in lines[1:]]
    for row in rows:
        sample = row.get("sample") or row.get("Sample") or ""
        sp = row.get("species") or row.get("Species") or ""
        if not sample or not sp:
            continue
        species[sample] = sp
        ref = ""
        for k in ("reference_key", "mt_reference", "chrM_reference", "chrm_reference", "species_fasta", "reference_fasta"):
            if row.get(k):
                ref = str(row[k])
                break
        reference[sample] = ref or sp
    return species, reference


def read_sample_qc(path: Path) -> tuple[set[str], int]:
    """Return the sample-QC PASS set; qc_status must be present and non-empty."""
    if not path.is_file():
        raise FileNotFoundError(
            f"sample QC report not found: {path}. Run sample_variant_filtering before local_heteroplasmy_qc."
        )
    rows = read_tsv(path)
    if not rows:
        raise RuntimeError(f"sample QC report is empty: {path}")
    if "sample" not in rows[0] or "qc_status" not in rows[0]:
        raise RuntimeError(f"sample QC report must contain sample and qc_status columns: {path}")
    passed = {
        str(row.get("sample", "")).strip()
        for row in rows
        if str(row.get("qc_status", "")).strip().upper() == "PASS" and str(row.get("sample", "")).strip()
    }
    if not passed:
        raise RuntimeError(f"sample QC report contains zero PASS samples: {path}")
    return passed, len(rows)


def patterns(value, defaults: list[str]) -> list[str]:
    if isinstance(value, str):
        configured = [x.strip() for x in value.split(",") if x.strip()]
    elif value:
        configured = [str(x) for x in value]
    else:
        configured = []
    out = []
    for p in configured + defaults:
        if p not in out:
            out.append(p)
    return out


def find_sample_file(directory: Path, pats: list[str], sample: str) -> Path:
    for pat in pats:
        p = directory / pat.format(sample=sample)
        if p.is_file():
            return p
    return directory / pats[0].format(sample=sample)


def count_sample_files(directory: Path, pats: list[str], samples: list[str]) -> int:
    if not directory.is_dir():
        return 0
    return sum(any((directory / p.format(sample=s)).is_file() for p in pats) for s in samples)


def choose_numt_dir(configured, pats: list[str], samples: list[str], fallbacks: list[Path]) -> Path:
    candidates = []
    if configured:
        candidates.append(resolve(configured))
    candidates.extend(fallbacks)
    seen = set()
    for d in candidates:
        if str(d) in seen:
            continue
        seen.add(str(d))
        if count_sample_files(d, pats, samples) > 0:
            return d
    return candidates[0]


def eligible_het(row: dict, sec: dict) -> bool:
    af, dp = number(row.get("source_af")), number(row.get("source_dp"))
    return (
        row.get("source_call_class") == "HET"
        and row.get("source_filter") == "PASS"
        and row.get("source_variant_qc") == "PASS"
        and dp is not None
        and dp > float(sec.get("dp_min", 100))
        and af is not None
        and af > float(sec.get("heteroplasmy_af_min", 0.10))
        and af < float(sec.get("heteroplasmy_af_max", 0.95))
    )


def classify(any_numt: bool, recurrent: bool) -> str:
    if any_numt and recurrent:
        return "NUMT_RECURRENT"
    if any_numt:
        return "NUMT_ONLY"
    if recurrent:
        return "RECURRENT_ONLY"
    return "UNRESOLVED"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, required=True)
    args = ap.parse_args()
    cfg = read_simple_yaml(args.config)
    sec = cfg.get("local_heteroplasmy_qc", {})
    if sec.get("enabled", True) is False:
        return 0

    source_report = resolve(sec.get("source_report", "results/qc/pre_liftover_variant_qc/reports/source_variant_qc.tsv"))
    sample_qc_report = resolve(sec.get("sample_qc_report", "results/qc/sample_variant_filtering/reports/sample_qc.tsv"))
    report_dir = resolve(sec.get("output_dir", "results/qc/local_heteroplasmy_qc")) / "reports"
    report_dir.mkdir(parents=True, exist_ok=True)

    metadata_path = resolve(sec.get("sample_ref_file", "config/sample_ref_file.tsv"))
    sample_to_species, sample_to_reference = read_metadata(metadata_path)
    qc_pass, n_sample_qc_rows = read_sample_qc(sample_qc_report)

    raw_rows = read_tsv(source_report)
    source_samples = sorted({str(r.get("sample", "")).strip() for r in raw_rows if str(r.get("sample", "")).strip()})

    # The analysis cohort is defined by sample_variant_filtering PASS, not by the
    # presence of a heteroplasmic call. This keeps cluster detection consistent
    # with the high-quality sample set used downstream.
    analysis_samples = sorted(s for s in qc_pass if sample_to_species.get(s))
    qc_pass_missing_metadata = sorted(s for s in qc_pass if not sample_to_species.get(s))
    if not analysis_samples:
        raise RuntimeError("zero sample-QC PASS samples have species metadata")

    strict_het_all = [dict(r) for r in raw_rows if eligible_het(r, sec)]
    strict_het_samples_before_qc = {r["sample"] for r in strict_het_all}
    eligible = [r for r in strict_het_all if r.get("sample") in qc_pass]
    for r in eligible:
        r["source_pos"] = int(r["source_pos"])
        r["source_af"] = float(r["source_af"])
        r["source_dp"] = float(r["source_dp"])

    eligible_samples = sorted({r["sample"] for r in eligible})
    missing_species = [s for s in eligible_samples if not sample_to_species.get(s)]
    if eligible_samples and len(missing_species) == len(eligible_samples):
        raise RuntimeError(f"metadata parsing failed: 0/{len(eligible_samples)} eligible samples have species in {metadata_path}")

    for s in analysis_samples:
        sample_to_reference.setdefault(s, sample_to_species.get(s, ""))

    by_sample: dict[str, list[dict]] = defaultdict(list)
    for r in eligible:
        if sample_to_species.get(r["sample"]):
            by_sample[r["sample"]].append(r)
    for vals in by_sample.values():
        vals.sort(key=lambda r: (r["source_pos"], r.get("source_ref", ""), r.get("source_alt", "")))

    cfg_cluster = ClusterConfig(
        mt_length=int(sec.get("mt_length", 16569)),
        window_bp=int(sec.get("window_bp", 250)),
        af_span_max=float(sec.get("af_span_max", 0.06)),
        min_seed_variants=int(sec.get("min_seed_variants", 3)),
        simulations=int(sec.get("simulations", 2000)),
        empirical_p_max=float(sec.get("empirical_p_max", 0.01)),
        random_seed=int(sec.get("random_seed", 20260904)),
    )

    besthit_pats = patterns(sec.get("numt_besthit_pattern"), [
        "{sample}.numt_vs_chrM.besthit.tsv",
        "{sample}.numt_besthit.tsv",
        "{sample}.besthit.tsv",
        "{sample}.tsv",
    ])
    highconf_pats = patterns(sec.get("numt_highconf_pattern"), ["{sample}.highconf_numt.bed"])
    repo_numt_dir = ROOT / "data/numt_besthit"

    # Build the species NUMT catalogue from the full QC-PASS cohort, including
    # PASS samples that have zero strict HETs. They can still provide valid
    # species-level NUMT support for another sample.
    besthit_dir = choose_numt_dir(sec.get("numt_besthit_dir"), besthit_pats, analysis_samples, [repo_numt_dir])
    highconf_dir = choose_numt_dir(sec.get("numt_highconf_dir"), highconf_pats, analysis_samples, [besthit_dir, repo_numt_dir])
    n_besthit_files = count_sample_files(besthit_dir, besthit_pats, analysis_samples)
    n_highconf_files = count_sample_files(highconf_dir, highconf_pats, analysis_samples)
    if analysis_samples and n_besthit_files == 0:
        raise RuntimeError(
            f"no NUMT besthit files found; checked {besthit_dir}; expected e.g. SAMPLE.numt_vs_chrM.besthit.tsv"
        )

    print(
        f"[local_heteroplasmy_qc] sample_qc PASS={len(qc_pass)}/{n_sample_qc_rows} path={sample_qc_report}",
        file=sys.stderr,
    )
    print(
        f"[local_heteroplasmy_qc] analysis_cohort={len(analysis_samples)} strict_HET_samples={len(eligible_samples)}",
        file=sys.stderr,
    )
    print(
        f"[local_heteroplasmy_qc] metadata species={len(analysis_samples)}/{len(qc_pass)} path={metadata_path}",
        file=sys.stderr,
    )
    print(f"[local_heteroplasmy_qc] besthit={n_besthit_files}/{len(analysis_samples)} dir={besthit_dir}", file=sys.stderr)
    print(f"[local_heteroplasmy_qc] highconf={n_highconf_files}/{len(analysis_samples)} dir={highconf_dir}", file=sys.stderr)

    sample_numts, all_numts = {}, []
    for sample in analysis_samples:
        ivals = load_sample_numts(
            sample,
            sample_to_species.get(sample, ""),
            sample_to_reference.get(sample, sample_to_species.get(sample, "")),
            find_sample_file(besthit_dir, besthit_pats, sample),
            find_sample_file(highconf_dir, highconf_pats, sample),
        )
        sample_numts[sample] = ivals
        all_numts.extend(ivals)
    if n_besthit_files and not all_numts:
        raise RuntimeError("NUMT besthit files were found but zero intervals parsed; inspect besthit TSV format")

    species_intervals = merge_species_intervals(all_numts, int(sec.get("species_numt_merge_gap_bp", 0)))
    species_rows = []
    for i, iv in enumerate(species_intervals, 1):
        support = sorted(iv["support_samples"])
        loci = sorted(iv["nuclear_loci"])
        species_rows.append({
            "species_interval_id": f"SNUMT{i}",
            "species": iv["species"],
            "reference_key": iv["reference_key"],
            "tier": iv["tier"],
            "chrm_start": iv["chrm_start"],
            "chrm_end": iv["chrm_end"],
            "n_support_samples": len(support),
            "support_samples": ",".join(support),
            "n_nuclear_loci": len(loci),
            "nuclear_loci": ";".join(f"{c}:{s}-{e}" for c, s, e in loci),
        })
    write_tsv(
        report_dir / "species_numt_intervals.tsv",
        species_rows,
        ["species_interval_id", "species", "reference_key", "tier", "chrm_start", "chrm_end", "n_support_samples", "support_samples", "n_nuclear_loci", "nuclear_loci"],
    )

    recurrence_index = build_species_sample_index(eligible, sample_to_species)
    min_ov_n = int(sec.get("numt_min_overlap_variants", 2))
    min_ov_frac = float(sec.get("numt_min_overlap_fraction", 0.50))

    cluster_rows, variant_rows, sample_rows, removal_rows, numt_sample_rows = [], [], [], [], []
    class_counts, class_vars = Counter(), Counter()
    sample_runtime = {}

    for sample in analysis_samples:
        srows = by_sample.get(sample, [])
        species = sample_to_species.get(sample, "")
        refkey = sample_to_reference.get(sample, species)
        pos = [int(r["source_pos"]) for r in srows]
        afs = [float(r["source_af"]) for r in srows]
        assessed = len(srows) >= cfg_cluster.min_seed_variants
        observed = max_af_coherent_count(pos, afs, cfg_cluster) if assessed else 0
        null = permutation_null(afs, cfg_cluster, zlib.crc32(sample.encode())) if assessed else []
        pval = empirical_p(observed, null)
        crit = critical_count(null, cfg_cluster.empirical_p_max, cfg_cluster.min_seed_variants) if assessed else cfg_cluster.min_seed_variants
        detected = pval is not None and pval <= cfg_cluster.empirical_p_max
        seeds = independent_seeds(pos, afs, crit, cfg_cluster) if detected else []
        expanded = expand_clusters(seeds, pos, afs, cfg_cluster) if seeds else []

        ann = {i: {"cluster_id": "", "cluster_class": "", "filter_action": "KEEP", "filter_reason": ""} for i in range(len(srows))}
        scount, removed = Counter(), 0

        for cnum, indices in enumerate(expanded, 1):
            cid = f"{sample}_C{cnum}"
            members = [srows[i] for i in indices]
            sm = summarize_cluster(indices, srows, cfg_cluster)
            seed_n = len(seeds[cnum - 1]) if cnum - 1 < len(seeds) else len(indices)

            sample_hit = best_cluster_overlap(members, sample_numts.get(sample, []))
            sample_ok = bool(sample_hit and sample_hit["overlap_n"] >= min_ov_n and sample_hit["overlap_fraction"] >= min_ov_frac)

            species_hit, species_ok = None, False
            if not sample_ok and species:
                species_hit = best_species_overlap(members, sample, species, refkey, species_intervals)
                species_ok = bool(species_hit and species_hit["overlap_n"] >= min_ov_n and species_hit["overlap_fraction"] >= min_ov_frac)

            rec = assess_recurrence(
                members,
                sample,
                species,
                recurrence_index,
                min_shared_variants=int(sec.get("recurrence_min_shared_variants", 2)),
                min_shared_fraction=float(sec.get("recurrence_min_shared_fraction", 0.50)),
                max_median_af_difference=float(sec.get("recurrence_max_median_af_difference", 0.05)),
            )
            recurrent = bool(rec["recurrent"])
            any_numt = sample_ok or species_ok
            cls = classify(any_numt, recurrent)

            if sample_ok:
                action, reason, scope = str(sec.get("sample_numt_action", "REMOVE")).upper(), "SAMPLE_NUMT_OVERLAP", "SAMPLE"
            elif species_ok and recurrent:
                action, reason, scope = str(sec.get("species_numt_recurrent_action", "REMOVE")).upper(), "SPECIES_NUMT_AND_RECURRENCE", "SPECIES"
            elif species_ok:
                action, reason, scope = str(sec.get("species_numt_only_action", "FLAG")).upper(), "SPECIES_NUMT_OVERLAP", "SPECIES"
            elif recurrent:
                action, reason, scope = str(sec.get("recurrence_only_action", "FLAG")).upper(), "SAME_SPECIES_RECURRENCE", "NONE"
            else:
                action, reason, scope = str(sec.get("unresolved_action", "KEEP")).upper(), "UNRESOLVED_LOCAL_CLUSTER", "NONE"

            hit = sample_hit if sample_ok else species_hit if species_ok else None
            interval = hit["interval"] if hit else None
            if sample_ok:
                tier = interval.tier
                numt_start, numt_end = interval.chrm_start, interval.chrm_end
                support_n, support_samples = 0, ""
            elif species_ok:
                tier = interval["tier"]
                numt_start, numt_end = interval["chrm_start"], interval["chrm_end"]
                support_samples = ",".join(species_hit["other_support_samples"])
                support_n = len(species_hit["other_support_samples"])
            else:
                tier, numt_start, numt_end, support_n, support_samples = "NONE", "NA", "NA", 0, ""

            cluster_rows.append({
                "sample": sample,
                "species": species,
                "reference_key": refkey,
                "cluster_id": cid,
                "cluster_start": sm["cluster_start"],
                "cluster_end": sm["cluster_end"],
                "cluster_span_bp": sm["cluster_span_bp"],
                "wraps_origin": sm["wraps_origin"],
                "seed_n_variants": seed_n,
                "expanded_n_variants": sm["n_variants"],
                "median_af": f"{sm['median_af']:.6g}",
                "af_min": f"{sm['af_min']:.6g}",
                "af_max": f"{sm['af_max']:.6g}",
                "af_span": f"{sm['af_span']:.6g}",
                "observed_max": observed,
                "critical_count": crit,
                "empirical_p": "NA" if pval is None else f"{pval:.6g}",
                "sample_numt_evidence": "YES" if sample_ok else "NO",
                "species_numt_evidence": "YES" if species_ok else "NO",
                "any_numt_evidence": "YES" if any_numt else "NO",
                "numt_scope": scope,
                "numt_tier": tier,
                "numt_chrm_start": numt_start,
                "numt_chrm_end": numt_end,
                "numt_overlap_n": hit["overlap_n"] if hit else 0,
                "numt_overlap_fraction": f"{hit['overlap_fraction']:.6g}" if hit else "0",
                "species_numt_support_n": support_n,
                "species_numt_support_samples": support_samples,
                "recurrence": "YES" if recurrent else "NO",
                "best_recurrent_sample": rec["best_recurrent_sample"],
                "shared_variants": rec["shared_variants"],
                "shared_fraction": f"{rec['shared_fraction']:.6g}",
                "median_af_difference": rec["median_af_difference"] if rec["median_af_difference"] == "NA" else f"{rec['median_af_difference']:.6g}",
                "cluster_class": cls,
                "filter_action": action,
                "filter_reason": reason,
            })
            class_counts[cls] += 1
            class_vars[cls] += len(indices)
            scount[cls] += 1

            for i in indices:
                ann[i] = {
                    "cluster_id": cid,
                    "cluster_class": cls,
                    "filter_action": action,
                    "filter_reason": reason,
                    "numt_scope": scope,
                    "numt_tier": tier,
                    "recurrence": "YES" if recurrent else "NO",
                }
                if action == "REMOVE":
                    removed += 1
                    m = srows[i]
                    removal_rows.append({
                        "sample": sample,
                        "species": species,
                        "source_chrom": m.get("source_chrom", ""),
                        "source_pos": m["source_pos"],
                        "source_ref": m.get("source_ref", ""),
                        "source_alt": m.get("source_alt", ""),
                        "source_af": f"{float(m['source_af']):.6g}",
                        "source_dp": f"{float(m['source_dp']):.6g}",
                        "cluster_id": cid,
                        "cluster_class": cls,
                        "numt_scope": scope,
                        "numt_tier": tier,
                        "filter_action": action,
                        "filter_reason": reason,
                    })

        for i, m in enumerate(srows):
            a = ann[i]
            variant_rows.append({
                "sample": sample,
                "species": species,
                "source_chrom": m.get("source_chrom", ""),
                "source_pos": m["source_pos"],
                "source_ref": m.get("source_ref", ""),
                "source_alt": m.get("source_alt", ""),
                "source_af": f"{float(m['source_af']):.6g}",
                "source_dp": f"{float(m['source_dp']):.6g}",
                "clustered": "YES" if a["cluster_id"] else "NO",
                "cluster_id": a["cluster_id"],
                "cluster_class": a["cluster_class"],
                "numt_scope": a.get("numt_scope", "NONE"),
                "numt_tier": a.get("numt_tier", "NONE"),
                "recurrence": a.get("recurrence", "NO"),
                "filter_action": a["filter_action"],
                "filter_reason": a["filter_reason"],
            })

        n_clusters = len(expanded)
        n_clustered_het = sum(len(c) for c in expanded)
        n_numt_clusters = scount["NUMT_RECURRENT"] + scount["NUMT_ONLY"]
        summary = {
            "sample": sample,
            "species": species,
            "n_het": len(srows),
            "assessed": str(assessed).lower(),
            "observed_max": observed,
            "critical_count": crit if assessed else "NA",
            "empirical_p": "NA" if pval is None else f"{pval:.6g}",
            "n_clusters": n_clusters,
            "n_clustered_het": n_clustered_het,
            "n_numt_clusters": n_numt_clusters,
            "n_numt_recurrent_clusters": scount["NUMT_RECURRENT"],
            "n_numt_only_clusters": scount["NUMT_ONLY"],
            "n_recurrent_only_clusters": scount["RECURRENT_ONLY"],
            "n_unresolved_clusters": scount["UNRESOLVED"],
            "n_numt_variants_to_remove": removed,
            "fraction_het_removed": f"{removed / len(srows):.6g}" if srows else "0",
            "numt_sample": "YES" if n_numt_clusters > 0 else "NO",
            "sample_filter_action": "KEEP",
        }
        sample_rows.append(summary)
        sample_runtime[sample] = summary
        if n_numt_clusters > 0:
            numt_sample_rows.append(summary.copy())

    cluster_fields = [
        "sample", "species", "reference_key", "cluster_id", "cluster_start", "cluster_end", "cluster_span_bp", "wraps_origin",
        "seed_n_variants", "expanded_n_variants", "median_af", "af_min", "af_max", "af_span", "observed_max", "critical_count", "empirical_p",
        "sample_numt_evidence", "species_numt_evidence", "any_numt_evidence", "numt_scope", "numt_tier", "numt_chrm_start", "numt_chrm_end",
        "numt_overlap_n", "numt_overlap_fraction", "species_numt_support_n", "species_numt_support_samples", "recurrence", "best_recurrent_sample",
        "shared_variants", "shared_fraction", "median_af_difference", "cluster_class", "filter_action", "filter_reason",
    ]
    variant_fields = [
        "sample", "species", "source_chrom", "source_pos", "source_ref", "source_alt", "source_af", "source_dp", "clustered", "cluster_id",
        "cluster_class", "numt_scope", "numt_tier", "recurrence", "filter_action", "filter_reason",
    ]
    sample_fields = [
        "sample", "species", "n_het", "assessed", "observed_max", "critical_count", "empirical_p", "n_clusters", "n_clustered_het",
        "n_numt_clusters", "n_numt_recurrent_clusters", "n_numt_only_clusters", "n_recurrent_only_clusters", "n_unresolved_clusters",
        "n_numt_variants_to_remove", "fraction_het_removed", "numt_sample", "sample_filter_action",
    ]
    removal_fields = [
        "sample", "species", "source_chrom", "source_pos", "source_ref", "source_alt", "source_af", "source_dp", "cluster_id", "cluster_class",
        "numt_scope", "numt_tier", "filter_action", "filter_reason",
    ]

    write_tsv(report_dir / "local_heteroplasmy_cluster_summary.tsv", cluster_rows, cluster_fields)
    write_tsv(report_dir / "local_heteroplasmy_variant_detail.tsv", variant_rows, variant_fields)
    write_tsv(report_dir / "local_heteroplasmy_sample_summary.tsv", sample_rows, sample_fields)
    write_tsv(report_dir / "numt_samples.tsv", numt_sample_rows, sample_fields)
    write_tsv(report_dir / "numt_variants_to_remove.tsv", removal_rows, removal_fields)

    total_clusters = sum(class_counts.values())
    total_clustered = sum(class_vars.values())
    cohort_rows = []
    for cls in ("NUMT_RECURRENT", "NUMT_ONLY", "RECURRENT_ONLY", "UNRESOLVED"):
        ncl = class_counts[cls]
        nvar = class_vars[cls]
        cohort_rows.append({
            "cluster_class": cls,
            "n_clusters": ncl,
            "fraction_clusters": f"{ncl / total_clusters:.6g}" if total_clusters else "0",
            "n_clustered_het": nvar,
            "fraction_clustered_het": f"{nvar / total_clustered:.6g}" if total_clustered else "0",
        })
    write_tsv(
        report_dir / "heteroplasmy_analysis_summary.tsv",
        cohort_rows,
        ["cluster_class", "n_clusters", "fraction_clusters", "n_clustered_het", "fraction_clustered_het"],
    )

    diagnostic_rows = [{
        "source_report": str(source_report),
        "n_samples_in_source_report": len(source_samples),
        "sample_qc_report": str(sample_qc_report),
        "n_samples_in_sample_qc_report": n_sample_qc_rows,
        "n_samples_sample_qc_pass": len(qc_pass),
        "n_sample_qc_pass_missing_metadata": len(qc_pass_missing_metadata),
        "n_analysis_samples": len(analysis_samples),
        "n_strict_het_samples_before_sample_qc": len(strict_het_samples_before_qc),
        "n_samples_with_eligible_het": len(eligible_samples),
        "n_samples_assessed_for_clusters": sum(str(r["assessed"]).lower() == "true" for r in sample_rows),
        "n_samples_with_clusters": sum(int(r["n_clusters"]) > 0 for r in sample_rows),
        "metadata_path": str(metadata_path),
        "numt_besthit_dir": str(besthit_dir),
        "samples_with_besthit_file": n_besthit_files,
        "numt_highconf_dir": str(highconf_dir),
        "samples_with_highconf_file": n_highconf_files,
        "parsed_numt_intervals": len(all_numts),
        "species_numt_intervals": len(species_rows),
        "detected_clusters": total_clusters,
    }]
    write_tsv(
        report_dir / "heteroplasmy_input_diagnostics.tsv",
        diagnostic_rows,
        [
            "source_report", "n_samples_in_source_report", "sample_qc_report", "n_samples_in_sample_qc_report",
            "n_samples_sample_qc_pass", "n_sample_qc_pass_missing_metadata", "n_analysis_samples",
            "n_strict_het_samples_before_sample_qc", "n_samples_with_eligible_het", "n_samples_assessed_for_clusters",
            "n_samples_with_clusters", "metadata_path", "numt_besthit_dir", "samples_with_besthit_file",
            "numt_highconf_dir", "samples_with_highconf_file", "parsed_numt_intervals", "species_numt_intervals", "detected_clusters",
        ],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
