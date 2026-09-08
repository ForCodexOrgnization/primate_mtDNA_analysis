#!/usr/bin/env python3
"""Native-coordinate heteroplasmy clustering, NUMT annotation and recurrence QC.

This step is report-only. It detects local AF-coherent HET clusters, annotates
sample- and species-level NUMT evidence, evaluates same-species recurrence, and
emits an explicit source-coordinate variant-removal table for final filtering.
"""
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
from qc_analysis.lib.numt_annotation import (
    best_cluster_overlap,
    best_species_overlap,
    load_sample_numts,
    merge_species_intervals,
)
from qc_analysis.lib.simple_yaml import read_simple_yaml


def resolve(value: str | Path) -> Path:
    path = Path(str(value)).expanduser()
    return path if path.is_absolute() else ROOT / path


def read_tsv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def read_metadata(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    """Read either a headered metadata table or the project's headerless 2-col file.

    The production config/sample_ref_file.tsv is normally:
        SAMPLE<TAB>Species_name
    with no header. In that case the species name is also used as reference_key,
    because all samples in a species are in the same native chrM coordinate system.
    Headered files remain supported and may provide an explicit reference key.
    """
    species: dict[str, str] = {}
    reference: dict[str, str] = {}
    if not path.is_file():
        raise FileNotFoundError(f"sample metadata not found: {path}")

    raw_lines = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line and not line.startswith("#"):
                raw_lines.append(line)
    if not raw_lines:
        raise RuntimeError(f"sample metadata is empty: {path}")

    first = raw_lines[0].split()
    header_tokens = {x.lower() for x in first}
    has_header = bool(header_tokens.intersection({"sample", "species", "reference_key", "mt_reference", "chrm_reference"}))

    if not has_header:
        for line_no, line in enumerate(raw_lines, 1):
            fields = line.split()
            if len(fields) < 2:
                raise RuntimeError(f"invalid headerless sample metadata at {path}:{line_no}: {line}")
            sample, sp = fields[0], fields[1]
            species[sample] = sp
            reference[sample] = sp
        return species, reference

    delimiter = "\t" if "\t" in raw_lines[0] else None
    if delimiter == "\t":
        reader = csv.DictReader(raw_lines, delimiter="\t")
        rows = list(reader)
    else:
        names = first
        rows = [dict(zip(names, line.split())) for line in raw_lines[1:]]

    for row in rows:
        sample = row.get("sample") or row.get("Sample") or ""
        if not sample:
            continue
        sp = row.get("species") or row.get("Species") or ""
        if not sp:
            continue
        species[sample] = sp
        ref = ""
        for name in ("reference_key", "mt_reference", "chrM_reference", "chrm_reference", "species_fasta", "reference_fasta"):
            if row.get(name):
                ref = str(row[name])
                break
        reference[sample] = ref or sp
    return species, reference


def normalize_patterns(patterns, defaults: list[str]) -> list[str]:
    if patterns is None:
        return defaults
    if isinstance(patterns, str):
        parsed = [item.strip() for item in patterns.split(",") if item.strip()]
        return parsed or defaults
    return [str(x) for x in patterns] or defaults


def find_sample_file(directory: Path, patterns: list[str], sample: str) -> Path:
    for pattern in patterns:
        candidate = directory / pattern.format(sample=sample)
        if candidate.is_file():
            return candidate
    return directory / patterns[0].format(sample=sample)


def directory_has_sample_file(directory: Path, patterns: list[str], samples: list[str]) -> bool:
    if not directory.is_dir():
        return False
    for sample in samples[:100]:
        if any((directory / pattern.format(sample=sample)).is_file() for pattern in patterns):
            return True
    return False


def select_numt_directory(configured, patterns: list[str], samples: list[str], fallbacks: list[Path]) -> Path:
    candidates: list[Path] = []
    if configured:
        candidates.append(resolve(configured))
    candidates.extend(fallbacks)
    seen = set()
    for candidate in candidates:
        key = str(candidate)
        if key in seen:
            continue
        seen.add(key)
        if directory_has_sample_file(candidate, patterns, samples):
            return candidate
    # Return the first candidate for a useful error message.
    return candidates[0] if candidates else ROOT / "data/numt_besthit"


def eligible_het(row: dict, sec: dict) -> bool:
    af = number(row.get("source_af"))
    dp = number(row.get("source_dp"))
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


def cluster_class(any_numt: bool, recurrent: bool) -> str:
    if any_numt and recurrent:
        return "NUMT_RECURRENT"
    if any_numt:
        return "NUMT_ONLY"
    if recurrent:
        return "RECURRENT_ONLY"
    return "UNRESOLVED"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()

    cfg = read_simple_yaml(args.config)
    sec = cfg.get("local_heteroplasmy_qc", {})
    if sec.get("enabled", True) is False:
        return 0

    source_report = resolve(sec.get("source_report", "results/qc/pre_liftover_variant_qc/reports/source_variant_qc.tsv"))
    output_dir = resolve(sec.get("output_dir", "results/qc/local_heteroplasmy_qc"))
    report_dir = output_dir / "reports"
    report_dir.mkdir(parents=True, exist_ok=True)

    metadata_path = resolve(sec.get("sample_ref_file", "config/sample_ref_file.tsv"))
    sample_to_species, sample_to_reference = read_metadata(metadata_path)

    rows = read_tsv(source_report)
    eligible = [dict(row) for row in rows if eligible_het(row, sec)]
    for row in eligible:
        row["source_pos"] = int(row["source_pos"])
        row["source_af"] = float(row["source_af"])
        row["source_dp"] = float(row["source_dp"])

    eligible_samples = sorted({row["sample"] for row in eligible})
    missing_metadata = [sample for sample in eligible_samples if not sample_to_species.get(sample)]
    if eligible_samples and len(missing_metadata) == len(eligible_samples):
        raise RuntimeError(
            f"metadata parsing failed: none of {len(eligible_samples)} eligible samples were assigned a species from {metadata_path}"
        )
    if missing_metadata:
        print(
            f"[local_heteroplasmy_qc] WARNING: {len(missing_metadata)}/{len(eligible_samples)} eligible samples lack species metadata",
            file=sys.stderr,
        )

    for sample in eligible_samples:
        sample_to_species.setdefault(sample, "")
        sample_to_reference.setdefault(sample, sample_to_species.get(sample, ""))

    cluster_cfg = ClusterConfig(
        mt_length=int(sec.get("mt_length", 16569)),
        window_bp=int(sec.get("window_bp", 250)),
        af_span_max=float(sec.get("af_span_max", 0.06)),
        min_seed_variants=int(sec.get("min_seed_variants", 3)),
        simulations=int(sec.get("simulations", 2000)),
        empirical_p_max=float(sec.get("empirical_p_max", 0.01)),
        random_seed=int(sec.get("random_seed", 20260904)),
    )

    by_sample: dict[str, list[dict]] = defaultdict(list)
    for row in eligible:
        by_sample[row["sample"]].append(row)
    for values in by_sample.values():
        values.sort(key=lambda x: (x["source_pos"], x.get("source_ref", ""), x.get("source_alt", "")))

    besthit_patterns = normalize_patterns(
        sec.get("numt_besthit_pattern"),
        [
            "{sample}.numt_vs_chrM.besthit.tsv",
            "{sample}.numt_besthit.tsv",
            "{sample}.besthit.tsv",
            "{sample}.tsv",
        ],
    )
    highconf_patterns = normalize_patterns(
        sec.get("numt_highconf_pattern"),
        ["{sample}.highconf_numt.bed"],
    )
    fallback_numt_dir = ROOT / "data/numt_besthit"
    besthit_dir = select_numt_directory(
        sec.get("numt_besthit_dir"), besthit_patterns, eligible_samples, [fallback_numt_dir]
    )
    highconf_dir = select_numt_directory(
        sec.get("numt_highconf_dir"), highconf_patterns, eligible_samples, [besthit_dir, fallback_numt_dir]
    )

    found_besthit = sum(
        1 for sample in eligible_samples
        if any((besthit_dir / pattern.format(sample=sample)).is_file() for pattern in besthit_patterns)
    )
    found_highconf = sum(
        1 for sample in eligible_samples
        if any((highconf_dir / pattern.format(sample=sample)).is_file() for pattern in highconf_patterns)
    )
    if eligible_samples and found_besthit == 0:
        raise RuntimeError(
            "no NUMT besthit files were found for eligible samples; "
            f"directory={besthit_dir}; patterns={','.join(besthit_patterns)}"
        )

    print(
        f"[local_heteroplasmy_qc] metadata={metadata_path} species_assigned={len(eligible_samples)-len(missing_metadata)}/{len(eligible_samples)}",
        file=sys.stderr,
    )
    print(
        f"[local_heteroplasmy_qc] numt_besthit_dir={besthit_dir} besthit_files={found_besthit}/{len(eligible_samples)}",
        file=sys.stderr,
    )
    print(
        f"[local_heteroplasmy_qc] numt_highconf_dir={highconf_dir} highconf_files={found_highconf}/{len(eligible_samples)}",
        file=sys.stderr,
    )

    sample_numts: dict[str, list] = {}
    all_numts = []
    for sample in eligible_samples:
        intervals = load_sample_numts(
            sample=sample,
            species=sample_to_species.get(sample, ""),
            reference_key=sample_to_reference.get(sample, sample_to_species.get(sample, "")),
            besthit_path=find_sample_file(besthit_dir, besthit_patterns, sample),
            highconf_path=find_sample_file(highconf_dir, highconf_patterns, sample),
        )
        sample_numts[sample] = intervals
        all_numts.extend(intervals)

    if found_besthit > 0 and not all_numts:
        raise RuntimeError(
            f"NUMT besthit files were found in {besthit_dir}, but zero intervals were parsed; check besthit file format"
        )

    species_intervals = merge_species_intervals(all_numts, int(sec.get("species_numt_merge_gap_bp", 0)))
    species_rows = []
    for i, interval in enumerate(species_intervals, 1):
        support = sorted(interval["support_samples"])
        loci = sorted(interval["nuclear_loci"])
        species_rows.append({
            "species_interval_id": f"SNUMT{i}",
            "species": interval["species"],
            "reference_key": interval["reference_key"],
            "tier": interval["tier"],
            "chrm_start": interval["chrm_start"],
            "chrm_end": interval["chrm_end"],
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
    min_overlap_n = int(sec.get("numt_min_overlap_variants", 2))
    min_overlap_fraction = float(sec.get("numt_min_overlap_fraction", 0.50))

    cluster_rows: list[dict] = []
    variant_rows: list[dict] = []
    sample_rows: list[dict] = []
    removal_rows: list[dict] = []
    numt_sample_rows: list[dict] = []
    cohort_class_counts = Counter()
    cohort_class_variants = Counter()

    for sample, sample_rows_in in sorted(by_sample.items()):
        species = sample_to_species.get(sample, "")
        reference_key = sample_to_reference.get(sample, species)
        positions = [int(row["source_pos"]) for row in sample_rows_in]
        afs = [float(row["source_af"]) for row in sample_rows_in]
        assessed = len(sample_rows_in) >= cluster_cfg.min_seed_variants
        observed = max_af_coherent_count(positions, afs, cluster_cfg) if assessed else 0
        null = permutation_null(afs, cluster_cfg, zlib.crc32(sample.encode("utf-8"))) if assessed else []
        pvalue = empirical_p(observed, null)
        threshold = critical_count(null, cluster_cfg.empirical_p_max, cluster_cfg.min_seed_variants) if assessed else cluster_cfg.min_seed_variants
        detected = bool(pvalue is not None and pvalue <= cluster_cfg.empirical_p_max)
        seeds = independent_seeds(positions, afs, threshold, cluster_cfg) if detected else []
        expanded = expand_clusters(seeds, positions, afs, cluster_cfg) if seeds else []

        variant_annotation = {
            i: {"cluster_id": "", "cluster_class": "", "filter_action": "KEEP", "filter_reason": ""}
            for i in range(len(sample_rows_in))
        }
        sample_counts = Counter()
        sample_removed = 0

        for cluster_number, indices in enumerate(expanded, 1):
            cluster_id = f"{sample}_C{cluster_number}"
            seed_n = len(seeds[cluster_number - 1]) if cluster_number - 1 < len(seeds) else len(indices)
            members = [sample_rows_in[i] for i in indices]
            summary = summarize_cluster(indices, sample_rows_in, cluster_cfg)

            sample_hit = best_cluster_overlap(members, sample_numts.get(sample, []))
            sample_numt_ok = bool(
                sample_hit
                and sample_hit["overlap_n"] >= min_overlap_n
                and sample_hit["overlap_fraction"] >= min_overlap_fraction
            )

            species_hit = None
            species_numt_ok = False
            if not sample_numt_ok and species:
                species_hit = best_species_overlap(members, sample, species, reference_key, species_intervals)
                species_numt_ok = bool(
                    species_hit
                    and species_hit["overlap_n"] >= min_overlap_n
                    and species_hit["overlap_fraction"] >= min_overlap_fraction
                )

            recurrence = assess_recurrence(
                members,
                sample,
                species,
                recurrence_index,
                min_shared_variants=int(sec.get("recurrence_min_shared_variants", 2)),
                min_shared_fraction=float(sec.get("recurrence_min_shared_fraction", 0.50)),
                max_median_af_difference=float(sec.get("recurrence_max_median_af_difference", 0.05)),
            )
            recurrent = bool(recurrence["recurrent"])
            any_numt = sample_numt_ok or species_numt_ok
            classification = cluster_class(any_numt, recurrent)

            if sample_numt_ok:
                action = str(sec.get("sample_numt_action", "REMOVE")).upper()
                reason = "SAMPLE_NUMT_OVERLAP"
                numt_scope = "SAMPLE"
            elif species_numt_ok and recurrent:
                action = str(sec.get("species_numt_recurrent_action", "REMOVE")).upper()
                reason = "SPECIES_NUMT_AND_RECURRENCE"
                numt_scope = "SPECIES"
            elif species_numt_ok:
                action = str(sec.get("species_numt_only_action", "FLAG")).upper()
                reason = "SPECIES_NUMT_OVERLAP"
                numt_scope = "SPECIES"
            elif recurrent:
                action = str(sec.get("recurrence_only_action", "FLAG")).upper()
                reason = "SAME_SPECIES_RECURRENCE"
                numt_scope = "NONE"
            else:
                action = str(sec.get("unresolved_action", "KEEP")).upper()
                reason = "UNRESOLVED_LOCAL_CLUSTER"
                numt_scope = "NONE"

            hit = sample_hit if sample_numt_ok else species_hit if species_numt_ok else None
            interval = hit["interval"] if hit else None
            if sample_numt_ok:
                numt_tier = interval.tier
                numt_start, numt_end = interval.chrm_start, interval.chrm_end
                species_support_n = 0
                species_support_samples = ""
            elif species_numt_ok:
                numt_tier = interval["tier"]
                numt_start, numt_end = interval["chrm_start"], interval["chrm_end"]
                species_support_samples = ",".join(species_hit["other_support_samples"])
                species_support_n = len(species_hit["other_support_samples"])
            else:
                numt_tier = "NONE"
                numt_start = numt_end = "NA"
                species_support_n = 0
                species_support_samples = ""

            cluster_row = {
                "sample": sample,
                "species": species,
                "reference_key": reference_key,
                "cluster_id": cluster_id,
                "cluster_start": summary["cluster_start"],
                "cluster_end": summary["cluster_end"],
                "cluster_span_bp": summary["cluster_span_bp"],
                "wraps_origin": summary["wraps_origin"],
                "seed_n_variants": seed_n,
                "expanded_n_variants": summary["n_variants"],
                "median_af": f"{summary['median_af']:.6g}",
                "af_min": f"{summary['af_min']:.6g}",
                "af_max": f"{summary['af_max']:.6g}",
                "af_span": f"{summary['af_span']:.6g}",
                "observed_max": observed,
                "critical_count": threshold,
                "empirical_p": "NA" if pvalue is None else f"{pvalue:.6g}",
                "sample_numt_evidence": "YES" if sample_numt_ok else "NO",
                "species_numt_evidence": "YES" if species_numt_ok else "NO",
                "any_numt_evidence": "YES" if any_numt else "NO",
                "numt_scope": numt_scope,
                "numt_tier": numt_tier,
                "numt_chrm_start": numt_start,
                "numt_chrm_end": numt_end,
                "numt_overlap_n": hit["overlap_n"] if hit else 0,
                "numt_overlap_fraction": f"{hit['overlap_fraction']:.6g}" if hit else "0",
                "species_numt_support_n": species_support_n,
                "species_numt_support_samples": species_support_samples,
                "recurrence": "YES" if recurrent else "NO",
                "best_recurrent_sample": recurrence["best_recurrent_sample"],
                "shared_variants": recurrence["shared_variants"],
                "shared_fraction": f"{recurrence['shared_fraction']:.6g}",
                "median_af_difference": recurrence["median_af_difference"] if recurrence["median_af_difference"] == "NA" else f"{recurrence['median_af_difference']:.6g}",
                "cluster_class": classification,
                "filter_action": action,
                "filter_reason": reason,
            }
            cluster_rows.append(cluster_row)
            cohort_class_counts[classification] += 1
            cohort_class_variants[classification] += len(indices)
            sample_counts[classification] += 1

            for i in indices:
                variant_annotation[i] = {
                    "cluster_id": cluster_id,
                    "cluster_class": classification,
                    "filter_action": action,
                    "filter_reason": reason,
                    "numt_scope": numt_scope,
                    "numt_tier": numt_tier,
                    "recurrence": "YES" if recurrent else "NO",
                }
                if action == "REMOVE":
                    sample_removed += 1
                    member = sample_rows_in[i]
                    removal_rows.append({
                        "sample": sample,
                        "species": species,
                        "source_chrom": member.get("source_chrom", ""),
                        "source_pos": member["source_pos"],
                        "source_ref": member.get("source_ref", ""),
                        "source_alt": member.get("source_alt", ""),
                        "source_af": f"{float(member['source_af']):.6g}",
                        "source_dp": f"{float(member['source_dp']):.6g}",
                        "cluster_id": cluster_id,
                        "cluster_class": classification,
                        "numt_scope": numt_scope,
                        "numt_tier": numt_tier,
                        "filter_action": action,
                        "filter_reason": reason,
                    })

        for i, member in enumerate(sample_rows_in):
            ann = variant_annotation[i]
            variant_rows.append({
                "sample": sample,
                "species": species,
                "source_chrom": member.get("source_chrom", ""),
                "source_pos": member["source_pos"],
                "source_ref": member.get("source_ref", ""),
                "source_alt": member.get("source_alt", ""),
                "source_af": f"{float(member['source_af']):.6g}",
                "source_dp": f"{float(member['source_dp']):.6g}",
                "clustered": "YES" if ann["cluster_id"] else "NO",
                "cluster_id": ann["cluster_id"],
                "cluster_class": ann["cluster_class"],
                "numt_scope": ann.get("numt_scope", "NONE"),
                "numt_tier": ann.get("numt_tier", "NONE"),
                "recurrence": ann.get("recurrence", "NO"),
                "filter_action": ann["filter_action"],
                "filter_reason": ann["filter_reason"],
            })

        n_clusters = len(expanded)
        n_clustered_het = sum(len(cluster) for cluster in expanded)
        n_numt_clusters = sample_counts["NUMT_RECURRENT"] + sample_counts["NUMT_ONLY"]
        sample_summary = {
            "sample": sample,
            "species": species,
            "n_het": len(sample_rows_in),
            "assessed": str(assessed).lower(),
            "observed_max": observed,
            "critical_count": threshold if assessed else "NA",
            "empirical_p": "NA" if pvalue is None else f"{pvalue:.6g}",
            "n_clusters": n_clusters,
            "n_clustered_het": n_clustered_het,
            "n_numt_clusters": n_numt_clusters,
            "n_numt_recurrent_clusters": sample_counts["NUMT_RECURRENT"],
            "n_numt_only_clusters": sample_counts["NUMT_ONLY"],
            "n_recurrent_only_clusters": sample_counts["RECURRENT_ONLY"],
            "n_unresolved_clusters": sample_counts["UNRESOLVED"],
            "n_numt_variants_to_remove": sample_removed,
            "fraction_het_removed": f"{sample_removed / len(sample_rows_in):.6g}" if sample_rows_in else "0",
            "numt_sample": "YES" if n_numt_clusters > 0 else "NO",
            "sample_filter_action": "KEEP",
        }
        sample_rows.append(sample_summary)
        if n_numt_clusters > 0:
            numt_sample_rows.append(sample_summary.copy())

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

    cohort_rows = []
    total_clusters = sum(cohort_class_counts.values())
    total_clustered = sum(cohort_class_variants.values())
    for classification in ("NUMT_RECURRENT", "NUMT_ONLY", "RECURRENT_ONLY", "UNRESOLVED"):
        n_clusters = cohort_class_counts[classification]
        n_variants = cohort_class_variants[classification]
        cohort_rows.append({
            "cluster_class": classification,
            "n_clusters": n_clusters,
            "fraction_clusters": f"{n_clusters / total_clusters:.6g}" if total_clusters else "0",
            "n_clustered_het": n_variants,
            "fraction_clustered_het": f"{n_variants / total_clustered:.6g}" if total_clustered else "0",
        })
    write_tsv(
        report_dir / "heteroplasmy_analysis_summary.tsv",
        cohort_rows,
        ["cluster_class", "n_clusters", "fraction_clusters", "n_clustered_het", "fraction_clustered_het"],
    )

    diagnostic_rows = [{
        "metadata_path": str(metadata_path),
        "eligible_samples": len(eligible_samples),
        "samples_with_species": len(eligible_samples) - len(missing_metadata),
        "samples_missing_species": len(missing_metadata),
        "numt_besthit_dir": str(besthit_dir),
        "samples_with_besthit_file": found_besthit,
        "numt_highconf_dir": str(highconf_dir),
        "samples_with_highconf_file": found_highconf,
        "parsed_numt_intervals": len(all_numts),
        "species_numt_intervals": len(species_rows),
        "detected_clusters": total_clusters,
    }]
    write_tsv(
        report_dir / "heteroplasmy_input_diagnostics.tsv",
        diagnostic_rows,
        [
            "metadata_path", "eligible_samples", "samples_with_species", "samples_missing_species",
            "numt_besthit_dir", "samples_with_besthit_file", "numt_highconf_dir",
            "samples_with_highconf_file", "parsed_numt_intervals", "species_numt_intervals", "detected_clusters",
        ],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
