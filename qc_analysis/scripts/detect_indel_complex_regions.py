#!/usr/bin/env python3
"""Detect indel-centered complex regions among residual native-coordinate HETs.

This step runs after local NUMT filtering/seed expansion.  It does not require a
pre-existing AF-coherent HET cluster.  Instead, each raw-VFC indel is used as a
local center and residual strict-HET SNVs are evaluated in a configurable
circular window together with low-AF SNV and additional-indel context.

Default REMOVE-level evidence patterns (all require >=3 residual strict HETs in
an indel-centered +/-100-bp window):

* AF_MATCHED: at least one local indel has AF within 0.05 of the median HET AF.
* MULTI_INDEL: at least two local indels occur in the window.
* WATERFALL: at least two 1-10% SNVs, at least five total local SNVs
  (residual HET + low-AF), and SNV AF span >=0.10.

Thresholds are optional config keys under ``indel_complex_region``.  When that
section is absent, conservative defaults are used and the raw VCF directory is
inherited from ``pre_liftover_variant_qc.input_vcf_dir``.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median

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


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def add_fields(fields: list[str], extras: list[str]) -> list[str]:
    out = list(fields)
    for field in extras:
        if field not in out:
            out.append(field)
    return out


def fnum(value):
    try:
        out = float(value)
        return out if math.isfinite(out) else None
    except (TypeError, ValueError):
        return None


def inum(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def yes(value) -> bool:
    return str(value).strip().upper() in {"YES", "TRUE", "T", "1"}


def circular_distance(a: int, b: int, length: int) -> int:
    d = abs(a - b)
    return min(d, length - d)


def source_key(row: dict) -> tuple[str, str, str, str, str]:
    return (
        str(row.get("sample", "")),
        str(row.get("source_chrom", "")),
        str(row.get("source_pos", "")),
        str(row.get("source_ref", "")),
        str(row.get("source_alt", "")),
    )


def variant_identity(row: dict) -> tuple:
    return (
        str(row.get("sample", "")),
        str(row.get("pos", "")),
        str(row.get("ref", "")),
        str(row.get("alt", "")),
        str(row.get("af", "")),
        str(row.get("variant_type", "")),
    )


def find_sample_vcf(vcf_dir: Path, sample: str) -> Path | None:
    exact = [vcf_dir / f"{sample}.vcf.gz", vcf_dir / f"{sample}.vcf"]
    for path in exact:
        if path.is_file():
            return path

    matches = []
    for path in vcf_dir.glob(f"{sample}*"):
        name = path.name
        if name.endswith(".vcf") or name.endswith(".vcf.gz"):
            matches.append(path)
    if not matches:
        return None
    matches.sort(key=lambda p: (len(p.name), p.name))
    return matches[0]


def parse_info(info: str) -> dict[str, str]:
    out = {}
    for item in str(info).split(";"):
        if not item or item == ".":
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            out[key] = value
    return out


def parse_format(keys_text: str, values_text: str) -> dict[str, str]:
    keys = str(keys_text).split(":")
    values = str(values_text).split(":")
    return {key: values[i] for i, key in enumerate(keys) if i < len(values)}


def parse_float_list(value: str | None) -> list[float | None]:
    if value in {None, "", "."}:
        return []
    return [fnum(item) for item in str(value).split(",")]


def read_raw_context(
    path: Path,
    sample: str,
    dp_min: int = 100,
    pass_only: bool = True,
    low_snv_af_min: float = 0.01,
    low_snv_af_max: float = 0.10,
    indel_af_min: float = 0.01,
    indel_af_max: float = 0.95,
) -> list[dict]:
    """Read raw native-coordinate LOW_AF_SNV and INDEL context from one VCF."""
    opener = gzip.open if path.name.endswith(".gz") else open
    out: list[dict] = []
    with opener(path, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos_text, _, ref, alt_text, _, filt, info_text = fields[:8]
            pos = inum(pos_text)
            if pos is None or not alt_text or alt_text == ".":
                continue
            if pass_only and filt not in {"PASS", "."}:
                continue

            fmt = parse_format(fields[8], fields[9]) if len(fields) >= 10 else {}
            info = parse_info(info_text)
            dp = fnum(fmt.get("DP"))
            if dp is None:
                dp = fnum(info.get("DP"))
            if dp is not None and dp <= dp_min:
                continue

            af_values = parse_float_list(fmt.get("AF"))
            if not af_values or all(value is None for value in af_values):
                af_values = parse_float_list(info.get("AF"))
            ad_values = parse_float_list(fmt.get("AD"))
            alts = alt_text.split(",")

            for j, alt in enumerate(alts):
                if not alt or alt == "*" or (alt.startswith("<") and alt.endswith(">")):
                    continue
                af = af_values[j] if j < len(af_values) else None
                if af is None and len(ad_values) >= j + 2:
                    finite_ad = [value for value in ad_values if value is not None]
                    total_ad = sum(finite_ad) if finite_ad else 0.0
                    alt_ad = ad_values[j + 1]
                    if total_ad > 0 and alt_ad is not None:
                        af = alt_ad / total_ad
                if af is None:
                    continue

                is_snv = (
                    len(ref) == 1
                    and len(alt) == 1
                    and ref.upper() in {"A", "C", "G", "T"}
                    and alt.upper() in {"A", "C", "G", "T"}
                )
                is_indel = len(ref) != len(alt)
                variant_type = None
                if is_snv and low_snv_af_min <= af <= low_snv_af_max:
                    variant_type = "LOW_AF_SNV"
                elif is_indel and indel_af_min <= af < indel_af_max:
                    variant_type = "INDEL"
                if variant_type is None:
                    continue
                out.append(
                    {
                        "sample": sample,
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "af": af,
                        "dp": dp,
                        "filter": filt,
                        "variant_type": variant_type,
                    }
                )
    return out


def reset_previous_indel_calls(variants: list[dict]) -> None:
    """Make the step safely rerunnable after a fresh or repeated NUMT expansion."""
    for row in variants:
        if str(row.get("indel_complex_region", "")).upper() != "YES":
            continue
        row["filter_action"] = row.get("pre_indel_filter_action", "KEEP") or "KEEP"
        row["filter_reason"] = row.get("pre_indel_filter_reason", "")
        for field in (
            "indel_complex_region_id",
            "indel_complex_pattern",
            "indel_complex_nearest_indel_distance_bp",
            "indel_complex_best_delta_af",
        ):
            row[field] = ""
        row["indel_complex_region"] = "NO"


def detect_candidate_windows(
    residuals: list[dict],
    indels: list[dict],
    low_snvs: list[dict],
    mt_length: int = 16569,
    radius_bp: int = 100,
    min_residual_het: int = 3,
    af_match_max_delta: float = 0.05,
    multiple_indel_min: int = 2,
    waterfall_min_low_snv: int = 2,
    waterfall_min_total_snv: int = 5,
    waterfall_min_snv_af_span: float = 0.10,
) -> list[dict]:
    """Return indel-centered windows meeting at least one artifact pattern."""
    candidates = []
    for center in indels:
        center_pos = int(center["pos"])
        local_hets = [
            row for row in residuals
            if circular_distance(int(row["pos"]), center_pos, mt_length) <= radius_bp
        ]
        if len(local_hets) < min_residual_het:
            continue
        local_indels = [
            row for row in indels
            if circular_distance(int(row["pos"]), center_pos, mt_length) <= radius_bp
        ]
        local_low = [
            row for row in low_snvs
            if circular_distance(int(row["pos"]), center_pos, mt_length) <= radius_bp
        ]

        het_afs = [float(row["af"]) for row in local_hets]
        het_median = median(het_afs)
        best_delta = min(abs(float(row["af"]) - het_median) for row in local_indels)
        snv_afs = het_afs + [float(row["af"]) for row in local_low]
        snv_af_span = max(snv_afs) - min(snv_afs) if snv_afs else 0.0

        patterns = []
        if best_delta <= af_match_max_delta + 1e-12:
            patterns.append("AF_MATCHED")
        if len(local_indels) >= multiple_indel_min:
            patterns.append("MULTI_INDEL")
        if (
            len(local_low) >= waterfall_min_low_snv
            and len(local_hets) + len(local_low) >= waterfall_min_total_snv
            and snv_af_span >= waterfall_min_snv_af_span - 1e-12
        ):
            patterns.append("WATERFALL")
        if not patterns:
            continue

        candidates.append(
            {
                "sample": str(center.get("sample", "")),
                "center": center,
                "center_pos": center_pos,
                "hets": local_hets,
                "indels": local_indels,
                "low_snvs": local_low,
                "patterns": patterns,
                "het_median_af": het_median,
                "het_af_span": max(het_afs) - min(het_afs),
                "snv_context_af_span": snv_af_span,
                "best_indel_delta_af": best_delta,
            }
        )
    return candidates


def merge_candidate_windows(candidates: list[dict]) -> list[list[dict]]:
    """Merge candidate windows that share at least one residual HET variant."""
    n = len(candidates)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    het_sets = [set(int(row["variant_index"]) for row in item["hets"]) for item in candidates]
    for i in range(n):
        for j in range(i + 1, n):
            if candidates[i]["sample"] != candidates[j]["sample"]:
                continue
            if het_sets[i].intersection(het_sets[j]):
                union(i, j)

    groups = defaultdict(list)
    for i, item in enumerate(candidates):
        groups[find(i)].append(item)
    return list(groups.values())


def build_regions(candidate_groups: list[list[dict]], mt_length: int, radius_bp: int) -> list[dict]:
    regions = []
    per_sample_counter = Counter()
    for group in candidate_groups:
        sample = group[0]["sample"]
        per_sample_counter[sample] += 1
        region_id = f"{sample}_ICR{per_sample_counter[sample]:03d}"

        het_by_idx = {}
        indel_by_key = {}
        low_by_key = {}
        patterns = set()
        trigger_positions = set()
        for item in group:
            patterns.update(item["patterns"])
            trigger_positions.add(int(item["center_pos"]))
            for row in item["hets"]:
                het_by_idx[int(row["variant_index"])] = row
            for row in item["indels"]:
                indel_by_key[variant_identity(row)] = row
            for row in item["low_snvs"]:
                low_by_key[variant_identity(row)] = row

        hets = list(het_by_idx.values())
        indels = list(indel_by_key.values())
        low_snvs = list(low_by_key.values())
        het_afs = [float(row["af"]) for row in hets]
        low_afs = [float(row["af"]) for row in low_snvs]
        med = median(het_afs)
        best_delta = min(abs(float(row["af"]) - med) for row in indels) if indels else None
        nearest_dist = min(
            circular_distance(int(h["pos"]), int(ind["pos"]), mt_length)
            for h in hets for ind in indels
        ) if hets and indels else None
        snv_afs = het_afs + low_afs

        regions.append(
            {
                "region_id": region_id,
                "sample": sample,
                "patterns": sorted(patterns),
                "trigger_positions": sorted(trigger_positions),
                "hets": hets,
                "indels": indels,
                "low_snvs": low_snvs,
                "n_residual_het": len(hets),
                "n_preexisting_clustered_het": sum(yes(row.get("clustered")) for row in hets),
                "n_unclustered_het": sum(not yes(row.get("clustered")) for row in hets),
                "n_indels": len(indels),
                "n_low_af_snv": len(low_snvs),
                "het_median_af": med,
                "het_af_min": min(het_afs),
                "het_af_max": max(het_afs),
                "het_af_span": max(het_afs) - min(het_afs),
                "snv_context_af_span": max(snv_afs) - min(snv_afs) if snv_afs else 0.0,
                "best_indel_delta_af": best_delta,
                "nearest_het_indel_distance_bp": nearest_dist,
                "radius_bp": radius_bp,
            }
        )
    return regions


def annotate_variants(variants: list[dict], regions: list[dict], action: str, mt_length: int) -> list[dict]:
    by_index = {}
    for region in regions:
        for het in region["hets"]:
            by_index[int(het["variant_index"])] = region

    annotated = []
    for i, row in enumerate(variants):
        region = by_index.get(i)
        if region is None:
            if not row.get("indel_complex_region"):
                row["indel_complex_region"] = "NO"
            continue

        if not row.get("pre_indel_filter_action"):
            row["pre_indel_filter_action"] = row.get("filter_action", "KEEP") or "KEEP"
            row["pre_indel_filter_reason"] = row.get("filter_reason", "")
        pos = inum(row.get("source_pos"))
        af = fnum(row.get("source_af"))
        nearest = None
        best_delta = None
        if pos is not None:
            nearest = min(circular_distance(pos, int(ind["pos"]), mt_length) for ind in region["indels"])
        if af is not None:
            best_delta = min(abs(af - float(ind["af"])) for ind in region["indels"])

        row["indel_complex_region"] = "YES"
        row["indel_complex_region_id"] = region["region_id"]
        row["indel_complex_pattern"] = "+".join(region["patterns"])
        row["indel_complex_nearest_indel_distance_bp"] = "" if nearest is None else str(nearest)
        row["indel_complex_best_delta_af"] = "" if best_delta is None else f"{best_delta:.6g}"
        if action == "REMOVE":
            row["filter_action"] = "REMOVE"
            row["filter_reason"] = "INDEL_COMPLEX_REGION"
        elif action == "FLAG" and str(row.get("filter_action", "")).upper() == "KEEP":
            row["filter_action"] = "FLAG"
            row["filter_reason"] = "INDEL_COMPLEX_REGION"
        annotated.append(dict(row))
    return annotated


def update_samples(samples: list[dict], variants: list[dict], regions: list[dict]) -> list[dict]:
    region_n = Counter(region["sample"] for region in regions)
    indel_var_n = Counter()
    total_remove = Counter()
    remain = Counter()
    for row in variants:
        sample = str(row.get("sample", ""))
        if str(row.get("indel_complex_region", "")).upper() == "YES":
            indel_var_n[sample] += 1
        if str(row.get("filter_action", "")).upper() == "REMOVE":
            total_remove[sample] += 1
        else:
            remain[sample] += 1

    out = []
    for raw in samples:
        row = dict(raw)
        sample = str(row.get("sample", ""))
        n_het = inum(row.get("n_het")) or 0
        row["n_indel_complex_regions"] = region_n[sample]
        row["n_indel_complex_variants"] = indel_var_n[sample]
        row["n_heteroplasmy_variants_to_remove"] = total_remove[sample]
        row["n_het_after_local_artifact_filter"] = remain[sample]
        row["fraction_het_removed_all"] = f"{total_remove[sample] / n_het:.6g}" if n_het else "0"
        out.append(row)
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()

    cfg = read_simple_yaml(args.config)
    sec = cfg.get("indel_complex_region", {})
    local_sec = cfg.get("local_heteroplasmy_qc", {})
    pre_sec = cfg.get("pre_liftover_variant_qc", {})
    if sec.get("enabled", True) is False:
        print("[indel_complex_region] disabled", file=sys.stderr)
        return 0

    local_out = resolve(local_sec.get("output_dir", "results/qc/local_heteroplasmy_qc"))
    report_dir = local_out / "reports"
    variant_path = report_dir / "local_heteroplasmy_variant_detail.tsv"
    sample_path = report_dir / "local_heteroplasmy_sample_summary.tsv"
    variants = read_tsv(variant_path)
    samples = read_tsv(sample_path)
    if not variants:
        raise RuntimeError(f"No variant rows found: {variant_path}")

    reset_previous_indel_calls(variants)

    vcf_dir_value = sec.get("vcf_dir") or pre_sec.get("input_vcf_dir") or "results/qc/collected_variant_calling_results/collected_vcf"
    vcf_dir = resolve(vcf_dir_value)
    if not vcf_dir.is_dir():
        raise RuntimeError(f"Raw VCF directory does not exist: {vcf_dir}")

    mt_length = int(sec.get("mt_length", local_sec.get("mt_length", 16569)))
    radius_bp = int(sec.get("radius_bp", 100))
    dp_min = int(sec.get("dp_min", local_sec.get("dp_min", 100)))
    pass_only = bool(sec.get("pass_only", True))
    low_min = float(sec.get("low_snv_af_min", 0.01))
    low_max = float(sec.get("low_snv_af_max", 0.10))
    indel_min = float(sec.get("indel_af_min", 0.01))
    indel_max = float(sec.get("indel_af_max", 0.95))
    min_het = int(sec.get("min_residual_het", 3))
    af_delta = float(sec.get("af_match_max_delta", 0.05))
    multi_min = int(sec.get("multiple_indel_min", 2))
    water_low = int(sec.get("waterfall_min_low_snv", 2))
    water_total = int(sec.get("waterfall_min_total_snv", 5))
    water_span = float(sec.get("waterfall_min_snv_af_span", 0.10))
    action = str(sec.get("action", "REMOVE")).strip().upper()
    if action not in {"REMOVE", "FLAG"}:
        raise ValueError("indel_complex_region.action must be REMOVE or FLAG")

    residual_by_sample = defaultdict(list)
    for i, row in enumerate(variants):
        if str(row.get("filter_action", "")).upper() == "REMOVE":
            continue
        pos, af = inum(row.get("source_pos")), fnum(row.get("source_af"))
        if pos is None or af is None:
            continue
        residual_by_sample[str(row.get("sample", ""))].append(
            {
                "variant_index": i,
                "sample": str(row.get("sample", "")),
                "species": str(row.get("species", "")),
                "pos": pos,
                "af": af,
                "clustered": row.get("clustered", "NO"),
                "cluster_id": row.get("cluster_id", ""),
            }
        )

    all_candidates = []
    missing_vcf = []
    raw_context_counts = Counter()
    for sample, residuals in sorted(residual_by_sample.items()):
        if len(residuals) < min_het:
            continue
        vcf = find_sample_vcf(vcf_dir, sample)
        if vcf is None:
            missing_vcf.append(sample)
            continue
        context = read_raw_context(
            vcf,
            sample,
            dp_min=dp_min,
            pass_only=pass_only,
            low_snv_af_min=low_min,
            low_snv_af_max=low_max,
            indel_af_min=indel_min,
            indel_af_max=indel_max,
        )
        indels = [row for row in context if row["variant_type"] == "INDEL"]
        low_snvs = [row for row in context if row["variant_type"] == "LOW_AF_SNV"]
        raw_context_counts["indels"] += len(indels)
        raw_context_counts["low_snvs"] += len(low_snvs)
        if not indels:
            continue
        all_candidates.extend(
            detect_candidate_windows(
                residuals,
                indels,
                low_snvs,
                mt_length=mt_length,
                radius_bp=radius_bp,
                min_residual_het=min_het,
                af_match_max_delta=af_delta,
                multiple_indel_min=multi_min,
                waterfall_min_low_snv=water_low,
                waterfall_min_total_snv=water_total,
                waterfall_min_snv_af_span=water_span,
            )
        )

    regions = build_regions(merge_candidate_windows(all_candidates), mt_length, radius_bp)
    annotated = annotate_variants(variants, regions, action, mt_length)
    samples = update_samples(samples, variants, regions)

    variant_extra = [
        "pre_indel_filter_action",
        "pre_indel_filter_reason",
        "indel_complex_region",
        "indel_complex_region_id",
        "indel_complex_pattern",
        "indel_complex_nearest_indel_distance_bp",
        "indel_complex_best_delta_af",
    ]
    variant_fields = add_fields(list(variants[0]), variant_extra)
    sample_fields = add_fields(
        list(samples[0]) if samples else [],
        [
            "n_indel_complex_regions",
            "n_indel_complex_variants",
            "n_heteroplasmy_variants_to_remove",
            "n_het_after_local_artifact_filter",
            "fraction_het_removed_all",
        ],
    )
    write_tsv(variant_path, variants, variant_fields)
    if samples:
        write_tsv(sample_path, samples, sample_fields)

    region_rows = []
    for region in regions:
        species = ""
        if region["hets"]:
            idx = int(region["hets"][0]["variant_index"])
            species = str(variants[idx].get("species", ""))
        region_rows.append(
            {
                "region_id": region["region_id"],
                "sample": region["sample"],
                "species": species,
                "patterns": "+".join(region["patterns"]),
                "trigger_indel_positions": ",".join(str(x) for x in region["trigger_positions"]),
                "radius_bp": region["radius_bp"],
                "n_residual_het": region["n_residual_het"],
                "n_preexisting_clustered_het": region["n_preexisting_clustered_het"],
                "n_unclustered_het": region["n_unclustered_het"],
                "n_indels": region["n_indels"],
                "n_low_af_snv": region["n_low_af_snv"],
                "het_median_af": f"{region['het_median_af']:.6g}",
                "het_af_min": f"{region['het_af_min']:.6g}",
                "het_af_max": f"{region['het_af_max']:.6g}",
                "het_af_span": f"{region['het_af_span']:.6g}",
                "snv_context_af_span": f"{region['snv_context_af_span']:.6g}",
                "best_indel_delta_af": "" if region["best_indel_delta_af"] is None else f"{region['best_indel_delta_af']:.6g}",
                "nearest_het_indel_distance_bp": "" if region["nearest_het_indel_distance_bp"] is None else region["nearest_het_indel_distance_bp"],
                "action": action,
            }
        )
    region_fields = [
        "region_id", "sample", "species", "patterns", "trigger_indel_positions", "radius_bp",
        "n_residual_het", "n_preexisting_clustered_het", "n_unclustered_het", "n_indels", "n_low_af_snv",
        "het_median_af", "het_af_min", "het_af_max", "het_af_span", "snv_context_af_span",
        "best_indel_delta_af", "nearest_het_indel_distance_bp", "action",
    ]
    write_tsv(report_dir / "indel_complex_regions.tsv", region_rows, region_fields)
    write_tsv(report_dir / "indel_complex_variant_detail.tsv", annotated, variant_fields)

    removals = [dict(row) for row in variants if str(row.get("filter_action", "")).upper() == "REMOVE"]
    # New generic combined report plus legacy filename for the existing final-filter wrapper.
    write_tsv(report_dir / "heteroplasmy_variants_to_remove.tsv", removals, variant_fields)
    write_tsv(report_dir / "numt_variants_to_remove.tsv", removals, variant_fields)

    pattern_counts = Counter()
    for region in regions:
        for pattern in region["patterns"]:
            pattern_counts[pattern] += 1
    indel_removed = sum(str(row.get("indel_complex_region", "")).upper() == "YES" for row in variants)
    summary = [
        {
            "action": action,
            "radius_bp": radius_bp,
            "min_residual_het": min_het,
            "af_match_max_delta": f"{af_delta:.6g}",
            "multiple_indel_min": multi_min,
            "waterfall_min_low_snv": water_low,
            "waterfall_min_total_snv": water_total,
            "waterfall_min_snv_af_span": f"{water_span:.6g}",
            "n_candidate_windows": len(all_candidates),
            "n_complex_regions": len(regions),
            "n_indel_complex_het": indel_removed,
            "n_regions_af_matched": pattern_counts["AF_MATCHED"],
            "n_regions_multi_indel": pattern_counts["MULTI_INDEL"],
            "n_regions_waterfall": pattern_counts["WATERFALL"],
            "n_raw_indels_examined": raw_context_counts["indels"],
            "n_raw_low_snv_examined": raw_context_counts["low_snvs"],
            "n_samples_missing_vcf": len(missing_vcf),
            "n_total_heteroplasmy_removals": len(removals),
        }
    ]
    write_tsv(report_dir / "indel_complex_summary.tsv", summary, list(summary[0]))
    if missing_vcf:
        write_tsv(
            report_dir / "indel_complex_missing_vcf.tsv",
            [{"sample": sample} for sample in missing_vcf],
            ["sample"],
        )

    print(
        f"[indel_complex_region] regions={len(regions)} indel_complex_het={indel_removed} "
        f"patterns(AF_MATCHED={pattern_counts['AF_MATCHED']},MULTI_INDEL={pattern_counts['MULTI_INDEL']},"
        f"WATERFALL={pattern_counts['WATERFALL']}) total_remove={len(removals)} missing_vcf={len(missing_vcf)}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
