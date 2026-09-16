#!/usr/bin/env python3
"""Apply local indel-associated artifact filtering to residual heteroplasmy.

This production step runs after ``detect_indel_complex_regions.py`` and before
later residual-artifact sensitivity / recurrent-block analyses.

Three complementary rules are applied to the same pre-step residual HET universe:

0) DIRECT_INDEL_PROXIMITY
   - any residual strict-HET SNV whose position is <=10 bp from a raw PASS indel;
   - no minimum HET count and no AF-matching requirement.
   -> REMOVE that HET only.

   This deliberately targets the visually obvious case where an SNV and indel
   essentially overlap.  The default 10-bp distance is intentionally much more
   stringent than the group-level windows below.

1) COHERENT_INDEL_CLUSTER
   - >=3 residual HETs with minimum circular span <=250 bp;
   - group AF span <=0.10;
   - nearest indel <=100 bp;
   - |group median AF - nearby indel AF| <=0.05.
   -> REMOVE the whole group.

2) COMPLEX_INDEL_REGION
   - >=3 residual HETs with minimum circular span <=250 bp;
   - group AF span >0.10;
   - nearest indel <=100 bp;
   - and either >=2 indels within 250 bp of the group OR >=2 low-AF SNVs
     within 250 bp of the group.
   -> REMOVE the whole group.

A narrow high-AF band near a low-AF indel is therefore protected unless the indel
AF actually matches the HET group, while a single HET that nearly coincides with
an indel can still be removed by the much stricter DIRECT_INDEL_PROXIMITY rule.
All rules use each sample's native mitochondrial length from the raw VCF header.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from qc_analysis.lib.simple_yaml import read_simple_yaml
from qc_analysis.scripts.detect_indel_complex_regions import (
    choose_native_length,
    circular_distance,
    find_sample_vcf,
    read_raw_context,
)


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
        return float(value)
    except (TypeError, ValueError):
        return None


def inum(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def yes(value) -> bool:
    return str(value).strip().upper() in {"YES", "TRUE", "T", "1"}


def reset_previous_calls(variants: list[dict]) -> None:
    """Restore pre-step decisions from prior implementations for deterministic reruns."""
    for row in variants:
        previously_applied = (
            yes(row.get("dense_indel_filtered"))
            or yes(row.get("het_group_indel_filtered"))
            or yes(row.get("direct_indel_proximity"))
        )
        if not previously_applied:
            continue
        row["filter_action"] = (
            row.get("pre_het_group_indel_filter_action")
            or row.get("pre_dense_indel_filter_action")
            or "KEEP"
        )
        row["filter_reason"] = (
            row.get("pre_het_group_indel_filter_reason")
            if row.get("pre_het_group_indel_filter_reason") not in {None, ""}
            else row.get("pre_dense_indel_filter_reason", "")
        )
        for field in (
            "dense_indel_filtered",
            "dense_indel_overlap",
            "multi_indel_dense",
            "dense_indel_rules",
            "dense_indel_trigger_positions",
            "dense_indel_min_distance_bp",
            "het_group_indel_filtered",
            "het_group_indel_rule",
            "het_group_indel_group_ids",
            "het_group_span_bp",
            "het_group_af_span",
            "het_group_median_af",
            "het_group_nearest_indel_distance_bp",
            "het_group_best_indel_delta_af",
            "het_group_n_indels_250bp",
            "het_group_n_low_snv_250bp",
            "direct_indel_proximity",
            "direct_indel_nearest_distance_bp",
            "direct_indel_positions",
            "direct_indel_best_delta_af",
        ):
            row[field] = ""


def minimum_circular_span(positions: list[int], mt_length: int | None) -> int:
    if len(positions) <= 1:
        return 0
    positions = sorted(positions)
    linear_span = positions[-1] - positions[0]
    if mt_length is None or mt_length <= 0:
        return linear_span
    gaps = [positions[i + 1] - positions[i] for i in range(len(positions) - 1)]
    gaps.append(positions[0] + mt_length - positions[-1])
    return mt_length - max(gaps)


def local_het_groups(
    residuals: list[dict],
    mt_length: int | None,
    max_span_bp: int = 250,
    min_het: int = 3,
) -> list[dict]:
    """Enumerate maximal unique HET sets with minimum circular span <= max_span_bp."""
    if len(residuals) < min_het:
        return []

    ordered = sorted(residuals, key=lambda row: (int(row["pos"]), int(row["variant_index"])))
    n = len(ordered)
    by_index = {int(row["variant_index"]): row for row in residuals}
    candidate_sets = []

    if mt_length is None or mt_length <= 0:
        for i in range(n):
            members = []
            start = int(ordered[i]["pos"])
            for j in range(i, n):
                if int(ordered[j]["pos"]) - start > max_span_bp:
                    break
                members.append(ordered[j])
            if len(members) >= min_het:
                candidate_sets.append(members)
    else:
        doubled = ordered + [dict(row, pos=int(row["pos"]) + mt_length) for row in ordered]
        for i in range(n):
            start = int(doubled[i]["pos"])
            member_indices = []
            for j in range(i, i + n):
                if int(doubled[j]["pos"]) - start > max_span_bp:
                    break
                idx = int(doubled[j]["variant_index"])
                if idx not in member_indices:
                    member_indices.append(idx)
            if len(member_indices) >= min_het:
                candidate_sets.append([by_index[idx] for idx in member_indices])

    unique = {}
    for members in candidate_sets:
        key = frozenset(int(row["variant_index"]) for row in members)
        unique[key] = members

    ordered_candidates = sorted(
        unique.items(),
        key=lambda item: (-len(item[0]), min(int(row["pos"]) for row in item[1])),
    )
    accepted_sets = []
    accepted_groups = []
    for member_set, members in ordered_candidates:
        if any(member_set < existing for existing in accepted_sets):
            continue
        accepted_sets.append(member_set)
        accepted_groups.append(members)

    groups = []
    for members in accepted_groups:
        positions = [int(row["pos"]) for row in members]
        afs = [float(row["af"]) for row in members]
        groups.append(
            {
                "members": members,
                "positions": positions,
                "n_het": len(members),
                "span_bp": minimum_circular_span(positions, mt_length),
                "median_af": median(afs),
                "af_min": min(afs),
                "af_max": max(afs),
                "af_span": max(afs) - min(afs),
            }
        )
    return groups


def distance_variant_to_group(pos: int, group_positions: list[int], mt_length: int | None) -> int:
    return min(circular_distance(pos, gpos, mt_length) for gpos in group_positions)


def classify_het_group(
    group: dict,
    indels: list[dict],
    low_snvs: list[dict],
    mt_length: int | None,
    nearby_indel_bp: int = 100,
    context_bp: int = 250,
    coherent_max_af_span: float = 0.10,
    coherent_max_delta_af: float = 0.05,
    complex_min_indels: int = 2,
    complex_min_low_snvs: int = 2,
) -> dict:
    group_positions = group["positions"]

    indel_context = []
    for indel in indels:
        distance = distance_variant_to_group(int(indel["pos"]), group_positions, mt_length)
        if distance <= context_bp:
            item = dict(indel)
            item["distance_to_group"] = distance
            item["delta_group_median_af"] = abs(float(indel["af"]) - float(group["median_af"]))
            indel_context.append(item)

    low_context = []
    for snv in low_snvs:
        distance = distance_variant_to_group(int(snv["pos"]), group_positions, mt_length)
        if distance <= context_bp:
            item = dict(snv)
            item["distance_to_group"] = distance
            low_context.append(item)

    nearby_indels = [row for row in indel_context if int(row["distance_to_group"]) <= nearby_indel_bp]
    nearest_indel_distance = min((int(row["distance_to_group"]) for row in indel_context), default=None)
    best_delta_nearby = min((float(row["delta_group_median_af"]) for row in nearby_indels), default=None)

    coherent = (
        group["af_span"] <= coherent_max_af_span + 1e-12
        and bool(nearby_indels)
        and best_delta_nearby is not None
        and best_delta_nearby <= coherent_max_delta_af + 1e-12
    )

    complex_region = (
        group["af_span"] > coherent_max_af_span + 1e-12
        and bool(nearby_indels)
        and (
            len(indel_context) >= complex_min_indels
            or len(low_context) >= complex_min_low_snvs
        )
    )

    rule = "COHERENT_INDEL_CLUSTER" if coherent else ("COMPLEX_INDEL_REGION" if complex_region else "")
    return {
        "rule": rule,
        "nearest_indel_distance_bp": nearest_indel_distance,
        "best_indel_delta_af": best_delta_nearby,
        "n_indels_250bp": len(indel_context),
        "n_low_snv_250bp": len(low_context),
        "nearby_indel_positions": sorted({int(row["pos"]) for row in nearby_indels}),
        "context_indel_positions": sorted({int(row["pos"]) for row in indel_context}),
    }


def update_sample_summary(samples: list[dict], variants: list[dict]) -> list[dict]:
    direct_n = Counter()
    coherent_n = Counter()
    complex_n = Counter()
    filtered_n = Counter()
    total_remove = Counter()
    residual_n = Counter()
    for row in variants:
        sample = str(row.get("sample", ""))
        if yes(row.get("direct_indel_proximity")):
            direct_n[sample] += 1
        rule = str(row.get("het_group_indel_rule", ""))
        if "COHERENT_INDEL_CLUSTER" in rule:
            coherent_n[sample] += 1
        if "COMPLEX_INDEL_REGION" in rule:
            complex_n[sample] += 1
        if yes(row.get("dense_indel_filtered")):
            filtered_n[sample] += 1
        if str(row.get("filter_action", "")).upper() == "REMOVE":
            total_remove[sample] += 1
        else:
            residual_n[sample] += 1

    out = []
    for raw in samples:
        row = dict(raw)
        sample = str(row.get("sample", ""))
        n_het = inum(row.get("n_het")) or 0
        row["n_direct_indel_proximity_variants"] = direct_n[sample]
        row["n_coherent_indel_cluster_variants"] = coherent_n[sample]
        row["n_complex_indel_region_variants"] = complex_n[sample]
        row["n_het_group_indel_filtered_variants"] = coherent_n[sample] + complex_n[sample]
        row["n_dense_indel_filtered_variants"] = filtered_n[sample]
        row["n_heteroplasmy_variants_to_remove"] = total_remove[sample]
        row["n_het_after_local_artifact_filter"] = residual_n[sample]
        row["fraction_het_removed_all"] = f"{total_remove[sample] / n_het:.6g}" if n_het else "0"
        out.append(row)
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()

    cfg = read_simple_yaml(args.config)
    sec = cfg.get("dense_indel_overlap", {})
    local_sec = cfg.get("local_heteroplasmy_qc", {})
    pre_sec = cfg.get("pre_liftover_variant_qc", {})
    if sec.get("enabled", True) is False:
        print("[local_indel_artifact_filter] disabled", file=sys.stderr)
        return 0

    report_dir = resolve(local_sec.get("output_dir", "results/qc/local_heteroplasmy_qc")) / "reports"
    variant_path = report_dir / "local_heteroplasmy_variant_detail.tsv"
    sample_path = report_dir / "local_heteroplasmy_sample_summary.tsv"
    variants = read_tsv(variant_path)
    samples = read_tsv(sample_path)
    if not variants:
        raise RuntimeError(f"No variant rows found: {variant_path}")

    reset_previous_calls(variants)

    vcf_dir_value = sec.get("vcf_dir") or pre_sec.get("input_vcf_dir") or "results/qc/collected_variant_calling_results/collected_vcf"
    vcf_dir = resolve(vcf_dir_value)
    if not vcf_dir.is_dir():
        raise RuntimeError(f"Raw VCF directory does not exist: {vcf_dir}")

    dp_min = int(sec.get("dp_min", local_sec.get("dp_min", 100)))
    pass_only = bool(sec.get("pass_only", True))
    low_min = float(sec.get("low_snv_af_min", 0.01))
    low_max = float(sec.get("low_snv_af_max", 0.10))
    indel_min = float(sec.get("indel_af_min", 0.01))
    indel_max = float(sec.get("indel_af_max", 0.95))

    direct_max_distance = int(sec.get("direct_indel_max_distance_bp", 10))
    group_span = int(sec.get("het_group_max_span_bp", 250))
    group_min_het = int(sec.get("het_group_min_residual_het", 3))
    nearby_indel_bp = int(sec.get("het_group_nearby_indel_bp", 100))
    context_bp = int(sec.get("het_group_context_bp", 250))
    coherent_af_span = float(sec.get("coherent_max_af_span", 0.10))
    coherent_delta = float(sec.get("coherent_max_indel_delta_af", 0.05))
    complex_min_indels = int(sec.get("complex_min_indels", 2))
    complex_min_low_snvs = int(sec.get("complex_min_low_snvs", 2))

    residual_by_sample = defaultdict(list)
    for i, row in enumerate(variants):
        if str(row.get("filter_action", "")).upper() == "REMOVE":
            continue
        pos = inum(row.get("source_pos"))
        af = fnum(row.get("source_af"))
        if pos is None or af is None:
            continue
        residual_by_sample[str(row.get("sample", ""))].append({
            "variant_index": i,
            "sample": str(row.get("sample", "")),
            "species": str(row.get("species", "")),
            "chrom": str(row.get("source_chrom", "")),
            "pos": pos,
            "af": af,
        })

    group_report_rows = []
    direct_report_rows = []
    group_evidence_by_variant = defaultdict(list)
    direct_evidence_by_variant = {}
    missing_vcf = []
    group_counter = Counter()
    rule_groups = Counter()

    for sample, residuals in sorted(residual_by_sample.items()):
        vcf = find_sample_vcf(vcf_dir, sample)
        if vcf is None:
            missing_vcf.append(sample)
            continue

        context, contig_lengths = read_raw_context(
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
        if not indels:
            continue

        native_length, _, _ = choose_native_length(contig_lengths, residuals, context)

        # ----------------------------------------------------
        # Rule 0: direct SNV-indel proximity.
        # Evaluate every residual HET independently, including isolated HETs.
        # ----------------------------------------------------
        for het in residuals:
            distances = [
                (circular_distance(int(het["pos"]), int(indel["pos"]), native_length), indel)
                for indel in indels
            ]
            nearest_distance = min(distance for distance, _ in distances)
            if nearest_distance > direct_max_distance:
                continue
            nearest_indels = [indel for distance, indel in distances if distance == nearest_distance]
            best_delta = min(abs(float(het["af"]) - float(indel["af"])) for indel in nearest_indels)
            direct_evidence_by_variant[int(het["variant_index"])] = {
                "distance": nearest_distance,
                "indel_positions": sorted({int(indel["pos"]) for indel in nearest_indels}),
                "best_delta_af": best_delta,
            }
            direct_report_rows.append({
                "sample": sample,
                "species": het.get("species", ""),
                "het_pos": het["pos"],
                "het_af": f"{float(het['af']):.6g}",
                "nearest_indel_distance_bp": nearest_distance,
                "nearest_indel_positions": ",".join(str(int(indel["pos"])) for indel in nearest_indels),
                "best_indel_delta_af": f"{best_delta:.6g}",
                "action": "REMOVE",
            })

        # ----------------------------------------------------
        # Group rules are evaluated on the same pre-step residual universe.
        # Direct hits are not removed first, so they cannot break up groups.
        # ----------------------------------------------------
        if len(residuals) < group_min_het:
            continue
        groups = local_het_groups(residuals, native_length, group_span, group_min_het)

        for group in groups:
            classification = classify_het_group(
                group,
                indels,
                low_snvs,
                native_length,
                nearby_indel_bp,
                context_bp,
                coherent_af_span,
                coherent_delta,
                complex_min_indels,
                complex_min_low_snvs,
            )
            rule = classification["rule"]
            if not rule:
                continue

            group_counter[sample] += 1
            group_id = f"{sample}_HIG{group_counter[sample]:03d}"
            rule_groups[rule] += 1
            for member in group["members"]:
                group_evidence_by_variant[int(member["variant_index"])].append({
                    "group_id": group_id,
                    "rule": rule,
                    "span_bp": group["span_bp"],
                    "af_span": group["af_span"],
                    "median_af": group["median_af"],
                    **classification,
                })

            group_report_rows.append({
                "group_id": group_id,
                "sample": sample,
                "species": group["members"][0].get("species", ""),
                "rule": rule,
                "n_het": group["n_het"],
                "group_span_bp": group["span_bp"],
                "het_positions": ",".join(str(x) for x in sorted(group["positions"])),
                "het_median_af": f"{group['median_af']:.6g}",
                "het_af_min": f"{group['af_min']:.6g}",
                "het_af_max": f"{group['af_max']:.6g}",
                "het_af_span": f"{group['af_span']:.6g}",
                "nearest_indel_distance_bp": classification["nearest_indel_distance_bp"] if classification["nearest_indel_distance_bp"] is not None else "",
                "best_indel_delta_af": f"{classification['best_indel_delta_af']:.6g}" if classification["best_indel_delta_af"] is not None else "",
                "n_indels_250bp": classification["n_indels_250bp"],
                "n_low_snv_250bp": classification["n_low_snv_250bp"],
                "nearby_indel_positions": ",".join(str(x) for x in classification["nearby_indel_positions"]),
                "context_indel_positions": ",".join(str(x) for x in classification["context_indel_positions"]),
            })

    all_filtered_indices = set(group_evidence_by_variant) | set(direct_evidence_by_variant)
    for idx in all_filtered_indices:
        row = variants[idx]
        if not row.get("pre_het_group_indel_filter_action"):
            row["pre_het_group_indel_filter_action"] = row.get("filter_action", "KEEP") or "KEEP"
            row["pre_het_group_indel_filter_reason"] = row.get("filter_reason", "")

        group_evidence = group_evidence_by_variant.get(idx, [])
        direct_evidence = direct_evidence_by_variant.get(idx)
        group_rules = sorted({item["rule"] for item in group_evidence})

        if direct_evidence is not None:
            row["direct_indel_proximity"] = "YES"
            row["direct_indel_nearest_distance_bp"] = str(direct_evidence["distance"])
            row["direct_indel_positions"] = ",".join(str(x) for x in direct_evidence["indel_positions"])
            row["direct_indel_best_delta_af"] = f"{direct_evidence['best_delta_af']:.6g}"
        else:
            row["direct_indel_proximity"] = "NO"

        if group_evidence:
            row["het_group_indel_filtered"] = "YES"
            row["het_group_indel_rule"] = "+".join(group_rules)
            row["het_group_indel_group_ids"] = ";".join(sorted({item["group_id"] for item in group_evidence}))
            row["het_group_span_bp"] = str(min(int(item["span_bp"]) for item in group_evidence))
            row["het_group_af_span"] = f"{min(float(item['af_span']) for item in group_evidence):.6g}"
            row["het_group_median_af"] = f"{median(float(item['median_af']) for item in group_evidence):.6g}"
            distance_values = [item["nearest_indel_distance_bp"] for item in group_evidence if item["nearest_indel_distance_bp"] is not None]
            delta_values = [item["best_indel_delta_af"] for item in group_evidence if item["best_indel_delta_af"] is not None]
            row["het_group_nearest_indel_distance_bp"] = str(min(distance_values)) if distance_values else ""
            row["het_group_best_indel_delta_af"] = f"{min(delta_values):.6g}" if delta_values else ""
            row["het_group_n_indels_250bp"] = str(max(int(item["n_indels_250bp"]) for item in group_evidence))
            row["het_group_n_low_snv_250bp"] = str(max(int(item["n_low_snv_250bp"]) for item in group_evidence))
        else:
            row["het_group_indel_filtered"] = "NO"

        all_rules = (["DIRECT_INDEL_PROXIMITY"] if direct_evidence is not None else []) + group_rules
        row["filter_action"] = "REMOVE"
        if direct_evidence is not None:
            row["filter_reason"] = "DIRECT_INDEL_PROXIMITY"
        elif "COHERENT_INDEL_CLUSTER" in group_rules:
            row["filter_reason"] = "COHERENT_INDEL_CLUSTER"
        else:
            row["filter_reason"] = "COMPLEX_INDEL_REGION"

        # Legacy fields remain populated so older plotting/summary scripts still work.
        row["dense_indel_filtered"] = "YES"
        row["dense_indel_rules"] = "+".join(all_rules)
        distance_candidates = []
        if direct_evidence is not None:
            distance_candidates.append(int(direct_evidence["distance"]))
        distance_candidates.extend(
            int(item["nearest_indel_distance_bp"])
            for item in group_evidence
            if item["nearest_indel_distance_bp"] is not None
        )
        row["dense_indel_min_distance_bp"] = str(min(distance_candidates)) if distance_candidates else ""

    variant_extra = [
        "pre_het_group_indel_filter_action",
        "pre_het_group_indel_filter_reason",
        "direct_indel_proximity",
        "direct_indel_nearest_distance_bp",
        "direct_indel_positions",
        "direct_indel_best_delta_af",
        "het_group_indel_filtered",
        "het_group_indel_rule",
        "het_group_indel_group_ids",
        "het_group_span_bp",
        "het_group_af_span",
        "het_group_median_af",
        "het_group_nearest_indel_distance_bp",
        "het_group_best_indel_delta_af",
        "het_group_n_indels_250bp",
        "het_group_n_low_snv_250bp",
        "dense_indel_filtered",
        "dense_indel_rules",
        "dense_indel_min_distance_bp",
    ]
    variant_fields = add_fields(list(variants[0]), variant_extra)
    write_tsv(variant_path, variants, variant_fields)

    if samples:
        samples = update_sample_summary(samples, variants)
        sample_fields = add_fields(list(samples[0]), [
            "n_direct_indel_proximity_variants",
            "n_coherent_indel_cluster_variants",
            "n_complex_indel_region_variants",
            "n_het_group_indel_filtered_variants",
            "n_dense_indel_filtered_variants",
            "n_heteroplasmy_variants_to_remove",
            "n_het_after_local_artifact_filter",
            "fraction_het_removed_all",
        ])
        write_tsv(sample_path, samples, sample_fields)

    group_detail = [dict(row) for row in variants if yes(row.get("het_group_indel_filtered"))]
    all_detail = [dict(row) for row in variants if yes(row.get("dense_indel_filtered"))]
    direct_detail = [dict(row) for row in variants if yes(row.get("direct_indel_proximity"))]
    write_tsv(report_dir / "het_group_indel_variant_detail.tsv", group_detail, variant_fields)
    write_tsv(report_dir / "direct_indel_proximity_variant_detail.tsv", direct_detail, variant_fields)
    write_tsv(report_dir / "local_indel_artifact_variant_detail.tsv", all_detail, variant_fields)
    write_tsv(report_dir / "dense_indel_overlap_variant_detail.tsv", all_detail, variant_fields)

    group_fields = [
        "group_id", "sample", "species", "rule", "n_het", "group_span_bp", "het_positions",
        "het_median_af", "het_af_min", "het_af_max", "het_af_span",
        "nearest_indel_distance_bp", "best_indel_delta_af", "n_indels_250bp", "n_low_snv_250bp",
        "nearby_indel_positions", "context_indel_positions",
    ]
    write_tsv(report_dir / "het_group_indel_regions.tsv", group_report_rows, group_fields)
    write_tsv(report_dir / "dense_indel_overlap_windows.tsv", group_report_rows, group_fields)
    write_tsv(
        report_dir / "direct_indel_proximity.tsv",
        direct_report_rows,
        [
            "sample", "species", "het_pos", "het_af", "nearest_indel_distance_bp",
            "nearest_indel_positions", "best_indel_delta_af", "action",
        ],
    )

    removals = [dict(row) for row in variants if str(row.get("filter_action", "")).upper() == "REMOVE"]
    write_tsv(report_dir / "heteroplasmy_variants_to_remove.tsv", removals, variant_fields)
    write_tsv(report_dir / "numt_variants_to_remove.tsv", removals, variant_fields)

    group_filtered_indices = [i for i, row in enumerate(variants) if yes(row.get("het_group_indel_filtered"))]
    direct_filtered_indices = [i for i, row in enumerate(variants) if yes(row.get("direct_indel_proximity"))]
    all_filtered_indices = [i for i, row in enumerate(variants) if yes(row.get("dense_indel_filtered"))]
    rule_variant_counts = Counter()
    for i in group_filtered_indices:
        for rule in str(variants[i].get("het_group_indel_rule", "")).split("+"):
            if rule:
                rule_variant_counts[rule] += 1

    summary = [{
        "direct_indel_max_distance_bp": direct_max_distance,
        "n_direct_indel_proximity_variants": len(direct_filtered_indices),
        "n_direct_indel_proximity_samples": len({str(variants[i].get("sample", "")) for i in direct_filtered_indices}),
        "het_group_max_span_bp": group_span,
        "het_group_min_residual_het": group_min_het,
        "nearby_indel_bp": nearby_indel_bp,
        "context_bp": context_bp,
        "coherent_max_af_span": f"{coherent_af_span:.6g}",
        "coherent_max_indel_delta_af": f"{coherent_delta:.6g}",
        "complex_min_indels": complex_min_indels,
        "complex_min_low_snvs": complex_min_low_snvs,
        "n_coherent_indel_groups": rule_groups["COHERENT_INDEL_CLUSTER"],
        "n_complex_indel_groups": rule_groups["COMPLEX_INDEL_REGION"],
        "n_coherent_indel_variant_memberships": rule_variant_counts["COHERENT_INDEL_CLUSTER"],
        "n_complex_indel_variant_memberships": rule_variant_counts["COMPLEX_INDEL_REGION"],
        "n_unique_group_filtered_variants": len(group_filtered_indices),
        "n_unique_local_indel_artifact_variants": len(all_filtered_indices),
        "n_samples_local_indel_artifact_filtered": len({str(variants[i].get("sample", "")) for i in all_filtered_indices}),
        "n_samples_missing_vcf": len(missing_vcf),
        "n_total_heteroplasmy_removals": len(removals),
    }]
    write_tsv(report_dir / "het_group_indel_summary.tsv", summary, list(summary[0]))
    write_tsv(report_dir / "dense_indel_overlap_summary.tsv", summary, list(summary[0]))

    if missing_vcf:
        write_tsv(report_dir / "dense_indel_overlap_missing_vcf.tsv", [{"sample": s} for s in missing_vcf], ["sample"])

    print(
        "[local_indel_artifact_filter] "
        f"direct={len(direct_filtered_indices)} "
        f"coherent_groups={rule_groups['COHERENT_INDEL_CLUSTER']} "
        f"complex_groups={rule_groups['COMPLEX_INDEL_REGION']} "
        f"unique_local_indel_removed={len(all_filtered_indices)} "
        f"samples={len({str(variants[i].get('sample', '')) for i in all_filtered_indices})} "
        f"total_remove={len(removals)}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
