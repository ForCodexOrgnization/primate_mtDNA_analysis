#!/usr/bin/env python3
"""Apply direct dense indel-overlap filtering to residual HETs.

This production step runs after ``detect_indel_complex_regions.py`` and before
later residual-artifact sensitivity / recurrent-block analyses.  It targets the
visually obvious local pattern where several residual strict-HET SNVs sit very
close to raw indels, without requiring AF concordance.

Default rules:

A) DENSE_INDEL_OVERLAP
   >=3 residual HETs within +/-50 bp of the same raw indel -> REMOVE those HETs.

B) MULTI_INDEL_DENSE
   >=3 residual HETs within +/-100 bp of an indel-centered window AND >=2 raw
   indels in that same +/-100-bp window -> REMOVE those HETs.

The rules use each sample's native mitochondrial length from the raw VCF header.
They are deliberately local and do not chain-expand.  A HET can satisfy both
rules; the more spatially stringent DENSE_INDEL_OVERLAP reason takes priority
for ``filter_reason`` while both evidence flags are retained.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

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


def reset_previous_dense_calls(variants: list[dict]) -> None:
    """Restore the pre-step decision so repeated full-pipeline runs are stable."""
    for row in variants:
        if not yes(row.get("dense_indel_filtered")):
            continue
        row["filter_action"] = row.get("pre_dense_indel_filter_action", "KEEP") or "KEEP"
        row["filter_reason"] = row.get("pre_dense_indel_filter_reason", "")
        row["dense_indel_filtered"] = "NO"
        for field in (
            "dense_indel_overlap",
            "multi_indel_dense",
            "dense_indel_rules",
            "dense_indel_trigger_positions",
            "dense_indel_min_distance_bp",
        ):
            row[field] = ""


def detect_dense_indel_windows(
    residuals: list[dict],
    indels: list[dict],
    mt_length: int | None,
    dense_radius_bp: int = 50,
    dense_min_het: int = 3,
    multi_radius_bp: int = 100,
    multi_min_het: int = 3,
    multi_min_indels: int = 2,
) -> list[dict]:
    """Return qualifying indel-centered windows without mutating variants."""
    windows = []
    for center in indels:
        center_pos = int(center["pos"])

        dense_hets = [
            row for row in residuals
            if circular_distance(int(row["pos"]), center_pos, mt_length) <= dense_radius_bp
        ]
        multi_hets = [
            row for row in residuals
            if circular_distance(int(row["pos"]), center_pos, mt_length) <= multi_radius_bp
        ]
        local_indels = [
            row for row in indels
            if circular_distance(int(row["pos"]), center_pos, mt_length) <= multi_radius_bp
        ]

        rules = []
        remove_indices: set[int] = set()
        if len(dense_hets) >= dense_min_het:
            rules.append("DENSE_INDEL_OVERLAP")
            remove_indices.update(int(row["variant_index"]) for row in dense_hets)

        if len(multi_hets) >= multi_min_het and len(local_indels) >= multi_min_indels:
            rules.append("MULTI_INDEL_DENSE")
            remove_indices.update(int(row["variant_index"]) for row in multi_hets)

        if not rules:
            continue

        windows.append(
            {
                "sample": str(center.get("sample", "")),
                "center_pos": center_pos,
                "center_af": center.get("af", ""),
                "rules": rules,
                "dense_hets": dense_hets,
                "multi_hets": multi_hets,
                "local_indels": local_indels,
                "remove_indices": remove_indices,
            }
        )
    return windows


def combine_window_evidence(
    windows: list[dict],
    residual_by_index: dict[int, dict],
    mt_length: int | None,
) -> tuple[dict[int, dict], list[dict]]:
    """Collapse overlapping qualifying windows to per-HET evidence and report rows."""
    evidence: dict[int, dict] = {}
    report_rows = []

    for window_i, window in enumerate(windows, start=1):
        rules_text = "+".join(window["rules"])
        removed_positions = sorted(
            int(residual_by_index[idx]["pos"])
            for idx in window["remove_indices"]
            if idx in residual_by_index
        )
        report_rows.append(
            {
                "sample": window["sample"],
                "window_id": f"{window['sample']}_DIW{window_i:03d}",
                "center_indel_pos": window["center_pos"],
                "center_indel_af": window["center_af"],
                "rules": rules_text,
                "n_dense_het_50bp": len(window["dense_hets"]),
                "n_het_100bp": len(window["multi_hets"]),
                "n_indels_100bp": len(window["local_indels"]),
                "n_unique_het_removed": len(window["remove_indices"]),
                "removed_het_positions": ",".join(str(pos) for pos in removed_positions),
            }
        )

        for idx in window["remove_indices"]:
            row = residual_by_index.get(idx)
            if row is None:
                continue
            item = evidence.setdefault(
                idx,
                {
                    "rules": set(),
                    "trigger_positions": set(),
                    "min_distance": None,
                },
            )
            for rule in window["rules"]:
                # Only attribute the rule if this HET actually belongs to that rule's window.
                distance = circular_distance(int(row["pos"]), int(window["center_pos"]), mt_length)
                if rule == "DENSE_INDEL_OVERLAP" and distance <= 50:
                    item["rules"].add(rule)
                elif rule == "MULTI_INDEL_DENSE" and distance <= 100:
                    item["rules"].add(rule)
            item["trigger_positions"].add(int(window["center_pos"]))
            distance = circular_distance(int(row["pos"]), int(window["center_pos"]), mt_length)
            if item["min_distance"] is None or distance < item["min_distance"]:
                item["min_distance"] = distance

    return evidence, report_rows


def update_sample_summary(samples: list[dict], variants: list[dict]) -> list[dict]:
    dense_n = Counter()
    multi_n = Counter()
    filtered_n = Counter()
    total_remove = Counter()
    residual_n = Counter()

    for row in variants:
        sample = str(row.get("sample", ""))
        if yes(row.get("dense_indel_overlap")):
            dense_n[sample] += 1
        if yes(row.get("multi_indel_dense")):
            multi_n[sample] += 1
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
        row["n_dense_indel_overlap_variants"] = dense_n[sample]
        row["n_multi_indel_dense_variants"] = multi_n[sample]
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
        print("[dense_indel_overlap] disabled", file=sys.stderr)
        return 0

    report_dir = resolve(local_sec.get("output_dir", "results/qc/local_heteroplasmy_qc")) / "reports"
    variant_path = report_dir / "local_heteroplasmy_variant_detail.tsv"
    sample_path = report_dir / "local_heteroplasmy_sample_summary.tsv"
    variants = read_tsv(variant_path)
    samples = read_tsv(sample_path)
    if not variants:
        raise RuntimeError(f"No variant rows found: {variant_path}")

    reset_previous_dense_calls(variants)

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

    dense_radius = int(sec.get("dense_radius_bp", 50))
    dense_min_het = int(sec.get("dense_min_residual_het", 3))
    multi_radius = int(sec.get("multi_indel_radius_bp", 100))
    multi_min_het = int(sec.get("multi_indel_min_residual_het", 3))
    multi_min_indels = int(sec.get("multi_indel_min_indels", 2))

    residual_by_sample = defaultdict(list)
    residual_index_lookup = {}
    for i, row in enumerate(variants):
        if str(row.get("filter_action", "")).upper() == "REMOVE":
            continue
        pos = inum(row.get("source_pos"))
        af = fnum(row.get("source_af"))
        if pos is None or af is None:
            continue
        item = {
            "variant_index": i,
            "sample": str(row.get("sample", "")),
            "species": str(row.get("species", "")),
            "chrom": str(row.get("source_chrom", "")),
            "pos": pos,
            "af": af,
        }
        residual_by_sample[item["sample"]].append(item)
        residual_index_lookup[i] = item

    all_window_rows = []
    missing_vcf = []
    total_rule_windows = Counter()
    total_rule_variants = Counter()

    for sample, residuals in sorted(residual_by_sample.items()):
        if len(residuals) < min(dense_min_het, multi_min_het):
            continue
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
        if not indels:
            continue
        native_length, _, _ = choose_native_length(contig_lengths, residuals, context)

        windows = detect_dense_indel_windows(
            residuals,
            indels,
            mt_length=native_length,
            dense_radius_bp=dense_radius,
            dense_min_het=dense_min_het,
            multi_radius_bp=multi_radius,
            multi_min_het=multi_min_het,
            multi_min_indels=multi_min_indels,
        )
        evidence, window_rows = combine_window_evidence(
            windows,
            residual_index_lookup,
            native_length,
        )
        all_window_rows.extend(window_rows)

        for window in windows:
            for rule in window["rules"]:
                total_rule_windows[rule] += 1

        for idx, item in evidence.items():
            row = variants[idx]
            if not row.get("pre_dense_indel_filter_action"):
                row["pre_dense_indel_filter_action"] = row.get("filter_action", "KEEP") or "KEEP"
                row["pre_dense_indel_filter_reason"] = row.get("filter_reason", "")

            rules = sorted(item["rules"])
            row["dense_indel_filtered"] = "YES"
            row["dense_indel_overlap"] = "YES" if "DENSE_INDEL_OVERLAP" in rules else "NO"
            row["multi_indel_dense"] = "YES" if "MULTI_INDEL_DENSE" in rules else "NO"
            row["dense_indel_rules"] = "+".join(rules)
            row["dense_indel_trigger_positions"] = ",".join(str(x) for x in sorted(item["trigger_positions"]))
            row["dense_indel_min_distance_bp"] = "" if item["min_distance"] is None else str(item["min_distance"])
            row["filter_action"] = "REMOVE"
            row["filter_reason"] = (
                "DENSE_INDEL_OVERLAP"
                if "DENSE_INDEL_OVERLAP" in rules
                else "MULTI_INDEL_DENSE"
            )
            for rule in rules:
                total_rule_variants[rule] += 1

    variant_extra = [
        "pre_dense_indel_filter_action",
        "pre_dense_indel_filter_reason",
        "dense_indel_filtered",
        "dense_indel_overlap",
        "multi_indel_dense",
        "dense_indel_rules",
        "dense_indel_trigger_positions",
        "dense_indel_min_distance_bp",
    ]
    variant_fields = add_fields(list(variants[0]), variant_extra)
    write_tsv(variant_path, variants, variant_fields)

    if samples:
        samples = update_sample_summary(samples, variants)
        sample_fields = add_fields(
            list(samples[0]),
            [
                "n_dense_indel_overlap_variants",
                "n_multi_indel_dense_variants",
                "n_dense_indel_filtered_variants",
                "n_heteroplasmy_variants_to_remove",
                "n_het_after_local_artifact_filter",
                "fraction_het_removed_all",
            ],
        )
        write_tsv(sample_path, samples, sample_fields)

    detail = [dict(row) for row in variants if yes(row.get("dense_indel_filtered"))]
    write_tsv(report_dir / "dense_indel_overlap_variant_detail.tsv", detail, variant_fields)
    write_tsv(
        report_dir / "dense_indel_overlap_windows.tsv",
        all_window_rows,
        [
            "sample", "window_id", "center_indel_pos", "center_indel_af", "rules",
            "n_dense_het_50bp", "n_het_100bp", "n_indels_100bp",
            "n_unique_het_removed", "removed_het_positions",
        ],
    )

    removals = [dict(row) for row in variants if str(row.get("filter_action", "")).upper() == "REMOVE"]
    write_tsv(report_dir / "heteroplasmy_variants_to_remove.tsv", removals, variant_fields)
    write_tsv(report_dir / "numt_variants_to_remove.tsv", removals, variant_fields)

    filtered_variants = {i for i, row in enumerate(variants) if yes(row.get("dense_indel_filtered"))}
    summary = [{
        "dense_radius_bp": dense_radius,
        "dense_min_residual_het": dense_min_het,
        "multi_indel_radius_bp": multi_radius,
        "multi_indel_min_residual_het": multi_min_het,
        "multi_indel_min_indels": multi_min_indels,
        "n_dense_rule_windows": total_rule_windows["DENSE_INDEL_OVERLAP"],
        "n_multi_indel_dense_windows": total_rule_windows["MULTI_INDEL_DENSE"],
        "n_dense_rule_variant_memberships": total_rule_variants["DENSE_INDEL_OVERLAP"],
        "n_multi_indel_dense_variant_memberships": total_rule_variants["MULTI_INDEL_DENSE"],
        "n_unique_dense_indel_filtered_variants": len(filtered_variants),
        "n_samples_dense_indel_filtered": len({str(variants[i].get("sample", "")) for i in filtered_variants}),
        "n_samples_missing_vcf": len(missing_vcf),
        "n_total_heteroplasmy_removals": len(removals),
    }]
    write_tsv(report_dir / "dense_indel_overlap_summary.tsv", summary, list(summary[0]))

    if missing_vcf:
        write_tsv(
            report_dir / "dense_indel_overlap_missing_vcf.tsv",
            [{"sample": sample} for sample in missing_vcf],
            ["sample"],
        )

    print(
        "[dense_indel_overlap] "
        f"dense_windows={total_rule_windows['DENSE_INDEL_OVERLAP']} "
        f"multi_dense_windows={total_rule_windows['MULTI_INDEL_DENSE']} "
        f"unique_removed={len(filtered_variants)} "
        f"samples={len({str(variants[i].get('sample', '')) for i in filtered_variants})} "
        f"total_remove={len(removals)}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
