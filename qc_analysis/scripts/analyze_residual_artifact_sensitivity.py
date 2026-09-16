#!/usr/bin/env python3
"""Report-only sensitivity analysis for residual local HET artifacts.

This module runs after NUMT filtering/expansion and indel-complex detection. It
NEVER changes ``filter_action``. It evaluates two proposed extensions before
promoting either to production filtering:

1) indel-seed AF-coherent expansion
   Existing REMOVE-level indel artifact regions are treated as fixed seeds.
   Residual pre-existing AF-coherent clusters are tested for proximity to the
   seed and similarity of cluster median AF. Expansion is one-step only; newly
   suggested clusters never become seeds, so there is no chain drift.

2) same-species recurrent residual blocks
   Residual HET alleles are searched within species + reference_key. A local
   block must contain multiple identical POS+REF+ALT alleles that are jointly
   present in the same other samples, with cross-sample AF-profile coherence.

Outputs are diagnostic TSVs under local_heteroplasmy_qc/reports. Nothing is
removed or flagged by this script.
"""
from __future__ import annotations

import argparse
import csv
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


def circular_distance(a: int, b: int, length: int | None) -> int:
    linear = abs(a - b)
    if length is None or length <= 0:
        return linear
    if a < 1 or b < 1 or a > length or b > length:
        return linear
    return min(linear, length - linear)


def allele_key(row: dict) -> tuple[str, str, str]:
    return (
        str(row.get("source_pos", "")),
        str(row.get("source_ref", "")),
        str(row.get("source_alt", "")),
    )


def sample_group_key(row: dict) -> tuple[str, str]:
    species = str(row.get("species", ""))
    reference_key = str(row.get("reference_key", ""))
    if not reference_key:
        reference_key = str(row.get("source_chrom", ""))
    return species, reference_key


def load_native_lengths(report_dir: Path) -> dict[str, int | None]:
    rows = read_tsv(report_dir / "indel_complex_native_mt_lengths.tsv")
    out = {}
    for row in rows:
        out[str(row.get("sample", ""))] = inum(row.get("native_mt_length"))
    return out


def residual_clusters(variants: list[dict]) -> list[dict]:
    """Summarize residual variants that still belong to original local clusters."""
    grouped = defaultdict(list)
    for row in variants:
        if str(row.get("filter_action", "")).upper() == "REMOVE":
            continue
        if not yes(row.get("clustered")):
            continue
        cid = str(row.get("cluster_id", ""))
        if not cid:
            continue
        pos = inum(row.get("source_pos"))
        af = fnum(row.get("source_af"))
        if pos is None or af is None:
            continue
        grouped[(str(row.get("sample", "")), cid)].append(row)

    out = []
    for (sample, cid), rows in grouped.items():
        afs = [float(row["source_af"]) for row in rows]
        poss = [int(float(row["source_pos"])) for row in rows]
        first = rows[0]
        out.append(
            {
                "sample": sample,
                "species": str(first.get("species", "")),
                "reference_key": str(first.get("reference_key", "")),
                "cluster_id": cid,
                "n_residual_het": len(rows),
                "positions": poss,
                "median_af": median(afs),
                "af_min": min(afs),
                "af_max": max(afs),
                "af_span": max(afs) - min(afs),
                "rows": rows,
            }
        )
    return out


def indel_seed_regions(variants: list[dict]) -> list[dict]:
    """Build immutable REMOVE-level indel seeds from the current variant detail."""
    grouped = defaultdict(list)
    for row in variants:
        rid = str(row.get("indel_complex_region_id", ""))
        if not rid:
            continue
        reason = str(row.get("filter_reason", ""))
        if str(row.get("filter_action", "")).upper() != "REMOVE":
            continue
        if reason not in {"INDEL_COMPLEX_REGION_HIGH_CONF", "INDEL_AF_MATCHED_VARIANT"}:
            continue
        pos = inum(row.get("source_pos"))
        af = fnum(row.get("source_af"))
        if pos is None or af is None:
            continue
        grouped[(str(row.get("sample", "")), rid)].append(row)

    out = []
    for (sample, rid), rows in grouped.items():
        poss = [int(float(row["source_pos"])) for row in rows]
        afs = [float(row["source_af"]) for row in rows]
        first = rows[0]
        out.append(
            {
                "sample": sample,
                "species": str(first.get("species", "")),
                "region_id": rid,
                "positions": poss,
                "median_af": median(afs),
                "af_min": min(afs),
                "af_max": max(afs),
                "n_removed_seed_het": len(rows),
            }
        )
    return out


def analyze_indel_seed_expansion(
    variants: list[dict],
    length_by_sample: dict[str, int | None],
    distances: list[int],
    deltas: list[float],
) -> tuple[list[dict], list[dict]]:
    clusters = residual_clusters(variants)
    seeds = indel_seed_regions(variants)
    seeds_by_sample = defaultdict(list)
    for seed in seeds:
        seeds_by_sample[seed["sample"]].append(seed)

    max_distance = max(distances)
    max_delta = max(deltas)
    pair_rows = []

    for cluster in clusters:
        sample = cluster["sample"]
        length = length_by_sample.get(sample)
        best = None
        for seed in seeds_by_sample.get(sample, []):
            distance = min(
                circular_distance(cpos, spos, length)
                for cpos in cluster["positions"]
                for spos in seed["positions"]
            )
            delta = abs(cluster["median_af"] - seed["median_af"])
            candidate = (distance, delta, seed["region_id"], seed)
            if best is None or candidate[:3] < best[:3]:
                best = candidate

        if best is None:
            continue
        distance, delta, _, seed = best
        if distance > max_distance or delta > max_delta + 1e-12:
            continue

        row = {
            "sample": sample,
            "species": cluster["species"],
            "cluster_id": cluster["cluster_id"],
            "n_residual_het": cluster["n_residual_het"],
            "cluster_positions": ",".join(str(x) for x in sorted(cluster["positions"])),
            "cluster_median_af": f"{cluster['median_af']:.6g}",
            "cluster_af_span": f"{cluster['af_span']:.6g}",
            "seed_region_id": seed["region_id"],
            "seed_n_removed_het": seed["n_removed_seed_het"],
            "seed_median_af": f"{seed['median_af']:.6g}",
            "nearest_seed_distance_bp": distance,
            "delta_median_af": f"{delta:.6g}",
            "baseline_250bp_delta005": "YES" if distance <= 250 and delta <= 0.05 + 1e-12 else "NO",
        }
        for dist in distances:
            for af_delta in deltas:
                label = f"hit_d{dist}_a{str(af_delta).replace('.', 'p')}"
                row[label] = "YES" if distance <= dist and delta <= af_delta + 1e-12 else "NO"
        pair_rows.append(row)

    sensitivity = []
    for dist in distances:
        for af_delta in deltas:
            hits = [
                row for row in pair_rows
                if int(row["nearest_seed_distance_bp"]) <= dist
                and float(row["delta_median_af"]) <= af_delta + 1e-12
            ]
            sensitivity.append(
                {
                    "distance_bp": dist,
                    "max_delta_median_af": f"{af_delta:.6g}",
                    "n_candidate_clusters": len(hits),
                    "n_candidate_samples": len({row["sample"] for row in hits}),
                    "n_candidate_het": sum(int(row["n_residual_het"]) for row in hits),
                }
            )
    return pair_rows, sensitivity


def residual_variant_universe(variants: list[dict]) -> list[dict]:
    out = []
    for row in variants:
        if str(row.get("filter_action", "")).upper() == "REMOVE":
            continue
        pos = inum(row.get("source_pos"))
        af = fnum(row.get("source_af"))
        if pos is None or af is None:
            continue
        copied = dict(row)
        copied["_pos"] = pos
        copied["_af"] = af
        out.append(copied)
    return out


def common_supporters_for_block(
    sample: str,
    block: list[dict],
    support_samples: dict[tuple[str, str, tuple[str, str, str]], set[str]],
) -> set[str]:
    group = sample_group_key(block[0])
    common = None
    for row in block:
        supporters = set(support_samples.get((group[0], group[1], allele_key(row)), set()))
        supporters.discard(sample)
        common = supporters if common is None else common.intersection(supporters)
        if not common:
            return set()
    return common or set()


def supporter_median_delta(
    sample: str,
    supporter: str,
    block: list[dict],
    af_lookup: dict[tuple[str, str, str, tuple[str, str, str]], float],
) -> float | None:
    group = sample_group_key(block[0])
    diffs = []
    for row in block:
        key = allele_key(row)
        af_other = af_lookup.get((group[0], group[1], supporter, key))
        if af_other is None:
            return None
        diffs.append(abs(float(row["_af"]) - af_other))
    return median(diffs) if diffs else None


def candidate_blocks_for_params(
    variants: list[dict],
    window_bp: int,
    min_shared_variants: int,
    min_other_samples: int,
    max_median_af_delta: float,
) -> list[dict]:
    """Find local recurrent blocks with joint same-sample support in other samples."""
    residual = residual_variant_universe(variants)
    support_samples = defaultdict(set)
    af_lookup = {}
    by_sample_group = defaultdict(list)

    for row in residual:
        species, reference_key = sample_group_key(row)
        sample = str(row.get("sample", ""))
        key = allele_key(row)
        support_samples[(species, reference_key, key)].add(sample)
        af_lookup[(species, reference_key, sample, key)] = float(row["_af"])
        by_sample_group[(species, reference_key, sample)].append(row)

    raw_candidates = []
    for (species, reference_key, sample), rows in by_sample_group.items():
        # Only alleles seen in at least one other sample can contribute.
        recurrent = [
            row for row in rows
            if len(support_samples[(species, reference_key, allele_key(row))] - {sample}) >= min_other_samples
        ]
        recurrent.sort(key=lambda row: (int(row["_pos"]), allele_key(row)))
        n = len(recurrent)
        for i in range(n):
            block = []
            start = int(recurrent[i]["_pos"])
            for j in range(i, n):
                pos = int(recurrent[j]["_pos"])
                if pos - start > window_bp:
                    break
                block.append(recurrent[j])
                if len(block) < min_shared_variants:
                    continue

                common = common_supporters_for_block(sample, block, support_samples)
                if len(common) < min_other_samples:
                    continue

                supporter_deltas = []
                passing_supporters = []
                for supporter in sorted(common):
                    delta = supporter_median_delta(sample, supporter, block, af_lookup)
                    if delta is None:
                        continue
                    supporter_deltas.append((supporter, delta))
                    if delta <= max_median_af_delta + 1e-12:
                        passing_supporters.append((supporter, delta))
                if len(passing_supporters) < min_other_samples:
                    continue

                afs = [float(row["_af"]) for row in block]
                raw_candidates.append(
                    {
                        "sample": sample,
                        "species": species,
                        "reference_key": reference_key,
                        "start": start,
                        "end": int(block[-1]["_pos"]),
                        "n_shared_variants": len(block),
                        "alleles": [allele_key(row) for row in block],
                        "afs": afs,
                        "passing_supporters": passing_supporters,
                    }
                )

    # Collapse nested/overlapping candidates per sample: keep largest block first.
    raw_candidates.sort(
        key=lambda x: (
            x["sample"],
            -x["n_shared_variants"],
            x["end"] - x["start"],
            x["start"],
        )
    )
    accepted = []
    accepted_sets = defaultdict(list)
    for cand in raw_candidates:
        aset = set(cand["alleles"])
        key = (cand["species"], cand["reference_key"], cand["sample"])
        if any(aset.issubset(existing) for existing in accepted_sets[key]):
            continue
        accepted_sets[key].append(aset)
        accepted.append(cand)
    return accepted


def recurrent_block_detail_rows(blocks: list[dict]) -> list[dict]:
    rows = []
    counter = Counter()
    for block in sorted(blocks, key=lambda x: (x["species"], x["sample"], x["start"])):
        counter[block["sample"]] += 1
        block_id = f"{block['sample']}_RRB{counter[block['sample']]:03d}"
        supporters = [sample for sample, _ in block["passing_supporters"]]
        deltas = [delta for _, delta in block["passing_supporters"]]
        allele_text = ",".join(f"{p}:{r}>{a}" for p, r, a in block["alleles"])
        rows.append(
            {
                "block_id": block_id,
                "sample": block["sample"],
                "species": block["species"],
                "reference_key": block["reference_key"],
                "block_start": block["start"],
                "block_end": block["end"],
                "block_span_bp": block["end"] - block["start"],
                "n_shared_variants": block["n_shared_variants"],
                "shared_alleles": allele_text,
                "sample_median_af": f"{median(block['afs']):.6g}",
                "sample_af_span": f"{max(block['afs']) - min(block['afs']):.6g}",
                "n_passing_other_samples": len(supporters),
                "passing_other_samples": ",".join(supporters),
                "median_supporter_profile_delta": f"{median(deltas):.6g}" if deltas else "",
                "max_supporter_profile_delta": f"{max(deltas):.6g}" if deltas else "",
            }
        )
    return rows


def analyze_recurrent_block_sensitivity(
    variants: list[dict],
    windows: list[int],
    min_shared_values: list[int],
    min_other_values: list[int],
    af_deltas: list[float],
) -> list[dict]:
    summary = []
    for window in windows:
        for min_shared in min_shared_values:
            for min_other in min_other_values:
                for delta in af_deltas:
                    blocks = candidate_blocks_for_params(
                        variants,
                        window_bp=window,
                        min_shared_variants=min_shared,
                        min_other_samples=min_other,
                        max_median_af_delta=delta,
                    )
                    unique_alleles = set()
                    for block in blocks:
                        for allele in block["alleles"]:
                            unique_alleles.add((block["sample"],) + allele)
                    summary.append(
                        {
                            "window_bp": window,
                            "min_shared_variants": min_shared,
                            "min_other_samples": min_other,
                            "max_median_af_delta": f"{delta:.6g}",
                            "n_candidate_blocks": len(blocks),
                            "n_candidate_samples": len({block["sample"] for block in blocks}),
                            "n_candidate_variant_memberships": sum(block["n_shared_variants"] for block in blocks),
                            "n_unique_sample_variants": len(unique_alleles),
                        }
                    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()

    cfg = read_simple_yaml(args.config)
    local_sec = cfg.get("local_heteroplasmy_qc", {})
    sec = cfg.get("residual_artifact_sensitivity", {})
    if sec.get("enabled", True) is False:
        print("[residual_artifact_sensitivity] disabled", file=sys.stderr)
        return 0

    report_dir = resolve(local_sec.get("output_dir", "results/qc/local_heteroplasmy_qc")) / "reports"
    variant_path = report_dir / "local_heteroplasmy_variant_detail.tsv"
    variants = read_tsv(variant_path)
    if not variants:
        raise RuntimeError(f"No variant rows found: {variant_path}")

    lengths = load_native_lengths(report_dir)

    # Report-only indel expansion grid.
    indel_distances = [int(x) for x in sec.get("indel_expansion_distance_grid", [100, 150, 250, 500])]
    indel_deltas = [float(x) for x in sec.get("indel_expansion_delta_af_grid", [0.03, 0.05, 0.07, 0.10])]
    indel_candidates, indel_sensitivity = analyze_indel_seed_expansion(
        variants,
        lengths,
        distances=indel_distances,
        deltas=indel_deltas,
    )

    indel_fields = [
        "sample", "species", "cluster_id", "n_residual_het", "cluster_positions",
        "cluster_median_af", "cluster_af_span", "seed_region_id", "seed_n_removed_het",
        "seed_median_af", "nearest_seed_distance_bp", "delta_median_af", "baseline_250bp_delta005",
    ]
    for dist in indel_distances:
        for delta in indel_deltas:
            indel_fields.append(f"hit_d{dist}_a{str(delta).replace('.', 'p')}")
    write_tsv(report_dir / "indel_seed_expansion_candidates_report_only.tsv", indel_candidates, indel_fields)
    write_tsv(
        report_dir / "indel_seed_expansion_sensitivity.tsv",
        indel_sensitivity,
        ["distance_bp", "max_delta_median_af", "n_candidate_clusters", "n_candidate_samples", "n_candidate_het"],
    )

    # Report-only recurrent residual block grid.
    recurrent_windows = [int(x) for x in sec.get("recurrent_block_window_grid", [100, 250, 500])]
    recurrent_shared = [int(x) for x in sec.get("recurrent_block_min_shared_grid", [3, 4, 5])]
    recurrent_other = [int(x) for x in sec.get("recurrent_block_min_other_samples_grid", [1, 2])]
    recurrent_deltas = [float(x) for x in sec.get("recurrent_block_af_delta_grid", [0.03, 0.05, 0.10])]

    recurrent_sensitivity = analyze_recurrent_block_sensitivity(
        variants,
        windows=recurrent_windows,
        min_shared_values=recurrent_shared,
        min_other_values=recurrent_other,
        af_deltas=recurrent_deltas,
    )
    write_tsv(
        report_dir / "species_recurrent_residual_block_sensitivity.tsv",
        recurrent_sensitivity,
        [
            "window_bp", "min_shared_variants", "min_other_samples", "max_median_af_delta",
            "n_candidate_blocks", "n_candidate_samples", "n_candidate_variant_memberships", "n_unique_sample_variants",
        ],
    )

    # Baseline diagnostic detail: deliberately strict joint recurrence.
    baseline_window = int(sec.get("recurrent_block_baseline_window_bp", 250))
    baseline_shared = int(sec.get("recurrent_block_baseline_min_shared", 3))
    baseline_other = int(sec.get("recurrent_block_baseline_min_other_samples", 2))
    baseline_delta = float(sec.get("recurrent_block_baseline_max_median_af_delta", 0.05))
    baseline_blocks = candidate_blocks_for_params(
        variants,
        window_bp=baseline_window,
        min_shared_variants=baseline_shared,
        min_other_samples=baseline_other,
        max_median_af_delta=baseline_delta,
    )
    recurrent_detail = recurrent_block_detail_rows(baseline_blocks)
    write_tsv(
        report_dir / "species_recurrent_residual_blocks_report_only.tsv",
        recurrent_detail,
        [
            "block_id", "sample", "species", "reference_key", "block_start", "block_end", "block_span_bp",
            "n_shared_variants", "shared_alleles", "sample_median_af", "sample_af_span",
            "n_passing_other_samples", "passing_other_samples", "median_supporter_profile_delta", "max_supporter_profile_delta",
        ],
    )

    baseline_indel = [row for row in indel_candidates if row["baseline_250bp_delta005"] == "YES"]
    summary = [{
        "analysis_mode": "REPORT_ONLY_NO_FILTER_ACTION_CHANGES",
        "indel_baseline_distance_bp": 250,
        "indel_baseline_max_delta_af": 0.05,
        "n_indel_expansion_candidate_clusters": len(baseline_indel),
        "n_indel_expansion_candidate_samples": len({row["sample"] for row in baseline_indel}),
        "n_indel_expansion_candidate_het": sum(int(row["n_residual_het"]) for row in baseline_indel),
        "recurrent_baseline_window_bp": baseline_window,
        "recurrent_baseline_min_shared": baseline_shared,
        "recurrent_baseline_min_other_samples": baseline_other,
        "recurrent_baseline_max_median_af_delta": f"{baseline_delta:.6g}",
        "n_recurrent_residual_blocks": len(recurrent_detail),
        "n_recurrent_residual_samples": len({row["sample"] for row in recurrent_detail}),
        "n_recurrent_residual_variant_memberships": sum(int(row["n_shared_variants"]) for row in recurrent_detail),
    }]
    write_tsv(report_dir / "residual_artifact_sensitivity_summary.tsv", summary, list(summary[0]))

    print(
        "[residual_artifact_sensitivity] REPORT ONLY "
        f"indel_baseline_clusters={len(baseline_indel)} "
        f"indel_baseline_het={sum(int(row['n_residual_het']) for row in baseline_indel)} "
        f"recurrent_baseline_blocks={len(recurrent_detail)} "
        f"recurrent_baseline_samples={len({row['sample'] for row in recurrent_detail})}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
