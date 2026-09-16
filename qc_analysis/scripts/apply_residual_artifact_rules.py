#!/usr/bin/env python3
"""Apply validated residual heteroplasmy artifact rules.

This production step runs after indel-complex detection. Two rules are applied:

1) Indel-seed cluster expansion (REMOVE)
   - fixed REMOVE-level indel artifact seeds only;
   - residual pre-existing local cluster with >=3 residual HETs;
   - nearest distance to seed <=150 bp;
   - |cluster median AF - seed median AF| <=0.03;
   - one expansion round only (no chain drift).

2) Strong same-species recurrent residual block (FLAG only)
   - same species + reference_key;
   - >=4 identical POS+REF+ALT variants within <=250 bp;
   - the whole block is jointly present in >=2 OTHER samples;
   - median across-variant AF-profile difference <=0.03 in those supporters;
   - annotate / FLAG only; never REMOVE by recurrence alone.

All decisions are made in immutable native SOURCE coordinates. The script updates
local heteroplasmy reports and the combined removal list consumed downstream.
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


def append_reason(existing: str, reason: str) -> str:
    parts = [item for item in str(existing or "").split(";") if item]
    if reason and reason not in parts:
        parts.append(reason)
    return ";".join(parts)


def circular_distance(a: int, b: int, length: int | None) -> int:
    linear = abs(a - b)
    if length is None or length <= 0:
        return linear
    if a < 1 or b < 1 or a > length or b > length:
        return linear
    return min(linear, length - linear)


def sample_group_key(row: dict) -> tuple[str, str]:
    species = str(row.get("species", ""))
    reference_key = str(row.get("reference_key", ""))
    if not reference_key:
        reference_key = str(row.get("source_chrom", ""))
    return species, reference_key


def allele_key(row: dict) -> tuple[str, str, str]:
    return (
        str(row.get("source_pos", "")),
        str(row.get("source_ref", "")),
        str(row.get("source_alt", "")),
    )


def load_native_lengths(report_dir: Path) -> dict[str, int | None]:
    out = {}
    for row in read_tsv(report_dir / "indel_complex_native_mt_lengths.tsv"):
        out[str(row.get("sample", ""))] = inum(row.get("native_mt_length"))
    return out


def reset_previous_production_calls(variants: list[dict]) -> None:
    """Restore prior actions from this script so reruns are deterministic."""
    # Recurrent block is applied second, so restore it first.
    for row in variants:
        if yes(row.get("species_recurrent_residual_block")):
            row["filter_action"] = row.get("pre_recurrent_block_filter_action", "KEEP") or "KEEP"
            row["filter_reason"] = row.get("pre_recurrent_block_filter_reason", "")
            row["species_recurrent_residual_block"] = "NO"
            for field in (
                "species_recurrent_residual_block_ids",
                "species_recurrent_support_samples",
                "species_recurrent_n_other_samples",
                "species_recurrent_profile_delta",
            ):
                row[field] = ""

    for row in variants:
        if yes(row.get("indel_seed_expanded")):
            row["filter_action"] = row.get("pre_indel_seed_filter_action", "KEEP") or "KEEP"
            row["filter_reason"] = row.get("pre_indel_seed_filter_reason", "")
            row["indel_seed_expanded"] = "NO"
            for field in (
                "indel_seed_expansion_seed_region_id",
                "indel_seed_expansion_distance_bp",
                "indel_seed_expansion_delta_af",
                "indel_seed_expansion_cluster_median_af",
            ):
                row[field] = ""


def fixed_indel_seeds(variants: list[dict]) -> list[dict]:
    """Build immutable seeds from pre-expansion REMOVE-level indel calls."""
    grouped = defaultdict(list)
    allowed_reasons = {"INDEL_COMPLEX_REGION_HIGH_CONF", "INDEL_AF_MATCHED_VARIANT"}
    for row in variants:
        if str(row.get("filter_action", "")).upper() != "REMOVE":
            continue
        if str(row.get("filter_reason", "")) not in allowed_reasons:
            continue
        region_id = str(row.get("indel_complex_region_id", ""))
        pos = inum(row.get("source_pos"))
        af = fnum(row.get("source_af"))
        if not region_id or pos is None or af is None:
            continue
        grouped[(str(row.get("sample", "")), region_id)].append(row)

    seeds = []
    for (sample, region_id), rows in grouped.items():
        poss = [int(float(row["source_pos"])) for row in rows]
        afs = [float(row["source_af"]) for row in rows]
        seeds.append(
            {
                "sample": sample,
                "region_id": region_id,
                "positions": poss,
                "median_af": median(afs),
                "n_seed_het": len(rows),
            }
        )
    return seeds


def residual_cluster_groups(variants: list[dict], min_residual_het: int) -> list[dict]:
    grouped = defaultdict(list)
    for i, row in enumerate(variants):
        if str(row.get("filter_action", "")).upper() == "REMOVE":
            continue
        if not yes(row.get("clustered")):
            continue
        cluster_id = str(row.get("cluster_id", ""))
        pos = inum(row.get("source_pos"))
        af = fnum(row.get("source_af"))
        if not cluster_id or pos is None or af is None:
            continue
        grouped[(str(row.get("sample", "")), cluster_id)].append((i, row))

    clusters = []
    for (sample, cluster_id), indexed_rows in grouped.items():
        if len(indexed_rows) < min_residual_het:
            continue
        rows = [row for _, row in indexed_rows]
        afs = [float(row["source_af"]) for row in rows]
        poss = [int(float(row["source_pos"])) for row in rows]
        first = rows[0]
        clusters.append(
            {
                "sample": sample,
                "species": str(first.get("species", "")),
                "cluster_id": cluster_id,
                "indices": [i for i, _ in indexed_rows],
                "positions": poss,
                "median_af": median(afs),
                "af_span": max(afs) - min(afs),
                "n_residual_het": len(rows),
            }
        )
    return clusters


def apply_indel_seed_expansion(
    variants: list[dict],
    length_by_sample: dict[str, int | None],
    max_distance_bp: int = 150,
    max_delta_af: float = 0.03,
    min_residual_het: int = 3,
) -> list[dict]:
    """One-round REMOVE expansion from fixed indel seeds; no chain drift."""
    seeds = fixed_indel_seeds(variants)
    seeds_by_sample = defaultdict(list)
    for seed in seeds:
        seeds_by_sample[seed["sample"]].append(seed)

    clusters = residual_cluster_groups(variants, min_residual_het=min_residual_het)
    promoted = []
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
        if distance > max_distance_bp or delta > max_delta_af + 1e-12:
            continue

        for i in cluster["indices"]:
            row = variants[i]
            if not row.get("pre_indel_seed_filter_action"):
                row["pre_indel_seed_filter_action"] = row.get("filter_action", "KEEP") or "KEEP"
                row["pre_indel_seed_filter_reason"] = row.get("filter_reason", "")
            row["indel_seed_expanded"] = "YES"
            row["indel_seed_expansion_seed_region_id"] = seed["region_id"]
            row["indel_seed_expansion_distance_bp"] = str(distance)
            row["indel_seed_expansion_delta_af"] = f"{delta:.6g}"
            row["indel_seed_expansion_cluster_median_af"] = f"{cluster['median_af']:.6g}"
            row["filter_action"] = "REMOVE"
            row["filter_reason"] = "INDEL_SEED_CLUSTER_EXPANSION"

        promoted.append(
            {
                "sample": sample,
                "species": cluster["species"],
                "cluster_id": cluster["cluster_id"],
                "n_removed_het": cluster["n_residual_het"],
                "cluster_positions": ",".join(str(x) for x in sorted(cluster["positions"])),
                "cluster_median_af": f"{cluster['median_af']:.6g}",
                "cluster_af_span": f"{cluster['af_span']:.6g}",
                "seed_region_id": seed["region_id"],
                "seed_n_het": seed["n_seed_het"],
                "seed_median_af": f"{seed['median_af']:.6g}",
                "nearest_seed_distance_bp": distance,
                "delta_median_af": f"{delta:.6g}",
            }
        )
    return promoted


def residual_variant_universe(variants: list[dict]) -> list[dict]:
    out = []
    for i, row in enumerate(variants):
        if str(row.get("filter_action", "")).upper() == "REMOVE":
            continue
        pos = inum(row.get("source_pos"))
        af = fnum(row.get("source_af"))
        if pos is None or af is None:
            continue
        copied = dict(row)
        copied["_index"] = i
        copied["_pos"] = pos
        copied["_af"] = af
        out.append(copied)
    return out


def common_supporters(
    sample: str,
    block: list[dict],
    support_samples: dict[tuple[str, str, tuple[str, str, str]], set[str]],
) -> set[str]:
    species, reference_key = sample_group_key(block[0])
    common = None
    for row in block:
        supporters = set(support_samples.get((species, reference_key, allele_key(row)), set()))
        supporters.discard(sample)
        common = supporters if common is None else common.intersection(supporters)
        if not common:
            return set()
    return common or set()


def supporter_profile_delta(
    sample: str,
    supporter: str,
    block: list[dict],
    af_lookup: dict[tuple[str, str, str, tuple[str, str, str]], float],
) -> float | None:
    species, reference_key = sample_group_key(block[0])
    diffs = []
    for row in block:
        af_other = af_lookup.get((species, reference_key, supporter, allele_key(row)))
        if af_other is None:
            return None
        diffs.append(abs(float(row["_af"]) - af_other))
    return median(diffs) if diffs else None


def find_recurrent_residual_blocks(
    variants: list[dict],
    window_bp: int = 250,
    min_shared_variants: int = 4,
    min_other_samples: int = 2,
    max_median_af_delta: float = 0.03,
) -> list[dict]:
    """Find strict joint recurrence blocks in the post-REMOVE residual universe."""
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

                common = common_supporters(sample, block, support_samples)
                if len(common) < min_other_samples:
                    continue

                passing = []
                for supporter in sorted(common):
                    delta = supporter_profile_delta(sample, supporter, block, af_lookup)
                    if delta is not None and delta <= max_median_af_delta + 1e-12:
                        passing.append((supporter, delta))
                if len(passing) < min_other_samples:
                    continue

                raw_candidates.append(
                    {
                        "sample": sample,
                        "species": species,
                        "reference_key": reference_key,
                        "start": start,
                        "end": int(block[-1]["_pos"]),
                        "rows": list(block),
                        "alleles": [allele_key(row) for row in block],
                        "afs": [float(row["_af"]) for row in block],
                        "passing_supporters": passing,
                    }
                )

    # Keep maximal blocks and remove nested duplicates within each sample.
    raw_candidates.sort(
        key=lambda x: (
            x["sample"],
            -len(x["alleles"]),
            x["end"] - x["start"],
            x["start"],
        )
    )
    accepted = []
    accepted_sets = defaultdict(list)
    for cand in raw_candidates:
        aset = set(cand["alleles"])
        group = (cand["species"], cand["reference_key"], cand["sample"])
        if any(aset.issubset(existing) for existing in accepted_sets[group]):
            continue
        accepted_sets[group].append(aset)
        accepted.append(cand)
    return accepted


def apply_recurrent_block_flags(variants: list[dict], blocks: list[dict]) -> list[dict]:
    """Annotate / FLAG recurrent blocks. Never convert a row to REMOVE."""
    by_index = defaultdict(list)
    block_rows = []
    sample_counter = Counter()

    for block in sorted(blocks, key=lambda x: (x["species"], x["sample"], x["start"])):
        sample_counter[block["sample"]] += 1
        block_id = f"{block['sample']}_RRB{sample_counter[block['sample']]:03d}"
        supporters = [sample for sample, _ in block["passing_supporters"]]
        deltas = [delta for _, delta in block["passing_supporters"]]
        for row in block["rows"]:
            by_index[int(row["_index"])].append(
                {
                    "block_id": block_id,
                    "supporters": supporters,
                    "n_supporters": len(supporters),
                    "profile_delta": median(deltas) if deltas else None,
                }
            )

        block_rows.append(
            {
                "block_id": block_id,
                "sample": block["sample"],
                "species": block["species"],
                "reference_key": block["reference_key"],
                "block_start": block["start"],
                "block_end": block["end"],
                "block_span_bp": block["end"] - block["start"],
                "n_shared_variants": len(block["alleles"]),
                "shared_alleles": ",".join(f"{p}:{r}>{a}" for p, r, a in block["alleles"]),
                "sample_median_af": f"{median(block['afs']):.6g}",
                "sample_af_span": f"{max(block['afs']) - min(block['afs']):.6g}",
                "n_passing_other_samples": len(supporters),
                "passing_other_samples": ",".join(supporters),
                "median_supporter_profile_delta": f"{median(deltas):.6g}" if deltas else "",
                "max_supporter_profile_delta": f"{max(deltas):.6g}" if deltas else "",
                "action": "FLAG",
            }
        )

    for i, memberships in by_index.items():
        row = variants[i]
        if str(row.get("filter_action", "")).upper() == "REMOVE":
            continue
        if not row.get("pre_recurrent_block_filter_action"):
            row["pre_recurrent_block_filter_action"] = row.get("filter_action", "KEEP") or "KEEP"
            row["pre_recurrent_block_filter_reason"] = row.get("filter_reason", "")
        row["species_recurrent_residual_block"] = "YES"
        row["species_recurrent_residual_block_ids"] = ";".join(sorted({m["block_id"] for m in memberships}))
        row["species_recurrent_support_samples"] = ";".join(sorted({s for m in memberships for s in m["supporters"]}))
        row["species_recurrent_n_other_samples"] = str(max(m["n_supporters"] for m in memberships))
        profile_values = [m["profile_delta"] for m in memberships if m["profile_delta"] is not None]
        row["species_recurrent_profile_delta"] = f"{min(profile_values):.6g}" if profile_values else ""
        if str(row.get("filter_action", "")).upper() == "KEEP":
            row["filter_action"] = "FLAG"
        row["filter_reason"] = append_reason(
            row.get("filter_reason", ""),
            "STRONG_SPECIES_RECURRENT_RESIDUAL_BLOCK",
        )
    return block_rows


def update_samples(samples: list[dict], variants: list[dict], expanded_clusters: list[dict], recurrent_blocks: list[dict]) -> list[dict]:
    expanded_variant_n = Counter()
    recurrent_variant_n = Counter()
    recurrent_block_n = Counter(row["sample"] for row in recurrent_blocks)
    total_remove = Counter()
    residual = Counter()

    for row in variants:
        sample = str(row.get("sample", ""))
        if yes(row.get("indel_seed_expanded")):
            expanded_variant_n[sample] += 1
        if yes(row.get("species_recurrent_residual_block")):
            recurrent_variant_n[sample] += 1
        if str(row.get("filter_action", "")).upper() == "REMOVE":
            total_remove[sample] += 1
        else:
            residual[sample] += 1

    expanded_cluster_n = Counter(row["sample"] for row in expanded_clusters)
    out = []
    for raw in samples:
        row = dict(raw)
        sample = str(row.get("sample", ""))
        n_het = inum(row.get("n_het")) or 0
        row["n_indel_seed_expanded_clusters"] = expanded_cluster_n[sample]
        row["n_indel_seed_expanded_variants"] = expanded_variant_n[sample]
        row["n_species_recurrent_residual_blocks"] = recurrent_block_n[sample]
        row["n_species_recurrent_residual_variants"] = recurrent_variant_n[sample]
        row["n_heteroplasmy_variants_to_remove"] = total_remove[sample]
        row["n_het_after_local_artifact_filter"] = residual[sample]
        row["fraction_het_removed_all"] = f"{total_remove[sample] / n_het:.6g}" if n_het else "0"
        out.append(row)
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()

    cfg = read_simple_yaml(args.config)
    local_sec = cfg.get("local_heteroplasmy_qc", {})
    sec = cfg.get("residual_artifact_production", {})
    if sec.get("enabled", True) is False:
        print("[residual_artifact_production] disabled", file=sys.stderr)
        return 0

    report_dir = resolve(local_sec.get("output_dir", "results/qc/local_heteroplasmy_qc")) / "reports"
    variant_path = report_dir / "local_heteroplasmy_variant_detail.tsv"
    sample_path = report_dir / "local_heteroplasmy_sample_summary.tsv"
    variants = read_tsv(variant_path)
    samples = read_tsv(sample_path)
    if not variants:
        raise RuntimeError(f"No variant rows found: {variant_path}")

    reset_previous_production_calls(variants)
    lengths = load_native_lengths(report_dir)

    # Validated production thresholds from cohort sensitivity analysis.
    indel_distance = int(sec.get("indel_seed_expansion_distance_bp", 150))
    indel_delta = float(sec.get("indel_seed_expansion_max_delta_af", 0.03))
    indel_min_het = int(sec.get("indel_seed_expansion_min_residual_het", 3))

    recurrent_window = int(sec.get("recurrent_block_window_bp", 250))
    recurrent_min_shared = int(sec.get("recurrent_block_min_shared_variants", 4))
    recurrent_min_other = int(sec.get("recurrent_block_min_other_samples", 2))
    recurrent_delta = float(sec.get("recurrent_block_max_median_af_delta", 0.03))

    expanded_clusters = apply_indel_seed_expansion(
        variants,
        lengths,
        max_distance_bp=indel_distance,
        max_delta_af=indel_delta,
        min_residual_het=indel_min_het,
    )

    recurrent_blocks_raw = find_recurrent_residual_blocks(
        variants,
        window_bp=recurrent_window,
        min_shared_variants=recurrent_min_shared,
        min_other_samples=recurrent_min_other,
        max_median_af_delta=recurrent_delta,
    )
    recurrent_blocks = apply_recurrent_block_flags(variants, recurrent_blocks_raw)

    variant_extra = [
        "pre_indel_seed_filter_action",
        "pre_indel_seed_filter_reason",
        "indel_seed_expanded",
        "indel_seed_expansion_seed_region_id",
        "indel_seed_expansion_distance_bp",
        "indel_seed_expansion_delta_af",
        "indel_seed_expansion_cluster_median_af",
        "pre_recurrent_block_filter_action",
        "pre_recurrent_block_filter_reason",
        "species_recurrent_residual_block",
        "species_recurrent_residual_block_ids",
        "species_recurrent_support_samples",
        "species_recurrent_n_other_samples",
        "species_recurrent_profile_delta",
    ]
    variant_fields = add_fields(list(variants[0]), variant_extra)

    sample_extra = [
        "n_indel_seed_expanded_clusters",
        "n_indel_seed_expanded_variants",
        "n_species_recurrent_residual_blocks",
        "n_species_recurrent_residual_variants",
        "n_heteroplasmy_variants_to_remove",
        "n_het_after_local_artifact_filter",
        "fraction_het_removed_all",
    ]
    sample_fields = add_fields(list(samples[0]) if samples else [], sample_extra)

    write_tsv(variant_path, variants, variant_fields)
    if samples:
        samples = update_samples(samples, variants, expanded_clusters, recurrent_blocks)
        write_tsv(sample_path, samples, sample_fields)

    write_tsv(
        report_dir / "indel_seed_expansion_applied.tsv",
        expanded_clusters,
        [
            "sample", "species", "cluster_id", "n_removed_het", "cluster_positions",
            "cluster_median_af", "cluster_af_span", "seed_region_id", "seed_n_het",
            "seed_median_af", "nearest_seed_distance_bp", "delta_median_af",
        ],
    )
    write_tsv(
        report_dir / "species_recurrent_residual_blocks.tsv",
        recurrent_blocks,
        [
            "block_id", "sample", "species", "reference_key", "block_start", "block_end", "block_span_bp",
            "n_shared_variants", "shared_alleles", "sample_median_af", "sample_af_span",
            "n_passing_other_samples", "passing_other_samples", "median_supporter_profile_delta",
            "max_supporter_profile_delta", "action",
        ],
    )

    removals = [dict(row) for row in variants if str(row.get("filter_action", "")).upper() == "REMOVE"]
    write_tsv(report_dir / "heteroplasmy_variants_to_remove.tsv", removals, variant_fields)
    # Legacy name retained for current downstream final-filter wrapper.
    write_tsv(report_dir / "numt_variants_to_remove.tsv", removals, variant_fields)

    summary = [{
        "indel_seed_expansion_distance_bp": indel_distance,
        "indel_seed_expansion_max_delta_af": f"{indel_delta:.6g}",
        "indel_seed_expansion_min_residual_het": indel_min_het,
        "n_indel_seed_expanded_clusters": len(expanded_clusters),
        "n_indel_seed_expanded_samples": len({row["sample"] for row in expanded_clusters}),
        "n_indel_seed_expanded_variants": sum(int(row["n_removed_het"]) for row in expanded_clusters),
        "recurrent_block_window_bp": recurrent_window,
        "recurrent_block_min_shared_variants": recurrent_min_shared,
        "recurrent_block_min_other_samples": recurrent_min_other,
        "recurrent_block_max_median_af_delta": f"{recurrent_delta:.6g}",
        "recurrent_block_action": "FLAG",
        "n_species_recurrent_residual_blocks": len(recurrent_blocks),
        "n_species_recurrent_residual_samples": len({row["sample"] for row in recurrent_blocks}),
        "n_species_recurrent_residual_variant_memberships": sum(int(row["n_shared_variants"]) for row in recurrent_blocks),
        "n_total_heteroplasmy_removals": len(removals),
    }]
    write_tsv(report_dir / "residual_artifact_production_summary.tsv", summary, list(summary[0]))

    print(
        "[residual_artifact_production] "
        f"indel_expanded_clusters={len(expanded_clusters)} "
        f"indel_expanded_variants={sum(int(row['n_removed_het']) for row in expanded_clusters)} "
        f"recurrent_flag_blocks={len(recurrent_blocks)} "
        f"recurrent_flag_samples={len({row['sample'] for row in recurrent_blocks})} "
        f"total_remove={len(removals)}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
