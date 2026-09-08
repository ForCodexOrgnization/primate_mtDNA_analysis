"""Same-species recurrence assessment for native-coordinate heteroplasmies."""
from __future__ import annotations

from collections import defaultdict
from statistics import median


def variant_key(row: dict) -> tuple[str, int, str, str]:
    return (
        row.get("source_chrom", ""),
        int(row["source_pos"]),
        row.get("source_ref", ""),
        row.get("source_alt", ""),
    )


def build_species_sample_index(rows: list[dict], sample_to_species: dict[str, str]) -> dict[str, dict[str, dict[tuple, float]]]:
    out: dict[str, dict[str, dict[tuple, float]]] = defaultdict(lambda: defaultdict(dict))
    for row in rows:
        sample = row["sample"]
        species = sample_to_species.get(sample, "")
        if not species:
            continue
        try:
            af = float(row["source_af"])
        except (TypeError, ValueError):
            continue
        out[species][sample][variant_key(row)] = af
    return out


def assess_recurrence(
    cluster_rows: list[dict],
    sample: str,
    species: str,
    index: dict[str, dict[str, dict[tuple, float]]],
    min_shared_variants: int = 2,
    min_shared_fraction: float = 0.50,
    max_median_af_difference: float = 0.05,
) -> dict:
    cluster = {variant_key(row): float(row["source_af"]) for row in cluster_rows}
    best = None
    for other_sample, other in index.get(species, {}).items():
        if other_sample == sample:
            continue
        shared = sorted(set(cluster).intersection(other))
        if not shared:
            continue
        fraction = len(shared) / len(cluster) if cluster else 0.0
        af_differences = [abs(cluster[key] - other[key]) for key in shared]
        af_difference = median(af_differences) if af_differences else float("inf")
        recurrent = (
            len(shared) >= min_shared_variants
            and fraction >= min_shared_fraction
            and af_difference <= max_median_af_difference
        )
        rank = (1 if recurrent else 0, len(shared), fraction, -af_difference)
        if best is None or rank > best[0]:
            best = (rank, other_sample, shared, fraction, af_difference, recurrent)
    if best is None:
        return {
            "recurrent": False,
            "best_recurrent_sample": "",
            "shared_variants": 0,
            "shared_fraction": 0.0,
            "median_af_difference": "NA",
        }
    _, other_sample, shared, fraction, af_difference, recurrent = best
    return {
        "recurrent": recurrent,
        "best_recurrent_sample": other_sample,
        "shared_variants": len(shared),
        "shared_fraction": fraction,
        "median_af_difference": af_difference,
    }
