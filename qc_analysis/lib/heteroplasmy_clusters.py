"""Native-coordinate local heteroplasmy cluster detection.

The implementation intentionally keeps cluster discovery independent from NUMT
annotation. It uses a circular mitochondrial coordinate system, AF-coherent
windows, sample-specific permutation nulls, greedy independent seed selection,
and post-hoc cluster expansion.
"""
from __future__ import annotations

import random
from dataclasses import dataclass
from statistics import median
from typing import Iterable, Sequence


@dataclass(frozen=True)
class ClusterConfig:
    mt_length: int = 16569
    window_bp: int = 250
    af_span_max: float = 0.06
    min_seed_variants: int = 3
    simulations: int = 2000
    empirical_p_max: float = 0.01
    random_seed: int = 20260904


def circular_distance(a: int, b: int, length: int) -> int:
    d = abs(a - b)
    return min(d, length - d)


def _window_members(positions: Sequence[int], anchor: int, length: int, window: int) -> list[int]:
    """Indices lying within a forward circular window [anchor, anchor + window]."""
    return [i for i, pos in enumerate(positions) if (pos - anchor) % length <= window]


def _af_coherent_subsets(
    indices: Sequence[int],
    afs: Sequence[float],
    max_span: float,
    minimum: int = 1,
) -> list[tuple[int, ...]]:
    """Enumerate maximal AF-sorted windows satisfying the AF-span constraint."""
    if not indices:
        return []
    ordered = sorted(indices, key=lambda i: (afs[i], i))
    out: list[tuple[int, ...]] = []
    seen: set[tuple[int, ...]] = set()
    right = 0
    for left in range(len(ordered)):
        right = max(right, left)
        while right + 1 < len(ordered) and afs[ordered[right + 1]] - afs[ordered[left]] <= max_span:
            right += 1
        candidate = tuple(sorted(ordered[left : right + 1]))
        if len(candidate) >= minimum and candidate not in seen:
            seen.add(candidate)
            out.append(candidate)
    return out


def _largest_af_coherent(indices: Sequence[int], afs: Sequence[float], max_span: float) -> tuple[int, ...]:
    candidates = _af_coherent_subsets(indices, afs, max_span, minimum=1)
    if not candidates:
        return ()
    return max(candidates, key=lambda x: (len(x), tuple(-i for i in x)))


def candidate_seeds(positions: Sequence[int], afs: Sequence[float], cfg: ClusterConfig) -> list[tuple[int, ...]]:
    """Collect all unique AF-coherent candidates from every circular positional window."""
    seen: set[tuple[int, ...]] = set()
    seeds: list[tuple[int, ...]] = []
    for anchor in positions:
        members = _window_members(positions, anchor, cfg.mt_length, cfg.window_bp)
        for seed in _af_coherent_subsets(members, afs, cfg.af_span_max, cfg.min_seed_variants):
            if seed not in seen:
                seen.add(seed)
                seeds.append(seed)
    seeds.sort(key=lambda x: (-len(x), min(positions[i] for i in x), x))
    return seeds


def max_af_coherent_count(positions: Sequence[int], afs: Sequence[float], cfg: ClusterConfig) -> int:
    if not positions:
        return 0
    best = 0
    for anchor in positions:
        members = _window_members(positions, anchor, cfg.mt_length, cfg.window_bp)
        best = max(best, len(_largest_af_coherent(members, afs, cfg.af_span_max)))
    return best


def permutation_null(afs: Sequence[float], cfg: ClusterConfig, seed_offset: int = 0) -> list[int]:
    n = len(afs)
    if n == 0:
        return []
    if n > cfg.mt_length:
        raise ValueError(f"cannot sample {n} unique mtDNA positions from length {cfg.mt_length}")
    rng = random.Random(cfg.random_seed + seed_offset)
    population = range(1, cfg.mt_length + 1)
    return [
        max_af_coherent_count(rng.sample(population, n), afs, cfg)
        for _ in range(cfg.simulations)
    ]


def empirical_p(observed: int, null: Sequence[int]) -> float | None:
    if not null:
        return None
    return (1 + sum(value >= observed for value in null)) / (len(null) + 1)


def critical_count(null: Sequence[int], p_max: float, minimum: int) -> int:
    if not null:
        return minimum
    upper = max(null) + 2
    for k in range(minimum, upper + 1):
        tail = (1 + sum(value >= k for value in null)) / (len(null) + 1)
        if tail <= p_max:
            return k
    return upper + 1


def independent_seeds(
    positions: Sequence[int],
    afs: Sequence[float],
    threshold: int,
    cfg: ClusterConfig,
) -> list[list[int]]:
    selected: list[list[int]] = []
    used: set[int] = set()
    for seed in candidate_seeds(positions, afs, cfg):
        if len(seed) < threshold or used.intersection(seed):
            continue
        selected.append(list(seed))
        used.update(seed)
    return selected


def expand_clusters(
    clusters: Sequence[Sequence[int]],
    positions: Sequence[int],
    afs: Sequence[float],
    cfg: ClusterConfig,
) -> list[list[int]]:
    expanded = [list(cluster) for cluster in clusters]
    assigned = {i for cluster in expanded for i in cluster}
    changed = True
    while changed:
        changed = False
        for idx in range(len(positions)):
            if idx in assigned:
                continue
            for cluster in expanded:
                close = any(
                    circular_distance(positions[idx], positions[j], cfg.mt_length) <= cfg.window_bp
                    for j in cluster
                )
                if not close:
                    continue
                trial = cluster + [idx]
                trial_af = [afs[j] for j in trial]
                if max(trial_af) - min(trial_af) <= cfg.af_span_max:
                    cluster.append(idx)
                    assigned.add(idx)
                    changed = True
                    break
    return [sorted(cluster, key=lambda i: positions[i]) for cluster in expanded]


def minimal_circular_span(values: Iterable[int], length: int) -> tuple[int, int, int, bool]:
    pos = sorted(set(values))
    if not pos:
        return 0, 0, 0, False
    if len(pos) == 1:
        return pos[0], pos[0], 0, False
    gaps = []
    for a, b in zip(pos, pos[1:]):
        gaps.append((b - a, a, b))
    gaps.append(((pos[0] + length) - pos[-1], pos[-1], pos[0] + length))
    largest_gap, gap_start, gap_end = max(gaps)
    start = ((gap_end - 1) % length) + 1
    end = gap_start
    span = length - largest_gap
    wraps = start > end
    return start, end, span, wraps


def summarize_cluster(indices: Sequence[int], rows: Sequence[dict], cfg: ClusterConfig) -> dict:
    positions = [int(rows[i]["source_pos"]) for i in indices]
    afs = [float(rows[i]["source_af"]) for i in indices]
    start, end, span, wraps = minimal_circular_span(positions, cfg.mt_length)
    return {
        "cluster_start": start,
        "cluster_end": end,
        "cluster_span_bp": span,
        "wraps_origin": str(wraps).lower(),
        "n_variants": len(indices),
        "median_af": median(afs),
        "af_min": min(afs),
        "af_max": max(afs),
        "af_span": max(afs) - min(afs),
    }
