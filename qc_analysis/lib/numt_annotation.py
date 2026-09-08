"""Sample- and species-level NUMT interval annotation helpers."""
from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

BESTHIT_COLUMNS = [
    "nuclear_chrom", "nuclear_start", "nuclear_end", "n_blocks", "sum_nreads",
    "max_nreads", "mean_nreads", "mean_MAPQ", "locus_id", "chrM_hit", "pident",
    "align_length", "mismatch", "gapopen", "qstart", "qend", "chrM_start",
    "chrM_end", "evalue", "bitscore",
]


@dataclass(frozen=True)
class NumtInterval:
    sample: str
    species: str
    reference_key: str
    nuclear_chrom: str
    nuclear_start: int
    nuclear_end: int
    chrm_start: int
    chrm_end: int
    tier: str
    locus_id: str = ""
    pident: str = ""
    align_length: str = ""
    sum_nreads: str = ""
    mean_mapq: str = ""


def _int(value) -> int:
    return int(float(str(value)))


def read_highconf_bed(path: Path) -> list[tuple[str, int, int]]:
    if not path.is_file():
        return []
    out = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 3:
                continue
            try:
                out.append((fields[0], _int(fields[1]), _int(fields[2])))
            except ValueError:
                continue
    return out


def nuclear_overlaps_highconf(chrom: str, start: int, end: int, highconf: Iterable[tuple[str, int, int]]) -> bool:
    lo, hi = sorted((start, end))
    for hchrom, hstart, hend in highconf:
        hlo, hhi = sorted((hstart, hend))
        if chrom == hchrom and max(lo, hlo) <= min(hi, hhi):
            return True
    return False


def read_besthit(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    with path.open(newline="") as handle:
        rows = [row for row in csv.reader(handle, delimiter="\t") if row]
    if not rows:
        return []
    first = rows[0]
    looks_headered = any(value in set(BESTHIT_COLUMNS) for value in first)
    if looks_headered:
        with path.open(newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))
    out = []
    for row in rows:
        if len(row) < len(BESTHIT_COLUMNS):
            continue
        out.append(dict(zip(BESTHIT_COLUMNS, row[: len(BESTHIT_COLUMNS)])))
    return out


def load_sample_numts(
    sample: str,
    species: str,
    reference_key: str,
    besthit_path: Path,
    highconf_path: Path,
) -> list[NumtInterval]:
    highconf = read_highconf_bed(highconf_path)
    out: list[NumtInterval] = []
    for row in read_besthit(besthit_path):
        try:
            nuclear_start = _int(row["nuclear_start"])
            nuclear_end = _int(row["nuclear_end"])
            chrm_start = _int(row["chrM_start"])
            chrm_end = _int(row["chrM_end"])
        except (KeyError, ValueError, TypeError):
            continue
        tier = "HIGH_CONF_NUMT" if nuclear_overlaps_highconf(
            row.get("nuclear_chrom", ""), nuclear_start, nuclear_end, highconf
        ) else "BESTHIT_ONLY_NUMT"
        out.append(
            NumtInterval(
                sample=sample,
                species=species,
                reference_key=reference_key,
                nuclear_chrom=row.get("nuclear_chrom", ""),
                nuclear_start=min(nuclear_start, nuclear_end),
                nuclear_end=max(nuclear_start, nuclear_end),
                chrm_start=min(chrm_start, chrm_end),
                chrm_end=max(chrm_start, chrm_end),
                tier=tier,
                locus_id=row.get("locus_id", ""),
                pident=row.get("pident", ""),
                align_length=row.get("align_length", ""),
                sum_nreads=row.get("sum_nreads", ""),
                mean_mapq=row.get("mean_MAPQ", ""),
            )
        )
    return out


def position_in_interval(pos: int, start: int, end: int) -> bool:
    return start <= pos <= end


def best_cluster_overlap(cluster_rows: list[dict], intervals: Iterable[NumtInterval]) -> dict | None:
    best = None
    cluster_positions = [int(row["source_pos"]) for row in cluster_rows]
    n = len(cluster_positions)
    for interval in intervals:
        overlap_positions = [pos for pos in cluster_positions if position_in_interval(pos, interval.chrm_start, interval.chrm_end)]
        overlap_n = len(overlap_positions)
        fraction = overlap_n / n if n else 0.0
        rank = (
            overlap_n,
            fraction,
            1 if interval.tier == "HIGH_CONF_NUMT" else 0,
            float(interval.pident) if str(interval.pident).replace(".", "", 1).isdigit() else -1.0,
            int(float(interval.align_length)) if str(interval.align_length).replace(".", "", 1).isdigit() else -1,
        )
        if best is None or rank > best[0]:
            best = (rank, interval, overlap_n, fraction)
    if best is None:
        return None
    _, interval, overlap_n, fraction = best
    return {
        "interval": interval,
        "overlap_n": overlap_n,
        "overlap_fraction": fraction,
    }


def merge_species_intervals(intervals: Iterable[NumtInterval], merge_gap_bp: int = 0) -> list[dict]:
    grouped: dict[tuple[str, str, str], list[NumtInterval]] = defaultdict(list)
    for interval in intervals:
        grouped[(interval.species, interval.reference_key, interval.tier)].append(interval)
    merged = []
    for (species, reference_key, tier), values in sorted(grouped.items()):
        values = sorted(values, key=lambda x: (x.chrm_start, x.chrm_end, x.sample))
        current = None
        for value in values:
            if current is None or value.chrm_start > current["chrm_end"] + merge_gap_bp:
                if current is not None:
                    merged.append(current)
                current = {
                    "species": species,
                    "reference_key": reference_key,
                    "tier": tier,
                    "chrm_start": value.chrm_start,
                    "chrm_end": value.chrm_end,
                    "support_samples": {value.sample},
                    "nuclear_loci": {(value.nuclear_chrom, value.nuclear_start, value.nuclear_end)},
                }
            else:
                current["chrm_end"] = max(current["chrm_end"], value.chrm_end)
                current["support_samples"].add(value.sample)
                current["nuclear_loci"].add((value.nuclear_chrom, value.nuclear_start, value.nuclear_end))
        if current is not None:
            merged.append(current)
    return merged


def best_species_overlap(
    cluster_rows: list[dict],
    sample: str,
    species: str,
    reference_key: str,
    species_intervals: Iterable[dict],
) -> dict | None:
    positions = [int(row["source_pos"]) for row in cluster_rows]
    n = len(positions)
    best = None
    for interval in species_intervals:
        if interval["species"] != species or interval["reference_key"] != reference_key:
            continue
        other_support = sorted(set(interval["support_samples"]) - {sample})
        if not other_support:
            continue
        overlap_n = sum(interval["chrm_start"] <= pos <= interval["chrm_end"] for pos in positions)
        fraction = overlap_n / n if n else 0.0
        rank = (overlap_n, fraction, len(other_support), 1 if interval["tier"] == "HIGH_CONF_NUMT" else 0)
        if best is None or rank > best[0]:
            best = (rank, interval, overlap_n, fraction, other_support)
    if best is None:
        return None
    _, interval, overlap_n, fraction, other_support = best
    return {
        "interval": interval,
        "overlap_n": overlap_n,
        "overlap_fraction": fraction,
        "other_support_samples": other_support,
    }
