#!/usr/bin/env python3
"""Discover conserved mtDNA coordinate-liftover anchors across primates.

The script scans human chrM windows, locally aligns each window against every
available species chrM FASTA, scores conservation/uniqueness/complexity, and
writes candidate, selected-anchor, and per-sample/species anchor tables.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

DNA = set("ACGTNacgtn")
FASTA_EXTENSIONS = (".fa", ".fasta", ".fna", ".fas")


@dataclass
class FastaRecord:
    name: str
    seq: str


@dataclass
class LocalHit:
    score: int
    start: int
    end: int
    strand: str
    identity: float
    aligned_length: int
    gaps: int


def read_fasta(path: Path, target: Optional[str] = None) -> FastaRecord:
    records: List[FastaRecord] = []
    name: Optional[str] = None
    chunks: List[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append(FastaRecord(name, "".join(chunks).upper()))
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        records.append(FastaRecord(name, "".join(chunks).upper()))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    if target:
        matches = [r for r in records if r.name == target]
        if len(matches) != 1:
            raise ValueError(f"Expected one target sequence {target!r} in {path}, found {len(matches)}")
        rec = matches[0]
    elif len(records) == 1:
        rec = records[0]
    else:
        mt = [r for r in records if r.name.lower() in {"chrm", "mt", "m", "mitochondrion", "mitochondrial"}]
        if len(mt) != 1:
            raise ValueError(f"{path} contains multiple records; set a target sequence")
        rec = mt[0]
    bad = set(rec.seq) - DNA
    if bad:
        raise ValueError(f"Unexpected FASTA bases in {path}: {''.join(sorted(bad))}")
    return rec


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


def local_align(query: str, target: str, strand: str, match: int = 2, mismatch: int = -1, gap: int = -3) -> LocalHit:
    """Return the best Smith-Waterman local alignment of query to target."""
    q = query.upper()
    t = target.upper()
    prev = [0] * (len(t) + 1)
    best = (0, 0, 0)
    pointer: Dict[Tuple[int, int], Tuple[int, int]] = {}
    for i, qb in enumerate(q, start=1):
        curr = [0] * (len(t) + 1)
        for j, tb in enumerate(t, start=1):
            diag = prev[j - 1] + (match if qb == tb and qb != "N" else mismatch)
            up = prev[j] + gap
            left = curr[j - 1] + gap
            score = max(0, diag, up, left)
            curr[j] = score
            if score == 0:
                continue
            if score == diag:
                pointer[(i, j)] = (i - 1, j - 1)
            elif score == up:
                pointer[(i, j)] = (i - 1, j)
            else:
                pointer[(i, j)] = (i, j - 1)
            if score > best[0]:
                best = (score, i, j)
        prev = curr
    score, i, j = best
    end = j
    matches = aligned = gaps = 0
    while (i, j) in pointer:
        pi, pj = pointer[(i, j)]
        if pi == i - 1 and pj == j - 1:
            aligned += 1
            matches += int(q[i - 1] == t[j - 1] and q[i - 1] != "N")
        else:
            aligned += 1
            gaps += 1
        i, j = pi, pj
        if i == 0 or j == 0:
            break
    start = j + 1
    identity = matches / aligned if aligned else 0.0
    return LocalHit(score, start, end, strand, identity, aligned, gaps)


def best_local_hits(query: str, species_seq: str) -> Tuple[LocalHit, LocalHit]:
    doubled = species_seq.upper() + species_seq.upper()[: max(0, len(query) - 1)]
    hits = [local_align(query, doubled, "+"), local_align(reverse_complement(query), doubled, "-")]
    hits.sort(key=lambda h: h.score, reverse=True)
    best = hits[0]
    # Mask a small neighborhood around the best hit and realign to estimate uniqueness.
    masked = list(doubled)
    flank = len(query)
    for idx in range(max(0, best.start - flank - 1), min(len(masked), best.end + flank)):
        masked[idx] = "N"
    second_candidates = [hits[1], local_align(query, "".join(masked), "+"), local_align(reverse_complement(query), "".join(masked), "-")]
    second_candidates.sort(key=lambda h: h.score, reverse=True)
    best.start = ((best.start - 1) % len(species_seq)) + 1 if species_seq else 0
    best.end = ((best.end - 1) % len(species_seq)) + 1 if species_seq else 0
    return best, second_candidates[0]


def low_complexity_score(seq: str) -> float:
    counts = Counter(seq.upper())
    entropy = -sum((n / len(seq)) * math.log2(n / len(seq)) for n in counts.values() if n) if seq else 0.0
    # Combine low entropy and dominant-base fraction as a simple repeat/composition proxy.
    dominant = max(counts.values()) / len(seq) if seq else 1.0
    return max(0.0, min(1.0, (2.0 - entropy) / 2.0)) * 0.7 + dominant * 0.3


def in_variable_region(start: int, end: int, human_len: int) -> bool:
    # Human chrM D-loop/control-region coordinates are comparatively variable.
    variable = [(16024, human_len), (1, 576)]
    return any(start <= b and end >= a for a, b in variable)


def read_metadata(path: Path) -> List[dict]:
    with path.open(newline="") as handle:
        sample_rows = list(csv.DictReader(handle, delimiter="\t"))
    if not sample_rows:
        raise ValueError(f"No metadata rows found in {path}")
    return sample_rows


def find_species_fasta(species: str, fasta_dir: Path) -> Optional[Path]:
    for ext in FASTA_EXTENSIONS:
        path = fasta_dir / f"{species}{ext}"
        if path.exists():
            return path
    matches = sorted(p for p in fasta_dir.iterdir() if p.is_file() and p.suffix.lower() in FASTA_EXTENSIONS and p.stem == species)
    return matches[0] if matches else None


def main(argv: Optional[Iterable[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--human-fasta", required=True, type=Path)
    ap.add_argument("--species-fasta-dir", required=True, type=Path)
    ap.add_argument("--metadata", required=True, type=Path, help="TSV with sample and species columns")
    ap.add_argument("--output-dir", required=True, type=Path)
    ap.add_argument("--window-length", type=int, default=150)
    ap.add_argument("--step-size", type=int, default=20)
    ap.add_argument("--min-identity", type=float, default=0.65)
    ap.add_argument("--min-aligned-fraction", type=float, default=0.75)
    args = ap.parse_args(argv)

    human = read_fasta(args.human_fasta)
    rows = read_metadata(args.metadata)
    species_to_samples: Dict[str, List[str]] = {}
    for row in rows:
        sample = row.get("sample", "").strip()
        species = (row.get("species") or sample).strip()
        if sample and species:
            species_to_samples.setdefault(species, []).append(sample)
    species_records: Dict[str, FastaRecord] = {}
    for species in species_to_samples:
        path = find_species_fasta(species, args.species_fasta_dir)
        if path:
            species_records[species] = read_fasta(path)
    if not species_records:
        raise ValueError("No species FASTAs could be resolved from metadata")

    outdir = args.output_dir / "anchor_discovery"
    outdir.mkdir(parents=True, exist_ok=True)
    candidate_rows: List[dict] = []
    per_window_hits: Dict[Tuple[int, int], Dict[str, Tuple[LocalHit, LocalHit]]] = {}
    for start0 in range(0, len(human.seq) - args.window_length + 1, args.step_size):
        start = start0 + 1
        end = start0 + args.window_length
        query = human.seq[start0:end]
        hits: Dict[str, Tuple[LocalHit, LocalHit]] = {}
        identities: List[float] = []
        lengths: List[int] = []
        uniqueness: List[float] = []
        gap_rates: List[float] = []
        found = 0
        for species, rec in species_records.items():
            best, second = best_local_hits(query, rec.seq)
            hits[species] = (best, second)
            found_here = best.identity >= args.min_identity and best.aligned_length >= args.window_length * args.min_aligned_fraction
            found += int(found_here)
            if found_here:
                identities.append(best.identity)
                lengths.append(best.aligned_length)
                uniqueness.append(best.score / second.score if second.score > 0 else float(best.score))
                gap_rates.append(best.gaps / best.aligned_length if best.aligned_length else 1.0)
        frac = found / len(species_records)
        mean_id = sum(identities) / len(identities) if identities else 0.0
        min_id = min(identities) if identities else 0.0
        mean_len = sum(lengths) / len(lengths) if lengths else 0.0
        uniq = sum(uniqueness) / len(uniqueness) if uniqueness else 0.0
        low_complex = low_complexity_score(query)
        gap_pen = sum(gap_rates) / len(gap_rates) if gap_rates else 1.0
        variable_penalty = 1.0 if in_variable_region(start, end, len(human.seq)) else 0.0
        composite = frac * 100 + mean_id * 35 + min_id * 20 + min(uniq, 3.0) * 8 - low_complex * 15 - gap_pen * 15 - variable_penalty * 10
        per_window_hits[(start, end)] = hits
        candidate_rows.append({
            "human_anchor_start": start, "human_anchor_end": end, "human_anchor_pos": start,
            "species_found_count": found, "species_found_fraction": f"{frac:.6f}",
            "mean_identity": f"{mean_id:.6f}", "min_identity": f"{min_id:.6f}",
            "mean_aligned_length": f"{mean_len:.2f}", "uniqueness_score": f"{uniq:.6f}",
            "low_complexity_score": f"{low_complex:.6f}", "gap_penalty": f"{gap_pen:.6f}",
            "variable_region_penalty": f"{variable_penalty:.1f}", "composite_score": f"{composite:.6f}",
        })
    candidate_rows.sort(key=lambda r: float(r["composite_score"]), reverse=True)
    fieldnames = list(candidate_rows[0].keys())
    with (outdir / "human_anchor_candidates.tsv").open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader(); writer.writerows(candidate_rows)
    selected = candidate_rows[0]
    with (outdir / "selected_human_anchor.tsv").open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader(); writer.writerow(selected)

    start, end = int(selected["human_anchor_start"]), int(selected["human_anchor_end"])
    hits = per_window_hits[(start, end)]
    position_fields = ["sample", "species", "human_anchor_pos", "human_anchor_start", "human_anchor_end", "species_anchor_pos", "strand", "anchor_identity", "anchor_aligned_length", "best_score", "second_best_score", "uniqueness_ratio", "anchor_status"]
    with (outdir / "species_anchor_positions.tsv").open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=position_fields, delimiter="\t")
        writer.writeheader()
        for species, samples in sorted(species_to_samples.items()):
            best, second = hits.get(species, (LocalHit(0, 0, 0, "+", 0.0, 0, 0), LocalHit(0, 0, 0, "+", 0.0, 0, 0)))
            ratio = best.score / second.score if second.score else float(best.score)
            status = "selected" if species in species_records and best.identity >= args.min_identity else "low_confidence"
            if species not in species_records:
                status = "missing_species_fasta"
            for sample in samples:
                writer.writerow({"sample": sample, "species": species, "human_anchor_pos": selected["human_anchor_pos"], "human_anchor_start": start, "human_anchor_end": end, "species_anchor_pos": best.start, "strand": best.strand, "anchor_identity": f"{best.identity:.6f}", "anchor_aligned_length": best.aligned_length, "best_score": best.score, "second_best_score": second.score, "uniqueness_ratio": f"{ratio:.6f}", "anchor_status": status})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
