"""Shared helpers for mitochondrial reference identity and circular anchors."""
from __future__ import annotations

import hashlib
import re
from pathlib import Path
from typing import Optional

IUPAC_DNA = set("ACGTRYSWKMBDHVN")
UNAMBIGUOUS_DNA = set("ACGT")
AMBIGUOUS_DNA = IUPAC_DNA - UNAMBIGUOUS_DNA


def normalize_sequence(seq: str) -> str:
    return "".join(str(seq).split()).upper()


def validate_iupac_sequence(seq: str) -> str:
    """Normalize and validate a sequence against the standard IUPAC alphabet."""
    seq = normalize_sequence(seq)
    bad = set(seq) - IUPAC_DNA
    if bad:
        raise ValueError(f"Unexpected FASTA bases: {''.join(sorted(bad))}")
    return seq


def mask_ambiguity_for_alignment(seq: str) -> str:
    """Replace IUPAC ambiguity codes with N without changing sequence length."""
    seq = validate_iupac_sequence(seq)
    return "".join(base if base in UNAMBIGUOUS_DNA else "N" for base in seq)


def sequence_sha256(seq: str) -> str:
    return hashlib.sha256(validate_iupac_sequence(seq).encode()).hexdigest()


def safe_token(value: object) -> str:
    text = str(value or "unknown").strip() or "unknown"
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text)
    return text.strip("_") or "unknown"


def derive_reference_id(species: Optional[str], species_fasta: Path | str, seq_hash: str) -> str:
    return f"{safe_token(species)}__{safe_token(Path(species_fasta).name)}__{seq_hash[:12]}"


def resolve_reference_id(row: dict, species_fasta: Path | str, seq_hash: str) -> str:
    rid = (row.get("reference_id") or "").strip()
    return rid if rid else derive_reference_id(row.get("species") or row.get("sample"), species_fasta, seq_hash)


def rotate_sequence(seq: str, anchor_1based: int) -> str:
    seq = normalize_sequence(seq)
    if anchor_1based < 1 or anchor_1based > len(seq):
        raise ValueError(f"Rotation anchor {anchor_1based} outside 1..{len(seq)}")
    i = anchor_1based - 1
    return seq[i:] + seq[:i]


def rotated_to_original(pos: int, anchor: int, length: int) -> int:
    return ((pos + anchor - 2) % length) + 1


def original_to_rotated(pos: int, anchor: int, length: int) -> int:
    return ((pos - anchor) % length) + 1
