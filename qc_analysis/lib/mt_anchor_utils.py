"""Shared helpers for mitochondrial reference identity and circular anchors."""
from __future__ import annotations

import hashlib
import re
from pathlib import Path
from typing import Optional

DNA = set("ACGTNacgtn")


def normalize_sequence(seq: str) -> str:
    return "".join(str(seq).split()).upper()


def sequence_sha256(seq: str) -> str:
    return hashlib.sha256(normalize_sequence(seq).encode()).hexdigest()


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
