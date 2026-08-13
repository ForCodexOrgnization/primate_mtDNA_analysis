#!/usr/bin/env python3
"""Normalize a complete PhyloTree/HaploGrep rCRS marker list for screening."""
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

MARKER_RE = re.compile(r"^([0-9]+)([ACGT])(!?)$")


def prepare(source: Path, output: Path, summary: Path | None = None) -> dict[str, int]:
    markers: dict[tuple[int, str], dict[str, object]] = {}
    qc = {"raw_marker_count": 0, "parsed_simple_snvs": 0,
          "back_mutation_markers": 0, "excluded_complex_markers": 0,
          "duplicate_pos_alt_rows_removed": 0}
    # Complete exports occur both one-marker-per-line and in delimited columns.
    tokens = []
    for line in source.read_text(encoding="utf-8-sig").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        tokens.extend(x for x in re.split(r"[\t,; ]+", line.strip()) if x)
    for marker in tokens:
        if marker.lower() in {"marker", "mutation", "polymorphism"}:
            continue
        qc["raw_marker_count"] += 1
        match = MARKER_RE.fullmatch(marker)
        if not match:
            qc["excluded_complex_markers"] += 1
            continue
        qc["parsed_simple_snvs"] += 1
        pos, alt, back = int(match[1]), match[2], bool(match[3])
        key = (pos, alt)
        if key in markers:
            qc["duplicate_pos_alt_rows_removed"] += 1
            markers[key]["is_back_mutation"] = bool(markers[key]["is_back_mutation"]) or back
        else:
            markers[key] = {"marker": marker, "pos": pos, "alt": alt, "is_back_mutation": back}
    qc["back_mutation_markers"] = sum(bool(x["is_back_mutation"]) for x in markers.values())
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("marker", "pos", "alt", "is_back_mutation"), delimiter="\t")
        writer.writeheader()
        for item in sorted(markers.values(), key=lambda x: (int(x["pos"]), str(x["alt"]))):
            writer.writerow({**item, "is_back_mutation": str(item["is_back_mutation"]).lower()})
    if summary:
        summary.parent.mkdir(parents=True, exist_ok=True)
        summary.write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8")
    return qc


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()
    if not args.input.is_file():
        parser.error(f"raw marker file does not exist: {args.input}")
    print(json.dumps(prepare(args.input, args.output, args.summary), sort_keys=True))


if __name__ == "__main__":
    main()
