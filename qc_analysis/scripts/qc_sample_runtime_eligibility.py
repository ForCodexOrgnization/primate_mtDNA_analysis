#!/usr/bin/env python3
"""Re-check sample eligibility when a deferred Slurm array worker actually starts.

Array manifests for codon/tRNA/rRNA are intentionally allowed to be planned before
all upstream products exist.  This helper turns later upstream exclusions (for
example a MITOS2 reference that fails production QC) into a successful per-sample
skip rather than an array failure.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from qc_analysis.lib.simple_yaml import read_simple_yaml


def resolve(value: object) -> Path:
    path = Path(str(value)).expanduser()
    return path if path.is_absolute() else ROOT / path


def table_samples(path: object) -> set[str]:
    if not path:
        return set()
    p = resolve(path)
    if not p.is_file():
        return set()
    with p.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows:
        return set()
    header = [value.strip().lower() for value in rows[0]]
    names = ("sample", "sample_id", "name")
    headered = any(name in header for name in names)
    column = next((header.index(name) for name in names if name in header), 0)
    start = 1 if headered else 0
    return {
        row[column].strip()
        for row in rows[start:]
        if len(row) > column and row[column].strip()
    }


def exists(path: Path) -> bool:
    return path.is_file() or Path(str(path) + ".gz").is_file()


def formatted(directory: object, pattern: object, sample: str) -> Path:
    return resolve(directory) / str(pattern).format(sample=sample)


def decision(step: str, sample: str, cfg: dict) -> tuple[bool, str]:
    if step == "coordinate_liftover":
        # Liftover itself already records skipped_missing_file and returns success,
        # which is more informative than silently dropping the sample here.
        return True, "liftover_handles_missing_inputs"

    if step == "codon_match":
        section = cfg.get("codon_match") or {}
        paths, settings = section.get("paths") or {}, section.get("settings") or {}
        mapping = paths.get("sample_reference_map")
        if not mapping or not resolve(mapping).is_file():
            return False, "production_codon_map_not_available"
        if sample not in table_samples(mapping):
            return False, "sample_not_in_pass_production_codon_map"
        inp = formatted(paths.get("input_vcf_dir", ""), settings.get("input_vcf_pattern", "{sample}.lifted.raw.vcf"), sample)
        if not exists(inp):
            return False, "liftover_vcf_not_available"
        return True, "eligible"

    if step == "trna_match":
        section = cfg.get("trna_match") or {}
        paths, settings = section.get("paths") or {}, section.get("settings") or {}
        mapping = paths.get("sample_reference_map")
        if mapping:
            if not resolve(mapping).is_file():
                return False, "trna_reference_map_not_available"
            if sample not in table_samples(mapping):
                return False, "sample_not_in_trna_reference_map"
        primary = formatted(paths.get("input_vcf_dir", ""), settings.get("input_vcf_pattern", "{sample}.lifted.codon.vcf"), sample)
        fallback = formatted(paths.get("fallback_input_vcf_dir", ""), settings.get("fallback_input_vcf_pattern", "{sample}.lifted.raw.vcf"), sample)
        if not (exists(primary) or exists(fallback)):
            return False, "no_trna_input_vcf_available"
        return True, "eligible"

    if step == "rrna_match":
        section = cfg.get("rrna_match") or {}
        paths, settings = section.get("paths") or {}, section.get("settings") or {}
        mapping = paths.get("sample_reference_map")
        if mapping:
            if not resolve(mapping).is_file():
                return False, "rrna_reference_map_not_available"
            if sample not in table_samples(mapping):
                return False, "sample_not_in_rrna_reference_map"
        choices = [
            formatted(paths.get("input_vcf_dir", ""), settings.get("input_vcf_pattern", "{sample}.lifted.codon.trna.vcf"), sample),
            formatted(paths.get("fallback_codon_vcf_dir", ""), settings.get("fallback_codon_vcf_pattern", "{sample}.lifted.codon.vcf"), sample),
            formatted(paths.get("fallback_raw_vcf_dir", ""), settings.get("fallback_raw_vcf_pattern", "{sample}.lifted.raw.vcf"), sample),
        ]
        if not any(exists(path) for path in choices):
            return False, "no_rrna_input_vcf_available"
        return True, "eligible"

    return True, "not_runtime_filtered"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", choices=("coordinate_liftover", "codon_match", "trna_match", "rrna_match"))
    parser.add_argument("sample")
    parser.add_argument("config", type=Path)
    args = parser.parse_args()
    cfg = read_simple_yaml(args.config)
    eligible, reason = decision(args.step, args.sample, cfg)
    print(f"ELIGIBLE={1 if eligible else 0}")
    print(f"REASON={reason}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
