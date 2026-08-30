#!/usr/bin/env python3
"""Re-check deferred sample eligibility immediately before a worker starts.

The submit-time manifest deliberately plans from current sample/reference inventory,
not from possibly stale downstream maps. This helper applies the latest runtime
maps, current collection status, and current liftover QC, then removes a scheduled
sample's prior managed annotation outputs before recomputation.
"""
from __future__ import annotations

import argparse
import csv
import os
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
    return {row[column].strip() for row in rows[start:] if len(row) > column and row[column].strip()}


def exists(path: Path) -> bool:
    return path.is_file() or Path(str(path) + ".gz").is_file()


def formatted(directory: object, pattern: object, sample: str) -> Path:
    return resolve(directory) / str(pattern).format(sample=sample)


def collection_status(sample: str, cfg: dict) -> str | None:
    sec = cfg.get("collect_variant_calling") or {}
    out = resolve(sec.get("outdir", "results/qc/collected_variant_calling_results"))
    report = out / "reports/variant_calling_collection_summary.tsv"
    if not report.is_file():
        return None
    with report.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            name = (row.get("sample") or row.get("Sample") or "").strip()
            if name == sample:
                return (row.get("status") or row.get("collection_status") or "").strip()
    return "MISSING_SAMPLE"


def liftover_status(sample: str, cfg: dict) -> str | None:
    section = cfg.get("coordinate_liftover") or {}
    paths = section.get("paths") or {}
    report = resolve(paths.get("output_dir", "results/qc/coordinate_liftover")) / "reports" / f"{sample}.coordinate_liftover_qc.tsv"
    if not report.is_file():
        return None
    with report.open(newline="", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) >= 2 and row[0].strip() == "status":
                return row[1].strip()
    return "MISSING_STATUS"


def managed_outputs(step: str, sample: str, cfg: dict) -> list[Path]:
    if step == "codon_match":
        sec = cfg.get("codon_match") or {}; p, s = sec.get("paths") or {}, sec.get("settings") or {}
        return [
            resolve(p.get("output_dir", "results/qc/codon_match")) / "vcf_codon" / f"{sample}{s.get('output_suffix', '.lifted.codon.vcf')}",
            resolve(p.get("reports_dir", "results/qc/codon_match/reports")) / f"{sample}.codon_match_summary.tsv",
        ]
    if step == "trna_match":
        sec = cfg.get("trna_match") or {}; p, s = sec.get("paths") or {}, sec.get("settings") or {}
        out = resolve(p.get("output_dir", "results/qc/trna_match")) / "vcf_trna"
        return [
            out / f"{sample}{s.get('output_suffix', '.lifted.codon.trna.vcf')}",
            out / f"{sample}.lifted.trna.vcf",
            resolve(p.get("reports_dir", "results/qc/trna_match/reports")) / f"{sample}.trna_match_summary.tsv",
        ]
    if step == "rrna_match":
        sec = cfg.get("rrna_match") or {}; p, s = sec.get("paths") or {}, sec.get("settings") or {}
        return [
            resolve(p.get("output_dir", "results/qc/rrna_match")) / "vcf_rrna" / f"{sample}{s.get('output_suffix', '.lifted.codon.trna.rrna.vcf')}",
            resolve(p.get("reports_dir", "results/qc/rrna_match/reports")) / f"{sample}.rrna_match_summary.tsv",
        ]
    return []


def cleanup(step: str, sample: str, cfg: dict) -> list[str]:
    removed = []
    for path in managed_outputs(step, sample, cfg):
        for candidate in (path, Path(str(path) + ".gz"), Path(str(path) + ".tbi")):
            if candidate.exists() or candidate.is_symlink():
                candidate.unlink(missing_ok=True)
                removed.append(str(candidate))
    return removed


def rrna_inputs(sample: str, cfg: dict) -> tuple[Path, Path, Path, Path]:
    section = cfg.get("rrna_match") or {}; paths, settings = section.get("paths") or {}, section.get("settings") or {}
    primary = formatted(paths.get("input_vcf_dir", ""), settings.get("input_vcf_pattern", "{sample}.lifted.codon.trna.vcf"), sample)
    trna = formatted(paths.get("fallback_trna_vcf_dir", paths.get("input_vcf_dir", "")), settings.get("fallback_trna_vcf_pattern", "{sample}.lifted.trna.vcf"), sample)
    codon = formatted(paths.get("fallback_codon_vcf_dir", ""), settings.get("fallback_codon_vcf_pattern", "{sample}.lifted.codon.vcf"), sample)
    raw = formatted(paths.get("fallback_raw_vcf_dir", ""), settings.get("fallback_raw_vcf_pattern", "{sample}.lifted.raw.vcf"), sample)
    return primary, trna, codon, raw


def prepare_rrna_trna_alias(sample: str, cfg: dict) -> str:
    """Expose raw-fallback tRNA output under the historical rRNA primary name."""
    primary, trna, _codon, _raw = rrna_inputs(sample, cfg)
    if exists(primary) or not exists(trna):
        return ""
    primary.parent.mkdir(parents=True, exist_ok=True)
    try:
        os.symlink(trna.resolve(), primary)
    except FileExistsError:
        pass
    return str(primary)


def decision(step: str, sample: str, cfg: dict) -> tuple[bool, str]:
    if step == "coordinate_liftover":
        return True, "liftover_handles_missing_inputs"

    current_collection = collection_status(sample, cfg)
    if current_collection is None:
        return False, "collection_summary_not_available"
    if current_collection != "OK":
        return False, "collection_status_" + (current_collection or "blank").lower()

    current_liftover = liftover_status(sample, cfg)
    if current_liftover is None:
        return False, "liftover_qc_not_available"
    if current_liftover != "completed":
        return False, "liftover_status_" + (current_liftover or "blank").lower()

    if step == "codon_match":
        section = cfg.get("codon_match") or {}; paths, settings = section.get("paths") or {}, section.get("settings") or {}
        mapping = paths.get("sample_reference_map")
        if not mapping or not resolve(mapping).is_file(): return False, "production_codon_map_not_available"
        if sample not in table_samples(mapping): return False, "sample_not_in_pass_production_codon_map"
        inp = formatted(paths.get("input_vcf_dir", ""), settings.get("input_vcf_pattern", "{sample}.lifted.raw.vcf"), sample)
        if not exists(inp): return False, "liftover_vcf_not_available"
        return True, "eligible"

    if step == "trna_match":
        section = cfg.get("trna_match") or {}; paths, settings = section.get("paths") or {}, section.get("settings") or {}
        mapping = paths.get("sample_reference_map")
        if mapping:
            if not resolve(mapping).is_file(): return False, "trna_reference_map_not_available"
            if sample not in table_samples(mapping): return False, "sample_not_in_trna_reference_map"
        primary = formatted(paths.get("input_vcf_dir", ""), settings.get("input_vcf_pattern", "{sample}.lifted.codon.vcf"), sample)
        fallback = formatted(paths.get("fallback_input_vcf_dir", ""), settings.get("fallback_input_vcf_pattern", "{sample}.lifted.raw.vcf"), sample)
        if not (exists(primary) or exists(fallback)): return False, "no_trna_input_vcf_available"
        return True, "eligible"

    if step == "rrna_match":
        section = cfg.get("rrna_match") or {}; paths = section.get("paths") or {}
        mapping = paths.get("sample_reference_map")
        if mapping:
            if not resolve(mapping).is_file(): return False, "rrna_reference_map_not_available"
            if sample not in table_samples(mapping): return False, "sample_not_in_rrna_reference_map"
        if not any(exists(path) for path in rrna_inputs(sample, cfg)): return False, "no_rrna_input_vcf_available"
        return True, "eligible"

    return True, "not_runtime_filtered"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", choices=("coordinate_liftover", "codon_match", "trna_match", "rrna_match"))
    parser.add_argument("sample")
    parser.add_argument("config", type=Path)
    args = parser.parse_args(); cfg = read_simple_yaml(args.config)
    eligible, reason = decision(args.step, args.sample, cfg)
    removed = cleanup(args.step, args.sample, cfg) if args.step in {"codon_match", "trna_match", "rrna_match"} else []
    alias = prepare_rrna_trna_alias(args.sample, cfg) if eligible and args.step == "rrna_match" else ""
    print(f"ELIGIBLE={1 if eligible else 0}")
    print(f"REASON={reason}")
    print("CLEANED=" + ";".join(removed))
    print(f"PREPARED_INPUT_ALIAS={alias}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
