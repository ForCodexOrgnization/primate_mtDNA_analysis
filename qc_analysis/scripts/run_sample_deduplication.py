#!/usr/bin/env python3
"""Build the canonical, deduplicated QC sample cohort from enriched metadata.

The input metadata contains one row per input accession. Duplicate biological
samples are represented by one ``canonical_id_row`` plus one
``archive_alias_row``. Production QC keeps ``unique_input`` and
``canonical_id_row`` rows and excludes archive aliases before any downstream
sample-level analysis.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from qc_analysis.lib.simple_yaml import read_simple_yaml


def resolve(value: str | Path) -> Path:
    p = Path(str(value)).expanduser()
    return p if p.is_absolute() else ROOT / p


def split_csv(value, default: str) -> set[str]:
    if value is None:
        value = default
    if isinstance(value, str):
        return {x.strip() for x in value.split(",") if x.strip()}
    return {str(x).strip() for x in value if str(x).strip()}


def read_rows(path: Path) -> tuple[list[dict], list[str]]:
    if not path.is_file():
        raise FileNotFoundError(
            f"sample deduplication metadata not found: {path}. "
            "Set sample_deduplication.metadata in config/qc_preprocessing.yaml."
        )
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise RuntimeError(f"metadata has no header: {path}")
        return list(reader), list(reader.fieldnames)


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def clean(value) -> str:
    value = "" if value is None else str(value).strip()
    return "" if value.upper() == "NA" else value


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, required=True)
    args = ap.parse_args()

    cfg = read_simple_yaml(args.config)
    sec = cfg.get("sample_deduplication") or {}
    if sec.get("enabled", True) is False:
        print("[sample_deduplication] disabled; skipping.")
        return 0

    metadata_path = resolve(sec.get("metadata", "data/metadata/primate_metadata_master_3332_enriched.tsv"))
    output_dir = resolve(sec.get("output_dir", "results/qc/sample_deduplication"))
    reports = output_dir / "reports"
    reports.mkdir(parents=True, exist_ok=True)

    accession_col = str(sec.get("input_accession_column", "input_accession"))
    canonical_col = str(sec.get("canonical_sample_id_column", "canonical_sample_id"))
    species_col = str(sec.get("species_column", "pipeline_species"))
    role_col = str(sec.get("duplicate_role_column", "duplicate_role"))
    duplicate_flag_col = str(sec.get("duplicate_flag_column", "duplicate_flag"))
    group_col = str(sec.get("duplicate_group_column", "duplicate_group"))
    group_size_col = str(sec.get("duplicate_group_size_column", "duplicate_group_size"))
    alias_col = str(sec.get("duplicate_alias_column", "duplicate_alias_accession"))
    keep_roles = split_csv(sec.get("keep_roles"), "unique_input,canonical_id_row")
    exclude_roles = split_csv(sec.get("exclude_roles"), "archive_alias_row")

    rows, fields = read_rows(metadata_path)
    required = {accession_col, canonical_col, species_col, role_col}
    missing = sorted(required - set(fields))
    if missing:
        raise RuntimeError(f"metadata missing required columns {missing}: {metadata_path}")

    seen_accessions: set[str] = set()
    duplicate_accession_rows: list[str] = []
    role_counts = Counter()
    by_group: dict[str, list[dict]] = defaultdict(list)
    kept: list[dict] = []
    excluded: list[dict] = []
    unknown_roles: list[tuple[int, str, str]] = []

    for line_no, row in enumerate(rows, 2):
        accession = clean(row.get(accession_col))
        canonical = clean(row.get(canonical_col))
        species = clean(row.get(species_col))
        role = clean(row.get(role_col))
        group = clean(row.get(group_col))

        if not accession:
            raise RuntimeError(f"empty {accession_col} at metadata line {line_no}")
        if accession in seen_accessions:
            duplicate_accession_rows.append(accession)
        seen_accessions.add(accession)
        if not canonical:
            raise RuntimeError(f"empty {canonical_col} for {accession} at metadata line {line_no}")
        if not species:
            raise RuntimeError(f"empty {species_col} for {accession} at metadata line {line_no}")
        if not role:
            raise RuntimeError(f"empty {role_col} for {accession} at metadata line {line_no}")

        role_counts[role] += 1
        if group:
            by_group[group].append(row)

        base = {
            "sample": accession,
            "species": species,
            "canonical_sample_id": canonical,
            "duplicate_flag": clean(row.get(duplicate_flag_col)) or "False",
            "duplicate_group": group or "NA",
            "duplicate_group_size": clean(row.get(group_size_col)) or "NA",
            "duplicate_role": role,
            "duplicate_alias_accession": clean(row.get(alias_col)) or "NA",
        }

        if role in keep_roles:
            kept.append(base)
        elif role in exclude_roles:
            excluded.append(base)
        else:
            unknown_roles.append((line_no, accession, role))

    if duplicate_accession_rows:
        examples = ",".join(sorted(set(duplicate_accession_rows))[:10])
        raise RuntimeError(f"metadata contains duplicate input_accession rows; examples: {examples}")
    if unknown_roles:
        examples = "; ".join(f"line {n} {s}:{r}" for n, s, r in unknown_roles[:10])
        raise RuntimeError(f"metadata contains unhandled duplicate_role values: {examples}")

    # Each declared duplicate group must have exactly one canonical representative
    # and one or more aliases. This avoids silently choosing an arbitrary accession.
    bad_groups = []
    for group, members in sorted(by_group.items()):
        roles = [clean(r.get(role_col)) for r in members]
        n_canonical = roles.count("canonical_id_row")
        n_alias = roles.count("archive_alias_row")
        if n_canonical != 1 or n_alias < 1:
            bad_groups.append((group, n_canonical, n_alias, len(members)))
    if bad_groups:
        examples = "; ".join(
            f"{g}:canonical={c},alias={a},rows={n}" for g, c, a, n in bad_groups[:10]
        )
        raise RuntimeError(f"invalid duplicate groups in metadata: {examples}")

    kept.sort(key=lambda r: r["sample"])
    excluded.sort(key=lambda r: r["sample"])

    # Guard against accidental many-to-one collapse among rows explicitly marked
    # unique_input. The enriched metadata can contain composite canonical IDs; the
    # duplicate_role/duplicate_group annotations, not canonical ID equality alone,
    # define which accessions are aliases.
    kept_samples = {r["sample"] for r in kept}
    if len(kept_samples) != len(kept):
        raise RuntimeError("deduplicated cohort still contains duplicate sample accessions")

    detail_fields = [
        "sample", "species", "canonical_sample_id", "duplicate_flag",
        "duplicate_group", "duplicate_group_size", "duplicate_role",
        "duplicate_alias_accession",
    ]
    write_tsv(reports / "deduplicated_samples.tsv", kept, detail_fields)
    write_tsv(reports / "excluded_duplicate_aliases.tsv", excluded, detail_fields)
    write_tsv(
        reports / "deduplicated_sample_ref_file.tsv",
        [{"sample": r["sample"], "species": r["species"]} for r in kept],
        ["sample", "species"],
    )

    summary = [{
        "metadata": str(metadata_path),
        "input_rows": len(rows),
        "unique_input_accessions": len(seen_accessions),
        "kept_accessions": len(kept),
        "excluded_alias_accessions": len(excluded),
        "duplicate_groups": len(by_group),
        "unique_input_rows": role_counts.get("unique_input", 0),
        "canonical_id_rows": role_counts.get("canonical_id_row", 0),
        "archive_alias_rows": role_counts.get("archive_alias_row", 0),
    }]
    write_tsv(
        reports / "sample_deduplication_summary.tsv",
        summary,
        [
            "metadata", "input_rows", "unique_input_accessions", "kept_accessions",
            "excluded_alias_accessions", "duplicate_groups", "unique_input_rows",
            "canonical_id_rows", "archive_alias_rows",
        ],
    )

    print(f"[sample_deduplication] metadata={metadata_path}")
    print(f"[sample_deduplication] input_rows={len(rows)}")
    print(f"[sample_deduplication] kept={len(kept)} excluded_aliases={len(excluded)} duplicate_groups={len(by_group)}")
    print(f"[sample_deduplication] sample_ref={reports / 'deduplicated_sample_ref_file.tsv'}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, ValueError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
