#!/usr/bin/env python3
"""Build the deduplicated QC sample cohort from the FINAL enriched metadata.

Archive-level aliases are removed before any downstream QC step. The FINAL
metadata also contains cross-BioSample biological-duplicate candidates. These
are reported for review but are NOT automatically collapsed, because the
metadata labels them candidates rather than confirmed duplicate accessions.
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
            "Set sample_deduplication.metadata in config/sample_deduplication.yaml."
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

    metadata_path = resolve(sec.get("metadata", "data/metadata/primate_metadata_master_3332_FINAL.tsv"))
    output_dir = resolve(sec.get("output_dir", "results/qc/sample_deduplication"))
    reports = output_dir / "reports"
    reports.mkdir(parents=True, exist_ok=True)

    accession_col = str(sec.get("input_accession_column", "input_accession"))
    canonical_col = str(sec.get("canonical_sample_id_column", "canonical_sample_id"))
    species_col = str(sec.get("species_column", "pipeline_species"))

    # FINAL metadata archive-level duplicate columns.
    role_col = str(sec.get("archive_duplicate_role_column", "archive_duplicate_role"))
    duplicate_flag_col = str(sec.get("archive_duplicate_flag_column", "archive_duplicate_flag"))
    group_col = str(sec.get("archive_duplicate_group_column", "archive_duplicate_group"))
    group_size_col = str(sec.get("archive_duplicate_group_size_column", "archive_duplicate_group_size"))
    alias_col = str(sec.get("archive_duplicate_alias_column", "archive_duplicate_alias_accession"))

    # FINAL metadata biological duplicate candidate columns.
    bio_flag_col = str(sec.get("biological_duplicate_candidate_flag_column", "biological_duplicate_candidate_flag"))
    bio_group_col = str(sec.get("biological_duplicate_candidate_group_column", "biological_duplicate_candidate_group"))
    bio_evidence_type_col = str(sec.get("biological_duplicate_evidence_type_column", "biological_duplicate_evidence_type"))
    bio_evidence_value_col = str(sec.get("biological_duplicate_evidence_value_column", "biological_duplicate_evidence_value"))
    bio_class_col = str(sec.get("biological_duplicate_candidate_class_column", "biological_duplicate_candidate_class"))

    keep_roles = split_csv(sec.get("keep_roles"), "unique_input,canonical_id_row")
    exclude_roles = split_csv(sec.get("exclude_roles"), "archive_alias_row")

    rows, fields = read_rows(metadata_path)
    required = {
        accession_col, canonical_col, species_col, role_col,
        duplicate_flag_col, group_col, group_size_col, alias_col,
        bio_flag_col, bio_group_col, bio_evidence_type_col,
        bio_evidence_value_col, bio_class_col,
    }
    missing = sorted(required - set(fields))
    if missing:
        raise RuntimeError(f"metadata missing required FINAL columns {missing}: {metadata_path}")

    seen_accessions: set[str] = set()
    duplicate_accession_rows: list[str] = []
    role_counts = Counter()
    bio_class_counts = Counter()
    by_group: dict[str, list[dict]] = defaultdict(list)
    kept: list[dict] = []
    excluded: list[dict] = []
    bio_candidates: list[dict] = []
    unknown_roles: list[tuple[int, str, str]] = []

    for line_no, row in enumerate(rows, 2):
        accession = clean(row.get(accession_col))
        canonical = clean(row.get(canonical_col))
        species = clean(row.get(species_col))
        role = clean(row.get(role_col))
        group = clean(row.get(group_col))
        bio_flag = clean(row.get(bio_flag_col)) or "False"
        bio_group = clean(row.get(bio_group_col))
        bio_class = clean(row.get(bio_class_col))

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
        if bio_class:
            bio_class_counts[bio_class] += 1
        if group:
            by_group[group].append(row)

        base = {
            "sample": accession,
            "species": species,
            "canonical_sample_id": canonical,
            "archive_duplicate_flag": clean(row.get(duplicate_flag_col)) or "False",
            "archive_duplicate_group": group or "NA",
            "archive_duplicate_group_size": clean(row.get(group_size_col)) or "NA",
            "archive_duplicate_role": role,
            "archive_duplicate_alias_accession": clean(row.get(alias_col)) or "NA",
            "biological_duplicate_candidate_flag": bio_flag,
            "biological_duplicate_candidate_group": bio_group or "NA",
            "biological_duplicate_evidence_type": clean(row.get(bio_evidence_type_col)) or "NA",
            "biological_duplicate_evidence_value": clean(row.get(bio_evidence_value_col)) or "NA",
            "biological_duplicate_candidate_class": bio_class or "NA",
        }

        if bio_flag.lower() == "true":
            bio_candidates.append(base.copy())

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
        raise RuntimeError(f"metadata contains unhandled archive_duplicate_role values: {examples}")

    # Every declared archive alias group must contain exactly one canonical row
    # and at least one archive alias row.
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
        raise RuntimeError(f"invalid archive duplicate groups in FINAL metadata: {examples}")

    kept.sort(key=lambda r: r["sample"])
    excluded.sort(key=lambda r: r["sample"])
    bio_candidates.sort(key=lambda r: (r["biological_duplicate_candidate_group"], r["sample"]))

    kept_samples = {r["sample"] for r in kept}
    if len(kept_samples) != len(kept):
        raise RuntimeError("deduplicated cohort still contains duplicate sample accessions")

    detail_fields = [
        "sample", "species", "canonical_sample_id",
        "archive_duplicate_flag", "archive_duplicate_group",
        "archive_duplicate_group_size", "archive_duplicate_role",
        "archive_duplicate_alias_accession",
        "biological_duplicate_candidate_flag",
        "biological_duplicate_candidate_group",
        "biological_duplicate_evidence_type",
        "biological_duplicate_evidence_value",
        "biological_duplicate_candidate_class",
    ]

    write_tsv(reports / "deduplicated_samples.tsv", kept, detail_fields)
    write_tsv(reports / "excluded_archive_aliases.tsv", excluded, detail_fields)
    write_tsv(reports / "biological_duplicate_candidates.tsv", bio_candidates, detail_fields)
    write_tsv(
        reports / "deduplicated_sample_ref_file.tsv",
        [{"sample": r["sample"], "species": r["species"]} for r in kept],
        ["sample", "species"],
    )

    # CROSS_BIOSAMPLE_CANDIDATE is intentionally report-only. It is evidence for
    # possible biological duplication, not a confirmed archive alias.
    cross_biosample = [
        r for r in bio_candidates
        if r["biological_duplicate_candidate_class"] == "CROSS_BIOSAMPLE_CANDIDATE"
    ]

    summary = [{
        "metadata": str(metadata_path),
        "input_rows": len(rows),
        "unique_input_accessions": len(seen_accessions),
        "kept_accessions": len(kept),
        "excluded_archive_alias_accessions": len(excluded),
        "archive_duplicate_groups": len(by_group),
        "unique_input_rows": role_counts.get("unique_input", 0),
        "canonical_id_rows": role_counts.get("canonical_id_row", 0),
        "archive_alias_rows": role_counts.get("archive_alias_row", 0),
        "biological_duplicate_candidate_rows": len(bio_candidates),
        "archive_alias_only_candidate_rows": bio_class_counts.get("ARCHIVE_ALIAS_ONLY", 0),
        "cross_biosample_candidate_rows": len(cross_biosample),
    }]
    write_tsv(
        reports / "sample_deduplication_summary.tsv",
        summary,
        [
            "metadata", "input_rows", "unique_input_accessions", "kept_accessions",
            "excluded_archive_alias_accessions", "archive_duplicate_groups",
            "unique_input_rows", "canonical_id_rows", "archive_alias_rows",
            "biological_duplicate_candidate_rows", "archive_alias_only_candidate_rows",
            "cross_biosample_candidate_rows",
        ],
    )

    print(f"[sample_deduplication] metadata={metadata_path}")
    print(f"[sample_deduplication] input_rows={len(rows)}")
    print(
        f"[sample_deduplication] kept={len(kept)} "
        f"excluded_archive_aliases={len(excluded)} archive_duplicate_groups={len(by_group)}"
    )
    print(
        f"[sample_deduplication] biological_duplicate_candidates={len(bio_candidates)} "
        f"cross_biosample_candidates={len(cross_biosample)} (report-only)"
    )
    print(f"[sample_deduplication] sample_ref={reports / 'deduplicated_sample_ref_file.tsv'}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, ValueError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
