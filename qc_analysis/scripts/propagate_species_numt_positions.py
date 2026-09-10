#!/usr/bin/env python3
"""Propagate NUMT-supported HET positions to same-species samples.

This is a post-processing step for ``run_local_heteroplasmy_qc.py``.

Rationale
---------
A NUMT-associated local cluster can be obvious in one sample but partially masked
in another sample from the same species because the mitochondrial haplotype differs.
Therefore, once a strict-HET position is observed inside a NUMT-supported cluster,
that native mtDNA position becomes species-level NUMT evidence.  A strict HET at the
same native position in another sample from the same species/reference coordinate
system is marked NUMT-associated and removed even if that sample does not itself
form a local cluster.

Propagation deliberately uses species + reference_key + source_pos, not REF/ALT.
It does not propagate across different reference coordinate systems and it does not
expand to every position inside a NUMT interval; only positions actually observed
inside a NUMT-supported cluster are propagated.
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from qc_analysis.lib.simple_yaml import read_simple_yaml


def resolve(value: str | Path) -> Path:
    p = Path(str(value)).expanduser()
    return p if p.is_absolute() else ROOT / p


def read_tsv(path: Path) -> list[dict]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def cluster_reference_map(cluster_rows: list[dict]) -> dict[str, str]:
    return {
        str(row.get("cluster_id", "")): str(row.get("reference_key", ""))
        for row in cluster_rows
        if str(row.get("cluster_id", ""))
    }


def build_species_numt_position_catalogue(
    cluster_rows: list[dict],
    variant_rows: list[dict],
) -> dict[tuple[str, str, int], dict[str, set[str]]]:
    """Build exact-position catalogue from NUMT-supported cluster members only."""
    numt_clusters = {
        str(row.get("cluster_id", ""))
        for row in cluster_rows
        if str(row.get("any_numt_evidence", "")).upper() == "YES"
    }
    ref_by_cluster = cluster_reference_map(cluster_rows)
    tier_by_cluster = {
        str(row.get("cluster_id", "")): str(row.get("numt_tier", "NONE"))
        for row in cluster_rows
        if str(row.get("cluster_id", ""))
    }

    catalogue: dict[tuple[str, str, int], dict[str, set[str]]] = defaultdict(
        lambda: {
            "source_samples": set(),
            "source_cluster_ids": set(),
            "numt_tiers": set(),
        }
    )
    for row in variant_rows:
        cid = str(row.get("cluster_id", ""))
        if cid not in numt_clusters:
            continue
        species = str(row.get("species", ""))
        refkey = ref_by_cluster.get(cid, "")
        try:
            pos = int(float(row.get("source_pos", "")))
        except (TypeError, ValueError):
            continue
        if not species or not refkey:
            continue
        key = (species, refkey, pos)
        catalogue[key]["source_samples"].add(str(row.get("sample", "")))
        catalogue[key]["source_cluster_ids"].add(cid)
        tier = tier_by_cluster.get(cid, "NONE")
        if tier:
            catalogue[key]["numt_tiers"].add(tier)
    return dict(catalogue)


def catalogue_rows(catalogue: dict[tuple[str, str, int], dict[str, set[str]]]) -> list[dict]:
    out = []
    for (species, refkey, pos), info in sorted(catalogue.items()):
        source_samples = sorted(x for x in info["source_samples"] if x)
        source_clusters = sorted(x for x in info["source_cluster_ids"] if x)
        tiers = sorted(x for x in info["numt_tiers"] if x)
        out.append({
            "species": species,
            "reference_key": refkey,
            "source_pos": pos,
            "n_source_samples": len(source_samples),
            "source_samples": ",".join(source_samples),
            "n_source_clusters": len(source_clusters),
            "source_cluster_ids": ",".join(source_clusters),
            "numt_tiers": ",".join(tiers),
        })
    return out


def apply_species_numt_position_propagation(
    cluster_rows: list[dict],
    variant_rows: list[dict],
    catalogue: dict[tuple[str, str, int], dict[str, set[str]]],
    remove_source_numt_clusters: bool = True,
) -> list[dict]:
    """Annotate and remove same-position strict HETs in other same-species samples."""
    numt_clusters = {
        str(row.get("cluster_id", ""))
        for row in cluster_rows
        if str(row.get("any_numt_evidence", "")).upper() == "YES"
    }
    ref_by_cluster = cluster_reference_map(cluster_rows)
    sample_ref: dict[str, str] = {}
    for row in cluster_rows:
        sample = str(row.get("sample", ""))
        refkey = str(row.get("reference_key", ""))
        if sample and refkey:
            sample_ref[sample] = refkey

    # Samples with no detected cluster do not appear in cluster_rows.  Infer their
    # reference key only when the species has exactly one reference key in this run.
    species_refs: dict[str, set[str]] = defaultdict(set)
    for row in cluster_rows:
        species = str(row.get("species", ""))
        refkey = str(row.get("reference_key", ""))
        if species and refkey:
            species_refs[species].add(refkey)

    out = []
    for raw in variant_rows:
        row = dict(raw)
        sample = str(row.get("sample", ""))
        species = str(row.get("species", ""))
        cid = str(row.get("cluster_id", ""))
        own_numt_cluster = cid in numt_clusters

        refkey = ref_by_cluster.get(cid, "") or sample_ref.get(sample, "")
        if not refkey and len(species_refs.get(species, set())) == 1:
            refkey = next(iter(species_refs[species]))

        try:
            pos = int(float(row.get("source_pos", "")))
        except (TypeError, ValueError):
            pos = None

        info = catalogue.get((species, refkey, pos)) if pos is not None and refkey else None
        row["species_numt_position_match"] = "YES" if info else "NO"
        row["numt_propagated"] = "NO"
        row["numt_source_samples"] = ""
        row["numt_source_cluster_ids"] = ""

        if info:
            source_samples = sorted(x for x in info["source_samples"] if x)
            source_clusters = sorted(x for x in info["source_cluster_ids"] if x)
            row["numt_source_samples"] = ",".join(source_samples)
            row["numt_source_cluster_ids"] = ",".join(source_clusters)

            if own_numt_cluster and remove_source_numt_clusters:
                row["filter_action"] = "REMOVE"
                if str(row.get("filter_reason", "")) not in {
                    "SAMPLE_NUMT_OVERLAP",
                    "SPECIES_NUMT_AND_RECURRENCE",
                }:
                    row["filter_reason"] = "NUMT_SUPPORTED_CLUSTER"

            # Propagation is intentionally restricted to another sample.  A sample
            # does not create additional evidence for itself merely by sharing the
            # same position with its own cluster.
            has_other_source = any(src and src != sample for src in source_samples)
            if (not own_numt_cluster) and has_other_source:
                row["numt_propagated"] = "YES"
                row["filter_action"] = "REMOVE"
                row["filter_reason"] = "SAME_SPECIES_NUMT_POSITION"
                row["numt_scope"] = "SPECIES"
                if not str(row.get("cluster_class", "")):
                    row["cluster_class"] = "NUMT_PROPAGATED"
                tiers = sorted(x for x in info["numt_tiers"] if x and x != "NONE")
                if tiers:
                    row["numt_tier"] = ",".join(tiers)

        out.append(row)
    return out


def rebuild_removal_rows(variant_rows: list[dict]) -> list[dict]:
    rows = []
    seen = set()
    for row in variant_rows:
        if str(row.get("filter_action", "")).upper() != "REMOVE":
            continue
        key = (
            row.get("sample", ""),
            row.get("source_chrom", ""),
            row.get("source_pos", ""),
            row.get("source_ref", ""),
            row.get("source_alt", ""),
        )
        if key in seen:
            continue
        seen.add(key)
        rows.append(dict(row))
    return rows


def update_sample_rows(sample_rows: list[dict], variant_rows: list[dict]) -> list[dict]:
    removed = Counter()
    propagated = Counter()
    for row in variant_rows:
        sample = str(row.get("sample", ""))
        if str(row.get("filter_action", "")).upper() == "REMOVE":
            removed[sample] += 1
        if str(row.get("numt_propagated", "")).upper() == "YES":
            propagated[sample] += 1

    out = []
    for raw in sample_rows:
        row = dict(raw)
        sample = str(row.get("sample", ""))
        try:
            n_het = int(float(row.get("n_het", 0) or 0))
        except (TypeError, ValueError):
            n_het = 0
        n_remove = removed[sample]
        n_prop = propagated[sample]
        row["n_numt_variants_to_remove"] = n_remove
        row["n_propagated_numt_variants"] = n_prop
        row["fraction_het_removed"] = f"{n_remove / n_het:.6g}" if n_het else "0"
        try:
            own_numt = int(float(row.get("n_numt_clusters", 0) or 0)) > 0
        except (TypeError, ValueError):
            own_numt = False
        row["numt_sample"] = "YES" if own_numt or n_prop > 0 else "NO"
        out.append(row)
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()

    cfg = read_simple_yaml(args.config)
    sec = cfg.get("local_heteroplasmy_qc", {})
    report_dir = resolve(sec.get("output_dir", "results/qc/local_heteroplasmy_qc")) / "reports"

    cluster_path = report_dir / "local_heteroplasmy_cluster_summary.tsv"
    variant_path = report_dir / "local_heteroplasmy_variant_detail.tsv"
    sample_path = report_dir / "local_heteroplasmy_sample_summary.tsv"

    cluster_rows = read_tsv(cluster_path)
    variant_rows = read_tsv(variant_path)
    sample_rows = read_tsv(sample_path)

    prop_cfg = sec.get("species_numt_position_propagation", {}) or {}
    enabled = prop_cfg.get("enabled", True)
    if enabled is False:
        print("[species_numt_position_propagation] disabled", file=sys.stderr)
        return 0
    remove_source = prop_cfg.get("remove_source_numt_clusters", True) is not False

    catalogue = build_species_numt_position_catalogue(cluster_rows, variant_rows)
    annotated = apply_species_numt_position_propagation(
        cluster_rows,
        variant_rows,
        catalogue,
        remove_source_numt_clusters=remove_source,
    )
    removal_rows = rebuild_removal_rows(annotated)
    updated_samples = update_sample_rows(sample_rows, annotated)
    numt_samples = [row for row in updated_samples if str(row.get("numt_sample", "")).upper() == "YES"]

    variant_fields = list(variant_rows[0].keys()) if variant_rows else []
    for field in (
        "species_numt_position_match",
        "numt_propagated",
        "numt_source_samples",
        "numt_source_cluster_ids",
    ):
        if field not in variant_fields:
            variant_fields.append(field)

    removal_fields = [
        "sample", "species", "source_chrom", "source_pos", "source_ref", "source_alt",
        "source_af", "source_dp", "clustered", "cluster_id", "cluster_class", "numt_scope",
        "numt_tier", "recurrence", "species_numt_position_match", "numt_propagated",
        "numt_source_samples", "numt_source_cluster_ids", "filter_action", "filter_reason",
    ]
    sample_fields = list(sample_rows[0].keys()) if sample_rows else []
    if "n_propagated_numt_variants" not in sample_fields:
        insert_at = sample_fields.index("n_numt_variants_to_remove") + 1 if "n_numt_variants_to_remove" in sample_fields else len(sample_fields)
        sample_fields.insert(insert_at, "n_propagated_numt_variants")

    write_tsv(
        report_dir / "species_numt_positions.tsv",
        catalogue_rows(catalogue),
        [
            "species", "reference_key", "source_pos", "n_source_samples", "source_samples",
            "n_source_clusters", "source_cluster_ids", "numt_tiers",
        ],
    )
    write_tsv(variant_path, annotated, variant_fields)
    write_tsv(report_dir / "numt_variants_to_remove.tsv", removal_rows, removal_fields)
    write_tsv(sample_path, updated_samples, sample_fields)
    write_tsv(report_dir / "numt_samples.tsv", numt_samples, sample_fields)

    n_prop = sum(str(row.get("numt_propagated", "")).upper() == "YES" for row in annotated)
    n_prop_samples = len({row.get("sample", "") for row in annotated if str(row.get("numt_propagated", "")).upper() == "YES"})
    print(
        "[species_numt_position_propagation] "
        f"catalogue_positions={len(catalogue)} propagated_variants={n_prop} "
        f"propagated_samples={n_prop_samples} total_remove={len(removal_rows)}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
