#!/usr/bin/env python3
"""Extend the native-coordinate artifact removal list to low-AF SNVs.

The historical local heteroplasmy workflow is intentionally retained for strict
HET discovery (0.10 < AF < 0.95).  This post-processing step adds only
high-specificity 0.01-0.10 SNV artifact calls so contamination scoring can use
an artifact-clean 0.01-0.20 AF universe without changing the validated strict-HET
cluster thresholds.

Low-AF SNVs are removed when at least one of the following is true:
  1. same-species/reference exact NUMT-supported native position;
  2. exact coordinate overlap with a raw PASS indel;
  3. <= configured distance from a raw PASS indel and AF-matched to that indel;
  4. member of a HIGH_CONFIDENCE indel-complex region.

The combined SOURCE-keyed removal list remains
reports/numt_variants_to_remove.tsv for downstream compatibility.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from qc_analysis.lib.simple_yaml import read_simple_yaml


def resolve(value) -> Path:
    p = Path(str(value)).expanduser()
    return p if p.is_absolute() else ROOT / p


def read_tsv(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def fnum(value):
    try:
        x = float(value)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def inum(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def yes(value) -> bool:
    return str(value or "").strip().upper() in {"YES", "TRUE", "T", "1"}


def parse_info(text: str) -> dict[str, str]:
    out = {}
    for item in str(text or "").split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
    return out


def parse_fmt(keys: str, vals: str) -> dict[str, str]:
    kk = str(keys).split(":")
    vv = str(vals).split(":")
    return {k: vv[i] for i, k in enumerate(kk) if i < len(vv)}


def float_list(value) -> list[float | None]:
    if value in {None, "", "."}:
        return []
    return [fnum(x) for x in str(value).split(",")]


def parse_contig(line: str) -> tuple[str | None, int | None]:
    if not line.startswith("##contig=<"):
        return None, None
    m_id = re.search(r"(?:^|[,<])ID=([^,>]+)", line)
    m_len = re.search(r"(?:^|[,<])(?:length|Length|LENGTH)=([0-9]+)", line)
    if not m_id or not m_len:
        return None, None
    return m_id.group(1).strip('"'), inum(m_len.group(1))


def circular_distance(a: int, b: int, length: int | None) -> int:
    d = abs(a - b)
    if length is None or length <= 0 or a < 1 or b < 1 or a > length or b > length:
        return d
    return min(d, length - d)


def find_vcf(vcf_dir: Path, patterns: list[str], sample: str) -> Path | None:
    for pat in patterns:
        p = vcf_dir / pat.format(sample=sample)
        if p.is_file():
            return p
    return None


def read_metadata(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    """Return sample->species and optional sample->reference_key."""
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows:
        return {}, {}

    header = rows[0]
    lower = [x.strip().lower() for x in header]
    headered = "sample" in lower and "species" in lower
    species, refkey = {}, {}

    if headered:
        si = lower.index("sample")
        spi = lower.index("species")
        ref_aliases = ("reference_key", "reference", "coordinate_reference", "ref_key")
        ri = next((lower.index(x) for x in ref_aliases if x in lower), None)
        for vals in rows[1:]:
            if len(vals) <= max(si, spi):
                continue
            s, sp = vals[si].strip(), vals[spi].strip()
            if not s or not sp:
                continue
            species[s] = sp
            if ri is not None and ri < len(vals) and vals[ri].strip():
                refkey[s] = vals[ri].strip()
    else:
        for vals in rows:
            if len(vals) >= 2 and vals[0].strip() and vals[1].strip():
                species[vals[0].strip()] = vals[1].strip()
    return species, refkey


def read_vcf_context(
    path: Path,
    sample: str,
    low_min: float,
    low_max: float,
    indel_min: float,
    indel_max: float,
    dp_min: float,
    pass_only: bool,
) -> tuple[list[dict], list[dict], dict[str, int]]:
    opener = gzip.open if path.name.endswith(".gz") else open
    low_snvs, indels, lengths = [], [], {}

    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##contig=<"):
                chrom, length = parse_contig(line)
                if chrom and length:
                    lengths[chrom] = length
                continue
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            chrom, pos_text, _, ref, alt_text, _, filt, info_text = f[:8]
            if pass_only and filt not in {"PASS", "."}:
                continue
            pos = inum(pos_text)
            if pos is None or not alt_text or alt_text == ".":
                continue

            info = parse_info(info_text)
            fmt = parse_fmt(f[8], f[9])
            dp = fnum(fmt.get("DP"))
            if dp is None:
                dp = fnum(info.get("DP"))
            if dp is None or dp < dp_min:
                continue

            afs = float_list(fmt.get("AF"))
            if not afs or all(x is None for x in afs):
                afs = float_list(info.get("AF"))
            ads = float_list(fmt.get("AD"))
            alts = alt_text.split(",")

            for j, alt in enumerate(alts):
                if not alt or alt == "*" or (alt.startswith("<") and alt.endswith(">")):
                    continue
                af = afs[j] if j < len(afs) else None
                if af is None and len(ads) >= j + 2:
                    finite = [x for x in ads if x is not None]
                    total = sum(finite) if finite else 0.0
                    alt_ad = ads[j + 1]
                    if total > 0 and alt_ad is not None:
                        af = alt_ad / total
                if af is None:
                    continue

                src_chrom = str(info.get("SOURCE_CHROM", chrom))
                src_pos = inum(info.get("SOURCE_POS")) or pos
                src_ref = str(info.get("SOURCE_REF", ref)).upper()
                src_alt = str(info.get("SOURCE_ALT", alt)).upper()

                is_snv = (
                    len(src_ref) == 1 and len(src_alt) == 1
                    and src_ref in "ACGT" and src_alt in "ACGT"
                )
                is_indel = len(src_ref) != len(src_alt)
                row = {
                    "sample": sample,
                    "source_chrom": src_chrom,
                    "source_pos": src_pos,
                    "source_ref": src_ref,
                    "source_alt": src_alt,
                    "source_af": af,
                    "source_dp": dp,
                }
                if is_snv and low_min <= af <= low_max:
                    low_snvs.append(row)
                elif is_indel and indel_min <= af < indel_max:
                    indels.append(row)

    return low_snvs, indels, lengths


def removal_key(row: dict) -> tuple[str, str, str, str, str]:
    return (
        str(row.get("sample", "")),
        str(row.get("source_chrom", "")),
        str(row.get("source_pos", "")),
        str(row.get("source_ref", "")),
        str(row.get("source_alt", "")),
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, required=True)
    args = ap.parse_args()

    cfg = read_simple_yaml(args.config)
    local = cfg.get("local_heteroplasmy_qc", {}) or {}
    ext = local.get("low_af_artifact_filter", {}) or {}
    if ext.get("enabled", True) is False:
        print("[low_af_artifact_filter] disabled", file=sys.stderr)
        return 0

    report_dir = resolve(local.get("output_dir", "results/qc/local_heteroplasmy_qc")) / "reports"
    combined_path = report_dir / "numt_variants_to_remove.tsv"
    existing = read_tsv(combined_path)
    if not combined_path.is_file():
        raise FileNotFoundError(
            f"strict-HET artifact removal list not found: {combined_path}; run local artifact QC first"
        )

    coord = cfg.get("coordinate_liftover", {}) or {}
    coord_paths = coord.get("paths", {}) or {}
    vcf_dir = resolve(coord_paths.get("vcf_dir", "results/qc/pre_liftover_variant_qc/vcf_source_qc"))
    raw_patterns = str(coord_paths.get("vcf_pattern", "{sample}.source_qc.vcf.gz,{sample}.source_qc.vcf"))
    patterns = [x.strip() for x in raw_patterns.split(",") if x.strip()]

    metadata_path = resolve(local.get("sample_ref_file", "config/sample_ref_file.tsv"))
    sample_species, sample_refkey = read_metadata(metadata_path)

    low_min = float(ext.get("low_af_min", 0.01))
    low_max = float(ext.get("low_af_max", 0.10))
    indel_min = float(ext.get("indel_af_min", 0.01))
    indel_max = float(ext.get("indel_af_max", 0.95))
    dp_min = float(ext.get("dp_min", local.get("dp_min", 100)))
    pass_only = ext.get("pass_only", True) is not False
    near_bp = int(ext.get("near_indel_distance_bp", 10))
    near_delta = float(ext.get("near_indel_max_delta_af", 0.05))
    use_numt = ext.get("apply_species_numt_positions", True) is not False
    use_high_conf_region = ext.get("apply_high_conf_indel_regions", True) is not False

    numt_rows = read_tsv(report_dir / "species_numt_positions.tsv")
    numt_catalog = defaultdict(set)
    species_refs = defaultdict(set)
    for row in numt_rows:
        sp = str(row.get("species", ""))
        rk = str(row.get("reference_key", ""))
        pos = inum(row.get("source_pos"))
        if sp and rk and pos is not None:
            numt_catalog[(sp, rk)].add(pos)
            species_refs[sp].add(rk)

    high_regions = defaultdict(list)
    for row in read_tsv(report_dir / "indel_complex_regions.tsv"):
        if not yes(row.get("high_confidence")):
            continue
        sample = str(row.get("sample", ""))
        triggers = [inum(x) for x in str(row.get("trigger_indel_positions", "")).split(",")]
        triggers = [x for x in triggers if x is not None]
        if not sample or not triggers:
            continue
        high_regions[sample].append({
            "region_id": str(row.get("region_id", "")),
            "triggers": triggers,
            "radius": inum(row.get("radius_bp")) or 100,
            "length": inum(row.get("native_mt_length")),
        })

    low_removals = []
    reason_counts = Counter()
    missing_vcf = []

    for sample in sorted(sample_species):
        vcf = find_vcf(vcf_dir, patterns, sample)
        if vcf is None:
            missing_vcf.append(sample)
            continue
        low_snvs, indels, lengths = read_vcf_context(
            vcf, sample, low_min, low_max, indel_min, indel_max, dp_min, pass_only
        )
        species = sample_species[sample]

        rk = sample_refkey.get(sample, "")
        if not rk and len(species_refs.get(species, set())) == 1:
            rk = next(iter(species_refs[species]))

        by_chrom_indels = defaultdict(list)
        for indel in indels:
            by_chrom_indels[str(indel["source_chrom"])].append(indel)

        for snv in low_snvs:
            reasons = []
            details = []
            pos = int(snv["source_pos"])
            chrom = str(snv["source_chrom"])
            af = float(snv["source_af"])
            length = lengths.get(chrom)

            if use_numt and rk and pos in numt_catalog.get((species, rk), set()):
                reasons.append("LOW_AF_SAME_SPECIES_NUMT_POSITION")
                details.append(f"numt_reference_key={rk}")

            local_indels = by_chrom_indels.get(chrom, [])
            nearest = None
            best_delta = None
            if local_indels:
                distance_pairs = [
                    (circular_distance(pos, int(ind["source_pos"]), length), ind)
                    for ind in local_indels
                ]
                nearest, nearest_indel = min(distance_pairs, key=lambda x: x[0])
                best_delta = min(abs(af - float(ind["source_af"])) for _, ind in distance_pairs)
                if nearest == 0:
                    reasons.append("LOW_AF_DIRECT_INDEL_OVERLAP")
                elif nearest <= near_bp:
                    matched = min(
                        (abs(af - float(ind["source_af"])) for d, ind in distance_pairs if d <= near_bp),
                        default=None,
                    )
                    if matched is not None and matched <= near_delta + 1e-12:
                        reasons.append("LOW_AF_NEAR_INDEL_AF_MATCHED")
                        best_delta = matched

            if use_high_conf_region:
                for region in high_regions.get(sample, []):
                    rlen = region["length"] or length
                    if min(circular_distance(pos, t, rlen) for t in region["triggers"]) <= region["radius"]:
                        reasons.append("LOW_AF_INDEL_COMPLEX_REGION_HIGH_CONF")
                        details.append(f"region_id={region['region_id']}")
                        break

            if not reasons:
                continue

            row = {
                "sample": sample,
                "species": species,
                "source_chrom": chrom,
                "source_pos": pos,
                "source_ref": snv["source_ref"],
                "source_alt": snv["source_alt"],
                "source_af": f"{af:.8g}",
                "source_dp": f"{float(snv['source_dp']):.8g}",
                "artifact_variant_class": "LOW_AF_SNV",
                "artifact_af_min": f"{low_min:.6g}",
                "artifact_af_max": f"{low_max:.6g}",
                "artifact_nearest_indel_distance_bp": "" if nearest is None else nearest,
                "artifact_best_indel_delta_af": "" if best_delta is None else f"{best_delta:.6g}",
                "artifact_details": ";".join(details),
                "filter_action": "REMOVE",
                "filter_reason": ";".join(sorted(set(reasons))),
            }
            low_removals.append(row)
            for reason in set(reasons):
                reason_counts[reason] += 1

    existing_by_key = {removal_key(row): dict(row) for row in existing}
    added = 0
    for row in low_removals:
        k = removal_key(row)
        if k in existing_by_key:
            old = existing_by_key[k]
            old_reasons = set(x for x in str(old.get("filter_reason", "")).split(";") if x)
            old_reasons.update(x for x in str(row.get("filter_reason", "")).split(";") if x)
            old["filter_reason"] = ";".join(sorted(old_reasons))
            if not old.get("artifact_variant_class"):
                old["artifact_variant_class"] = row["artifact_variant_class"]
            continue
        existing_by_key[k] = row
        added += 1

    combined = list(existing_by_key.values())
    fields = []
    for row in combined:
        for field in row:
            if field not in fields:
                fields.append(field)
    mandatory = [
        "sample", "species", "source_chrom", "source_pos", "source_ref", "source_alt",
        "source_af", "source_dp", "filter_action", "filter_reason",
    ]
    fields = mandatory + [x for x in fields if x not in mandatory]

    low_fields = [
        "sample", "species", "source_chrom", "source_pos", "source_ref", "source_alt",
        "source_af", "source_dp", "artifact_variant_class", "artifact_af_min", "artifact_af_max",
        "artifact_nearest_indel_distance_bp", "artifact_best_indel_delta_af", "artifact_details",
        "filter_action", "filter_reason",
    ]
    write_tsv(report_dir / "low_af_artifact_variants_to_remove.tsv", low_removals, low_fields)
    write_tsv(report_dir / "artifact_variants_to_remove.tsv", combined, fields)
    # Legacy filename retained because intra/inter/final filter already consume it.
    write_tsv(combined_path, combined, fields)

    summary = [{
        "low_af_min": low_min,
        "low_af_max": low_max,
        "dp_min": dp_min,
        "near_indel_distance_bp": near_bp,
        "near_indel_max_delta_af": near_delta,
        "n_existing_strict_artifact_removals": len(existing),
        "n_low_af_artifact_rows": len(low_removals),
        "n_new_unique_low_af_removals": added,
        "n_combined_artifact_removals": len(combined),
        "n_samples_missing_source_vcf": len(missing_vcf),
        **{f"n_{k.lower()}": v for k, v in sorted(reason_counts.items())},
    }]
    write_tsv(report_dir / "low_af_artifact_filter_summary.tsv", summary, list(summary[0]))

    if missing_vcf:
        write_tsv(
            report_dir / "low_af_artifact_missing_vcf.tsv",
            [{"sample": x} for x in missing_vcf],
            ["sample"],
        )

    print(
        "[low_af_artifact_filter] "
        f"low_af_rows={len(low_removals)} new_unique={added} "
        f"combined={len(combined)} missing_vcf={len(missing_vcf)}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
