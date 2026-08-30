#!/usr/bin/env python3
"""Low-memory implementation of terminal sample/variant filtering."""
from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import re
import shutil
import sqlite3
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from qc_analysis.lib.simple_yaml import read_simple_yaml

SAMPLE_COLUMNS = "sample species intraspecies_status human_contamination_status interspecies_status sample_level_qc_status final_sample_status final_sample_fail_reasons final_sample_warnings vcf_source".split()
VARIANT_COLUMNS = (
    "sample species human_chrom human_pos human_ref human_alt source_chrom source_pos source_ref source_alt "
    "AF DP vcf_filter variant_class call_class snv_type mt_median_coverage Percent_100 nuclear_median_coverage mtcn_median MAD "
    "sample_level_qc_status sample_failed_criteria intraspecies_status human_contamination_status interspecies_status "
    "region_type orthology_match_status orthology_fail_reason codon_match_status trna_match_status rrna_match_status "
    "final_variant_status final_variant_fail_reasons original_chrom original_pos original_ref original_alt liftover_status sample_variant_qc_status match_status"
).split()
MANAGED_OUTPUT_DIRS = ("reports", "final_vcf", "final_cov", "final_mtcn", ".work")
ANNOTATION_FIELDS = ("region_type", "orthology_match_status", "orthology_fail_reason")


def resolve(value):
    path = Path(str(value)).expanduser()
    return path if path.is_absolute() else ROOT / path


def iter_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def read_tsv(path):
    return list(iter_tsv(path))


def pick(row, names, default="NOT_AVAILABLE"):
    if isinstance(names, str):
        names = [value.strip() for value in names.split(",") if value.strip()]
    low = {key.lower(): value for key, value in row.items()}
    return next((low[name.lower()] for name in names if low.get(name.lower(), "") != ""), default)


def index(path):
    return {pick(row, ["sample", "Sample"]): row for row in iter_tsv(path)}


def names(value):
    return [item.strip() for item in value.split(",")] if isinstance(value, str) else list(value or [])


def is_fail(value, configured):
    return value.strip().lower() in {str(item).strip().lower() for item in names(configured)}


def sample_name(path):
    return path.name.split(".lifted", 1)[0].split(".vcf", 1)[0]


def open_vcf(path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open(encoding="utf-8")


def parse_info(value):
    return {
        item.split("=", 1)[0]: item.split("=", 1)[1] if "=" in item else True
        for item in value.split(";")
        if item and item != "."
    }


def number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def variant_evidence(fields):
    info = parse_info(fields[7])
    format_names = fields[8].split(":") if len(fields) > 8 else []
    format_values = fields[9].split(":") if len(fields) > 9 else []
    sample = dict(zip(format_names, format_values))
    af = number(sample.get("AF", "").split(",")[0])
    dp = number(sample.get("DP"))
    if dp is None:
        dp = number(info.get("DP"))
    if af is None and "," not in fields[4]:
        ad = [number(value) for value in sample.get("AD", "").split(",")]
        if len(ad) == 2 and None not in ad and sum(ad) > 0:
            af = ad[1] / sum(ad)
    return info, af, dp


def variant_classes(ref, alt):
    ref, alt = ref.upper(), alt.upper()
    simple = "," not in alt
    if simple and len(ref) == len(alt) == 1:
        variant_class = "SNV"
    elif simple and len(ref) != len(alt) and all(re.fullmatch(r"[ACGTN]+", allele) for allele in (ref, alt)):
        variant_class = "INDEL"
    else:
        variant_class = "OTHER"
    if variant_class == "SNV" and {ref, alt} in ({"A", "G"}, {"C", "T"}):
        snv_type = "SNV_transition"
    elif variant_class == "SNV" and ref in "ACGT" and alt in "ACGT":
        snv_type = "SNV_transversion"
    elif variant_class == "INDEL":
        snv_type = "indel"
    else:
        snv_type = "other"
    return variant_class, snv_type


def call_class(af):
    if af is None:
        return "UNKNOWN"
    if af >= 0.95:
        return "homoplasmic"
    if af >= 0.10:
        return "heteroplasmic"
    return "low_af"


def info_value(info, *aliases):
    return next(
        (str(info[name]) for name in aliases if info.get(name) not in (None, "", ".")),
        "NOT_AVAILABLE",
    )


LEGACY_SUFFIXES = (
    ".lifted.codon.trna.rrna.vcf",
    ".lifted.codon.trna.vcf",
    ".lifted.trna.vcf",
    ".lifted.codon.vcf",
    ".lifted.raw.vcf",
    ".vcf",
)


def source_specs(value):
    if isinstance(value, dict):
        return [
            (name, resolve(spec.get("dir", spec.get("directory"))), spec.get("pattern"))
            for name, spec in value.items()
            if isinstance(spec, dict)
        ]
    return [(str(directory), resolve(directory), None) for directory in names(value)]


def find_vcf(directory, sample, pattern=None):
    base = [directory / pattern.format(sample=sample)] if pattern else [directory / f"{sample}{suffix}" for suffix in LEGACY_SUFFIXES]
    candidates = []
    for path in base:
        candidates.extend(candidate for candidate in (path, Path(str(path) + ".gz")) if candidate.is_file())
    candidates = sorted(set(candidates))
    if len(candidates) > 1:
        raise ValueError(
            f"ambiguous VCF source for sample {sample} in {directory}: {', '.join(map(str, candidates))}"
        )
    return candidates[0] if candidates else None


def liftover_status(sample, cfg):
    section = cfg.get("coordinate_liftover") or {}
    paths = section.get("paths") or {}
    report = resolve(paths.get("output_dir", "results/qc/coordinate_liftover")) / "reports" / f"{sample}.coordinate_liftover_qc.tsv"
    if not report.is_file():
        return "NOT_AVAILABLE"
    with report.open(newline="", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) >= 2 and row[0].strip() == "status":
                return row[1].strip() or "NOT_AVAILABLE"
    return "NOT_AVAILABLE"


def report_variant_key(row, spec, source):
    sample = pick(row, ["sample", "Sample"], "")
    human = [pick(row, [f"human_{field}"], "") for field in ("chrom", "pos", "ref", "alt")]
    if all(human):
        return (sample, *human)
    system = str(spec.get("coordinate_system", "")).strip().lower()
    if system not in {"human", "post-liftover", "post_liftover"}:
        raise ValueError(
            f"variant report {source} uses generic coordinates but coordinate_system is unknown or incompatible: {system or 'not declared'}"
        )
    generic = [pick(row, [field.upper(), field], "") for field in ("chrom", "pos", "ref", "alt")]
    if not all(generic):
        raise ValueError(f"variant report {source} lacks a complete human-coordinate variant key")
    return (sample, *generic)


def reset_managed_outputs(out):
    """Remove derived products while preserving scheduler manifests and logs."""
    out.mkdir(parents=True, exist_ok=True)
    for name in MANAGED_OUTPUT_DIRS:
        path = out / name
        if path.is_symlink() or path.is_file():
            path.unlink()
        elif path.exists():
            shutil.rmtree(path)
    for name in MANAGED_OUTPUT_DIRS:
        (out / name).mkdir(parents=True, exist_ok=True)
    (out / "logs").mkdir(parents=True, exist_ok=True)


def build_variant_report_store(db_path, variant_reports, passing_samples):
    """Stream variant reports into a disk-backed store for low-memory lookup."""
    connection = sqlite3.connect(db_path)
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute("PRAGMA temp_store=FILE")
    connection.execute(
        """
        CREATE TABLE evidence (
            sample TEXT NOT NULL,
            chrom TEXT NOT NULL,
            pos TEXT NOT NULL,
            ref TEXT NOT NULL,
            alt TEXT NOT NULL,
            source_order INTEGER NOT NULL,
            row_order INTEGER NOT NULL,
            fail_reason TEXT NOT NULL,
            region_type TEXT NOT NULL,
            orthology_match_status TEXT NOT NULL,
            orthology_fail_reason TEXT NOT NULL
        )
        """
    )
    for source_order, (source, spec) in enumerate(variant_reports.items()):
        path = resolve(spec["path"] if isinstance(spec, dict) else spec)
        if not path.is_file():
            continue
        if not isinstance(spec, dict):
            raise ValueError(f"variant report {source} must declare path and coordinate_system")
        status_names = names(spec.get("status_columns", [])) + ["qc_status", "status", "match_status"]
        fail_status = spec.get("fail_status", ["FAIL"])
        batch = []
        for row_order, row in enumerate(iter_tsv(path)):
            key = report_variant_key(row, spec, source)
            if key[0] not in passing_samples:
                continue
            status = pick(row, status_names, "PASS")
            fail_reason = source + ":" + status if is_fail(status, fail_status) else ""
            annotation = {field: pick(row, [field], "") for field in ANNOTATION_FIELDS}
            batch.append(
                (
                    *key,
                    source_order,
                    row_order,
                    fail_reason,
                    annotation["region_type"],
                    annotation["orthology_match_status"],
                    annotation["orthology_fail_reason"],
                )
            )
            if len(batch) >= 10000:
                connection.executemany("INSERT INTO evidence VALUES (?,?,?,?,?,?,?,?,?,?,?)", batch)
                batch.clear()
        if batch:
            connection.executemany("INSERT INTO evidence VALUES (?,?,?,?,?,?,?,?,?,?,?)", batch)
        connection.commit()
    connection.execute("CREATE INDEX evidence_sample_idx ON evidence(sample)")
    connection.commit()
    return connection


def load_sample_variant_evidence(connection, sample):
    """Load only one sample's report-derived evidence into memory."""
    flags = defaultdict(list)
    annotations = defaultdict(dict)
    rows = connection.execute(
        """
        SELECT chrom,pos,ref,alt,fail_reason,region_type,orthology_match_status,orthology_fail_reason
        FROM evidence
        WHERE sample=?
        ORDER BY source_order,row_order
        """,
        (sample,),
    )
    for chrom, pos, ref, alt, fail_reason, region_type, orthology_status, orthology_reason in rows:
        key = (chrom, pos, ref, alt)
        if region_type:
            annotations[key]["region_type"] = region_type
        if orthology_status:
            annotations[key]["orthology_match_status"] = orthology_status
        if orthology_reason:
            annotations[key]["orthology_fail_reason"] = orthology_reason
        if fail_reason:
            flags[key].append(fail_reason)
    return flags, annotations


def bgzip_and_index(plain, dest):
    try:
        import pysam

        pysam.tabix_compress(str(plain), str(dest), force=True)
        pysam.tabix_index(str(dest), preset="vcf", force=True)
        return
    except ImportError:
        pass
    bgzip, tabix = shutil.which("bgzip"), shutil.which("tabix")
    if not bgzip or not tabix:
        raise RuntimeError("BGZF/index output requires pysam or both bgzip and tabix")
    with dest.open("wb") as handle:
        subprocess.run([bgzip, "-c", str(plain)], stdout=handle, check=True)
    subprocess.run([tabix, "-f", "-p", "vcf", str(dest)], check=True)


def sort_plain_vcf(input_vcf: Path, output_vcf: Path) -> None:
    contig_order = {}
    headers = []
    records = []
    with input_vcf.open("r", encoding="utf-8", newline="") as source:
        for line_number, line in enumerate(source, 1):
            if line.startswith("#"):
                headers.append(line)
                match = re.match(r"^##contig=<ID=([^,>]+)", line)
                if match:
                    contig = match.group(1).strip().strip('"')
                    if contig not in contig_order:
                        contig_order[contig] = len(contig_order)
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 5:
                raise ValueError(f"invalid VCF data line {line_number} in {input_vcf}")
            try:
                pos = int(fields[1])
            except ValueError as exc:
                raise ValueError(
                    f"non-integer VCF POS on line {line_number} in {input_vcf}: {fields[1]!r}"
                ) from exc
            chrom, ref, alt = fields[0], fields[3], fields[4]
            contig_key = (0, contig_order[chrom]) if chrom in contig_order else (1, chrom)
            records.append(((contig_key, pos, ref, alt, line), line))
    records.sort(key=lambda item: item[0])
    with output_vcf.open("w", encoding="utf-8", newline="") as target:
        target.writelines(headers)
        target.writelines(line for _key, line in records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    cfg = read_simple_yaml(args.config)
    sec = cfg.get("final_filter") or {}
    if sec.get("enabled", True) is False:
        print("[final_filter] disabled; skipping.")
        return 0

    out = resolve(sec.get("output_dir", "results/qc/final_filter"))
    collected = resolve(sec.get("collected_dir", "results/qc/collected_variant_calling_results"))
    collection_path = collected / "reports/variant_calling_collection_summary.tsv"
    if not collection_path.is_file():
        raise ValueError(f"missing collection summary: {collection_path}")

    sources = sec.get("sample_reports") or {}
    defaults = {
        "intraspecies": ("results/qc/intraspecies_contamination/reports/intraspecies_contamination_report.tsv", ["contamination_status", "qc_status"]),
        "human": ("results/qc/human_contamination/reports/human_contamination_report.tsv", ["human_contamination_status", "qc_status", "status"]),
        "interspecies": ("results/qc/interspecies_contamination/reports/interspecies_contamination_report.tsv", ["interspecies_status", "qc_status", "status"]),
        "sample_qc": ("results/qc/sample_variant_filtering/reports/sample_qc.tsv", ["qc_status"]),
    }
    required = set(names(sec.get("required_sample_reports", ["intraspecies", "sample_qc"])))
    optional = set(names(sec.get("optional_sample_reports", ["human", "interspecies"])))
    active = required | optional
    unknown = active - set(defaults)
    if unknown:
        raise ValueError(f"unknown sample reports: {sorted(unknown)}")

    indexed = {}
    fields = {}
    missing = []
    for name, (fallback, candidates) in defaults.items():
        spec = sources.get(name, {}) if isinstance(sources, dict) else {}
        fields[name] = spec.get("status_columns", candidates)
        if name not in active:
            indexed[name] = {}
            continue
        path = resolve(spec.get("path", fallback))
        if name in required and not path.is_file():
            missing.append(f"{name}={path}")
        indexed[name] = index(path) if path.is_file() else {}
    if missing:
        raise ValueError("missing required sample report(s): " + ", ".join(missing))

    collection = index(collection_path)
    strict = sec.get("strict_missing_samples", True) is not False
    for report in required:
        absent = sorted(set(collection) - set(indexed[report]))
        if absent and strict:
            raise ValueError(f"required report {report} is missing samples: {', '.join(absent)}")

    reset_managed_outputs(out)

    fail_cfg = sec.get("sample_fail_status") or {}
    fail_defaults = {
        "intraspecies": ["high_confidence_contaminated"],
        "human": ["FAIL"],
        "interspecies": ["FAIL"],
        "sample_qc": ["FAIL"],
    }
    vcf_sources = source_specs(
        sec.get(
            "vcf_sources",
            [
                "results/qc/rrna_match/vcf_rrna",
                "results/qc/trna_match/vcf_trna",
                "results/qc/codon_match/vcf_codon",
                "results/qc/coordinate_liftover/vcf_lifted_raw",
            ],
        )
    )

    sample_rows = []
    sample_context = {}
    passing = {}
    sample_qc_rows = indexed["sample_qc"]
    for sample, row in sorted(collection.items()):
        statuses = {name: pick(indexed[name].get(sample, {}), fields[name]) for name in defaults}
        reasons = []
        warnings = []
        for name, value in statuses.items():
            if name not in active:
                continue
            if is_fail(value, fail_cfg.get(name, fail_defaults[name])):
                if name == "sample_qc":
                    failed = pick(sample_qc_rows.get(sample, {}), ["failed_criteria"], "failed")
                    reasons.extend("sample_qc:" + criterion for criterion in failed.split(";") if criterion)
                else:
                    reasons.append(name + ":" + value)
            elif name == "intraspecies" and (value.startswith("insufficient_") or value == "candidate_contaminated"):
                warnings.append(name + ":" + value)
            elif name == "human" and value.upper() in {"CANDIDATE", "INSUFFICIENT_DATA"}:
                warnings.append(name + ":" + value)
            elif value == "NOT_AVAILABLE" and name in required:
                reasons.append(name + ":missing")
        lift = liftover_status(sample, cfg)
        if lift != "completed":
            reasons.append("liftover:" + lift)
        src = (
            next((found for _, directory, pattern in vcf_sources if (found := find_vcf(directory, sample, pattern))), None)
            if lift == "completed"
            else None
        )
        if not reasons and src is None:
            reasons.append("vcf:missing_downstream_source")
        status = "FAIL" if reasons else "PASS"
        if status == "PASS":
            passing[sample] = src
        sample_row = dict(
            sample=sample,
            species=pick(row, ["species", "Species"], ""),
            intraspecies_status=statuses["intraspecies"],
            human_contamination_status=statuses["human"],
            interspecies_status=statuses["interspecies"],
            sample_level_qc_status=statuses["sample_qc"],
            final_sample_status=status,
            final_sample_fail_reasons=";".join(reasons),
            final_sample_warnings=";".join(warnings),
            vcf_source=str(src or ""),
        )
        sample_rows.append(sample_row)
        sample_context[sample] = sample_row

    db_path = out / ".work" / "variant_reports.sqlite3"
    report_store = build_variant_report_store(db_path, sec.get("variant_reports") or {}, set(passing))
    kept = {}
    n_variants = 0
    n_variants_fail = 0
    final_variant_report = out / "reports" / "final_variant_qc.tsv"
    variant_report_tmp = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            newline="",
            encoding="utf-8",
            prefix="final_variant_qc.tsv.tmp.",
            delete=False,
            dir=final_variant_report.parent,
        ) as report_handle:
            variant_report_tmp = Path(report_handle.name)
            variant_writer = csv.DictWriter(report_handle, fieldnames=VARIANT_COLUMNS, delimiter="\t")
            variant_writer.writeheader()
            for sample, src in sorted(passing.items()):
                variant_flags, variant_annotations = load_sample_variant_evidence(report_store, sample)
                with tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False, dir=out) as target:
                    plain = Path(target.name)
                    n = 0
                    with open_vcf(src) as inp:
                        for line in inp:
                            if line.startswith("#"):
                                target.write(line)
                                continue
                            f = line.rstrip("\n").split("\t")
                            key = (f[0], f[1], f[3], f[4])
                            why = list(variant_flags.get(key, []))
                            vcf_filter = f[6]
                            if vcf_filter != "PASS":
                                why.append("vcf_filter:" + vcf_filter)
                            info, af, dp = variant_evidence(f)
                            variant_class, snv_type = variant_classes(f[3], f[4])
                            ref, alt = f[3].upper(), f[4].upper()
                            is_canonical_snv = len(ref) == len(alt) == 1 and ref in "ACGT" and alt in "ACGT" and "," not in alt
                            if not is_canonical_snv:
                                why.append("variant_class:" + variant_class)
                            final_status = "FAIL" if why else "PASS"
                            context = sample_context[sample]
                            qc = sample_qc_rows.get(sample, {})
                            annotation = variant_annotations.get(key, {})
                            source_chrom = info_value(info, "SRC_CHROM", "MTLIFT_ORIG_CHROM")
                            source_pos = info_value(info, "SRC_POS", "MTLIFT_ORIG_POS")
                            source_ref = info_value(info, "SRC_REF", "MTLIFT_ORIG_REF")
                            source_alt = info_value(info, "SRC_ALT", "MTLIFT_ORIG_ALT")
                            orthology_status = annotation.get("orthology_match_status", "NOT_AVAILABLE")
                            variant_writer.writerow(
                                dict(
                                    sample=sample,
                                    species=context["species"],
                                    human_chrom=f[0],
                                    human_pos=f[1],
                                    human_ref=f[3],
                                    human_alt=f[4],
                                    source_chrom=source_chrom,
                                    source_pos=source_pos,
                                    source_ref=source_ref,
                                    source_alt=source_alt,
                                    AF=af if af is not None else "NA",
                                    DP=dp if dp is not None else "NA",
                                    vcf_filter=vcf_filter,
                                    variant_class=variant_class,
                                    call_class=call_class(af),
                                    snv_type=snv_type,
                                    mt_median_coverage=pick(qc, ["mt_median_coverage"]),
                                    Percent_100=pick(qc, ["Percent_100"]),
                                    nuclear_median_coverage=pick(qc, ["nuclear_median_coverage"]),
                                    mtcn_median=pick(qc, ["mtcn_median"]),
                                    MAD=pick(qc, ["MAD"]),
                                    sample_level_qc_status=context["sample_level_qc_status"],
                                    sample_failed_criteria=pick(qc, ["failed_criteria"], ""),
                                    intraspecies_status=context["intraspecies_status"],
                                    human_contamination_status=context["human_contamination_status"],
                                    interspecies_status=context["interspecies_status"],
                                    region_type=annotation.get("region_type", "NOT_AVAILABLE"),
                                    orthology_match_status=orthology_status,
                                    orthology_fail_reason=annotation.get("orthology_fail_reason", "NOT_AVAILABLE"),
                                    codon_match_status=info_value(info, "MTCODON_STATUS"),
                                    trna_match_status=info_value(info, "MTTRNA_STATUS"),
                                    rrna_match_status=info_value(info, "MTRRNA_STATUS"),
                                    final_variant_status=final_status,
                                    final_variant_fail_reasons=";".join(why),
                                    original_chrom=source_chrom,
                                    original_pos=source_pos,
                                    original_ref=source_ref,
                                    original_alt=source_alt,
                                    liftover_status="PASS",
                                    sample_variant_qc_status="PASS" if not why else "FAIL",
                                    match_status=orthology_status,
                                )
                            )
                            n_variants += 1
                            if final_status == "FAIL":
                                n_variants_fail += 1
                            else:
                                target.write(line)
                                n += 1
                dest = out / "final_vcf" / f"{sample}.final.vcf.gz"
                sorted_plain = None
                try:
                    with tempfile.NamedTemporaryFile("w", suffix=".sorted.vcf", delete=False, dir=out) as sorted_target:
                        sorted_plain = Path(sorted_target.name)
                    sort_plain_vcf(plain, sorted_plain)
                    bgzip_and_index(sorted_plain, dest)
                finally:
                    plain.unlink(missing_ok=True)
                    if sorted_plain is not None:
                        sorted_plain.unlink(missing_ok=True)
                kept[sample] = n
                for kind in ("cov", "mtcn"):
                    for source in sorted((collected / f"collected_{kind}").glob(sample + ".*")):
                        shutil.copy2(source, out / f"final_{kind}" / source.name)
        os.replace(variant_report_tmp, final_variant_report)
        variant_report_tmp = None
    finally:
        report_store.close()
        if variant_report_tmp is not None:
            variant_report_tmp.unlink(missing_ok=True)

    with (out / "reports" / "final_sample_qc.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SAMPLE_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(sample_rows)
    with (out / "reports" / "final_filter_summary.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("metric", "value"))
        writer.writerows(
            (
                ("n_samples", len(collection)),
                ("n_samples_pass", len(passing)),
                ("n_samples_fail", len(collection) - len(passing)),
                ("n_variants_pass", sum(kept.values())),
                ("n_variants_fail", n_variants_fail),
            )
        )
    shutil.rmtree(out / ".work", ignore_errors=True)
    (out / "logs" / "final_filter.log").write_text(
        f"samples={len(collection)} pass={len(passing)} variants={n_variants}\n",
        encoding="utf-8",
    )
    print(f"[final_filter] output={out} samples_pass={len(passing)}/{len(collection)}")
    return 0
