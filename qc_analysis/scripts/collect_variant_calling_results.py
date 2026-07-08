#!/usr/bin/env python3
"""Collect variant calling outputs for preprocessing QC before coordinate liftover.

This script searches one directory per sample under an input root, collects mtCN and
VCF outputs, merges two per-base coverage files by taking the maximum depth per
(chrom, pos, target), and writes a summary table for downstream liftover inputs.
"""

from __future__ import annotations

import argparse
import configparser
import csv
import gzip
import logging
import os
import shutil
import sys
from pathlib import Path
from statistics import median
from typing import Dict, List, Optional, Sequence, Tuple

SUMMARY_COLUMNS = [
    "sample",
    "species",
    "mt_median_coverage",
    "nuclear_median_coverage",
    "mtcn_median",
    "Percent_100",
    "MAD",
    "n_hetero",
    "n_homo",
    "vcf_file",
    "cov_file",
    "mtcn_file",
    "status",
    "missing_files",
    "notes",
]

MTCN_ALIASES = {
    "mt_median_coverage": [
        "mt_median_coverage",
        "mitochondrial_median_coverage",
        "mito_median_coverage",
        "mt_median_depth",
        "mitochondrial_median_depth",
        "median_mt_coverage",
        "median_mito_coverage",
    ],
    "nuclear_median_coverage": [
        "nuclear_median_coverage",
        "nuc_median_coverage",
        "autosomal_median_coverage",
        "nuclear_median_depth",
        "nuc_median_depth",
        "median_nuclear_coverage",
    ],
    "mtcn_median": [
        "mtcn_median",
        "mtcn",
        "mt_cn",
        "mtcopy_number",
        "mt_copy_number",
        "mitochondrial_copy_number",
        "median_mtcn",
    ],
}


def norm_name(name: str) -> str:
    return "".join(ch.lower() for ch in name.strip() if ch.isalnum())


def parse_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def parse_scalar(value: str) -> object:
    value = value.strip()
    if value == "" or value.lower() in {"null", "none", "~"}:
        return None
    if (value.startswith('"') and value.endswith('"')) or (value.startswith("'") and value.endswith("'")):
        return value[1:-1]
    lowered = value.lower()
    if lowered == "true":
        return True
    if lowered == "false":
        return False
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value


def strip_yaml_comment(line: str) -> str:
    in_single = False
    in_double = False
    for idx, char in enumerate(line):
        if char == "'" and not in_double:
            in_single = not in_single
        elif char == '"' and not in_single:
            in_double = not in_double
        elif char == "#" and not in_single and not in_double:
            return line[:idx]
    return line


def read_simple_yaml(path: Path) -> Dict[str, object]:
    root: Dict[str, object] = {}
    stack: List[Tuple[int, Dict[str, object]]] = [(-1, root)]
    with path.open() as handle:
        for raw_line in handle:
            line = strip_yaml_comment(raw_line.rstrip("\n"))
            if not line.strip():
                continue
            indent = len(line) - len(line.lstrip(" "))
            text = line.strip()
            if ":" not in text:
                raise SystemExit(f"Unsupported YAML line in {path}: {raw_line.rstrip()}")
            key, value = text.split(":", 1)
            key = key.strip().strip('"\'')
            value = value.strip()
            while stack and indent <= stack[-1][0]:
                stack.pop()
            parent = stack[-1][1]
            if value == "":
                child: Dict[str, object] = {}
                parent[key] = child
                stack.append((indent, child))
            elif value == "{}":
                parent[key] = {}
            else:
                parent[key] = parse_scalar(value)
    return root


def flatten_mapping(mapping: Dict[str, object]) -> Dict[str, object]:
    return {key: value for key, value in mapping.items() if not isinstance(value, dict)}


def read_config(path: Path) -> Dict[str, object]:
    if path.suffix.lower() in {".yaml", ".yml"}:
        data = read_simple_yaml(path)
        section = data.get("collect_variant_calling", data)
        if not isinstance(section, dict):
            raise SystemExit("collect_variant_calling config must be a mapping")
        return flatten_mapping(section)

    parser = configparser.ConfigParser()
    parser.read(path)
    values: Dict[str, object] = {}
    for section in parser.sections():
        values.update(dict(parser.items(section)))
    return values


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--config", type=Path)
    known, _ = pre.parse_known_args(argv)
    defaults = {}
    if known.config:
        defaults = read_config(known.config)

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, help="YAML or INI configuration file")
    ap.add_argument("--input-root", type=Path, default=defaults.get("input_root"))
    ap.add_argument("--outdir", type=Path, default=defaults.get("outdir"))
    ap.add_argument("--metadata", type=Path, default=defaults.get("metadata"), help="Optional sample/species metadata TSV/CSV")
    ap.add_argument("--metadata-sample-column", default=defaults.get("metadata_sample_column", "sample"))
    ap.add_argument("--metadata-species-column", default=defaults.get("metadata_species_column", "species"))
    ap.add_argument("--coverage-threshold", type=float, default=float(defaults.get("coverage_threshold") or 100))
    ap.add_argument("--low-hetero", type=float, default=float(defaults.get("low_hetero") or 0.05))
    ap.add_argument("--low-homo", type=float, default=float(defaults.get("low_homo") or 0.95))
    ap.add_argument("--include-filtered", action="store_true", default=parse_bool(defaults.get("include_filtered", False)))
    ap.add_argument("--copy-files", action="store_true", default=parse_bool(defaults.get("copy_files", False)))
    ap.add_argument("--allow-single-cov", action="store_true", default=parse_bool(defaults.get("allow_single_cov", False)))
    ap.add_argument("--mtcn-mt-column", default=defaults.get("mtcn_mt_column"))
    ap.add_argument("--mtcn-nuclear-column", default=defaults.get("mtcn_nuclear_column"))
    ap.add_argument("--mtcn-mtcn-column", default=defaults.get("mtcn_mtcn_column"))
    args = ap.parse_args(argv)
    for attr in ("input_root", "outdir", "metadata"):
        value = getattr(args, attr)
        if value is not None and not isinstance(value, Path):
            setattr(args, attr, Path(value) if str(value).strip() else None)
    if not args.input_root or not args.outdir:
        ap.error("--input-root and --outdir are required, either on the CLI or in --config")
    if args.low_hetero > args.low_homo:
        ap.error("--low-hetero must be <= --low-homo")
    return args


def setup_outputs(outdir: Path) -> Dict[str, Path]:
    dirs = {
        "mtcn": outdir / "collected_mtcn",
        "vcf": outdir / "collected_vcf",
        "cov": outdir / "collected_cov",
        "reports": outdir / "reports",
        "logs": outdir / "logs",
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        filename=dirs["logs"] / "collection_warnings.log",
        level=logging.WARNING,
        format="%(asctime)s\t%(levelname)s\t%(message)s",
    )
    return dirs


def find_one(sample_dir: Path, filename: str) -> Optional[Path]:
    matches = sorted(p for p in sample_dir.rglob(filename) if p.is_file())
    if len(matches) > 1:
        logging.warning("%s: multiple matches found; using %s", filename, matches[0])
    return matches[0] if matches else None


def link_or_copy(src: Path, dest: Path, copy_files: bool) -> None:
    if dest.exists() or dest.is_symlink():
        dest.unlink()
    if copy_files:
        shutil.copy2(src, dest)
    else:
        os.symlink(src.resolve(), dest)


def sniff_delimiter(path: Path) -> str:
    return "," if path.suffix.lower() == ".csv" else "\t"


def load_metadata(path: Optional[Path], sample_col: str, species_col: str) -> Dict[str, str]:
    if not path:
        return {}
    species: Dict[str, str] = {}
    delimiter = sniff_delimiter(path)
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames and sample_col in reader.fieldnames and species_col in reader.fieldnames:
            for row in reader:
                if row.get(sample_col):
                    species[row[sample_col]] = row.get(species_col, "NA") or "NA"
            return species

    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        for row in reader:
            if not row or row[0].lstrip().startswith("#"):
                continue
            if len(row) < 2:
                continue
            if row[0].strip().lower() == sample_col.lower() and row[1].strip().lower() == species_col.lower():
                continue
            sample = row[0].strip()
            if sample:
                species[sample] = row[1].strip() or "NA"
    if species:
        logging.info("Loaded headerless two-column metadata from %s", path)
    else:
        logging.warning("Metadata missing required columns %s and/or %s", sample_col, species_col)
    return species


def alias_lookup(fieldnames: Sequence[str], requested: Optional[str], aliases: Sequence[str]) -> Optional[str]:
    if requested and requested in fieldnames:
        return requested
    normalized = {norm_name(f): f for f in fieldnames}
    if requested and norm_name(requested) in normalized:
        return normalized[norm_name(requested)]
    for alias in aliases:
        if norm_name(alias) in normalized:
            return normalized[norm_name(alias)]
    return None


def parse_mtcn(path: Optional[Path], args: argparse.Namespace) -> Tuple[Dict[str, str], List[str]]:
    values = {"mt_median_coverage": "NA", "nuclear_median_coverage": "NA", "mtcn_median": "NA"}
    notes: List[str] = []
    if not path:
        return values, notes
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            return values, ["mtCN file has no header"]
        row = next(reader, None)
        if row is None:
            return values, ["mtCN file has no data rows"]
        requested = {
            "mt_median_coverage": args.mtcn_mt_column,
            "nuclear_median_coverage": args.mtcn_nuclear_column,
            "mtcn_median": args.mtcn_mtcn_column,
        }
        for key, aliases in MTCN_ALIASES.items():
            col = alias_lookup(reader.fieldnames, requested[key], aliases)
            if col is None:
                notes.append(f"mtCN column not found for {key}")
            else:
                values[key] = row.get(col, "NA") or "NA"
    return values, notes


def read_cov(path: Path) -> Dict[Tuple[str, str, str], float]:
    data: Dict[Tuple[str, str, str], float] = {}
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 4:
                raise ValueError(f"Coverage row has fewer than 4 columns in {path}: {row}")
            if row[0].lower() == "chrom" and row[1].lower() == "pos":
                continue
            data[(row[0], row[1], row[2])] = float(row[3])
    return data


def merge_coverage(paths: Sequence[Path], dest: Path) -> List[float]:
    merged: Dict[Tuple[str, str, str], float] = {}
    for path in paths:
        for key, cov in read_cov(path).items():
            merged[key] = max(merged.get(key, cov), cov)
    def sort_key(item: Tuple[Tuple[str, str, str], float]):
        chrom, pos, target = item[0]
        try:
            pos_key = int(pos)
        except ValueError:
            pos_key = pos
        return chrom, pos_key, target
    coverages: List[float] = []
    with dest.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["chrom", "pos", "target", "coverage"])
        for (chrom, pos, target), cov in sorted(merged.items(), key=sort_key):
            writer.writerow([chrom, pos, target, format_number(cov)])
            coverages.append(cov)
    return coverages


def format_number(value: float) -> str:
    return str(int(value)) if float(value).is_integer() else f"{value:.6g}"


def coverage_metrics(coverages: Sequence[float], threshold: float) -> Tuple[str, str]:
    if not coverages:
        return "NA", "NA"
    pct = 100.0 * sum(c > threshold for c in coverages) / len(coverages)
    med = median(coverages)
    if med == 0:
        mad = "NA"
    else:
        mad = f"{median(abs(c - med) for c in coverages) / med:.6g}"
    return f"{pct:.6g}", mad


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def parse_afs(format_keys: List[str], sample_value: str) -> List[float]:
    values = sample_value.split(":")
    if "AF" not in format_keys:
        return []
    idx = format_keys.index("AF")
    if idx >= len(values) or values[idx] in {".", ""}:
        return []
    afs = []
    for val in values[idx].replace(";", ",").split(","):
        try:
            afs.append(float(val))
        except ValueError:
            continue
    return afs


def count_vcf(path: Optional[Path], low_hetero: float, low_homo: float, include_filtered: bool) -> Tuple[str, str, List[str]]:
    if not path:
        return "NA", "NA", []
    hetero = homo = 0
    notes: List[str] = []
    saw_af = False
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            filt = fields[6]
            if not include_filtered and filt not in {"PASS", "."}:
                continue
            afs = parse_afs(fields[8].split(":"), fields[9])
            if afs:
                saw_af = True
            for af in afs:
                if af >= low_homo:
                    homo += 1
                elif af >= low_hetero:
                    hetero += 1
    if not saw_af:
        notes.append("No AF values found in VCF sample FORMAT fields")
    return str(hetero), str(homo), notes


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    dirs = setup_outputs(args.outdir)
    species_map = load_metadata(args.metadata, args.metadata_sample_column, args.metadata_species_column)
    sample_dirs = sorted(p for p in args.input_root.iterdir() if p.is_dir())
    summary_path = dirs["reports"] / "variant_calling_collection_summary.tsv"

    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for sample_dir in sample_dirs:
            sample = sample_dir.name
            notes: List[str] = []
            missing: List[str] = []
            mtcn = find_one(sample_dir, f"{sample}.round2.mtcn.tsv")
            vcf = find_one(sample_dir, f"{sample}.round2.original_coords.clean.final.split.vcf.gz") or find_one(sample_dir, f"{sample}.round2.original_coords.clean.final.split.vcf")
            decoy_cov = find_one(sample_dir, f"{sample}.numt_decoy.clean.realigned.per_base_coverage.tsv")
            round2_cov = find_one(sample_dir, f"{sample}.round2.original_coords.per_base_coverage.tsv")

            if not mtcn: missing.append("mtcn")
            if not vcf: missing.append("vcf")
            cov_paths = [p for p in [decoy_cov, round2_cov] if p]
            if len(cov_paths) < 2 and not args.allow_single_cov:
                if not decoy_cov: missing.append("decoy_coverage")
                if not round2_cov: missing.append("round2_original_coords_coverage")
            elif len(cov_paths) == 1:
                notes.append("single coverage file used because --allow-single-cov is set")

            mtcn_dest = dirs["mtcn"] / f"{sample}.round2.mtcn.tsv"
            vcf_suffix = ".vcf.gz" if vcf and vcf.name.endswith(".vcf.gz") else ".vcf"
            vcf_dest = dirs["vcf"] / f"{sample}.round2.original_coords.clean.final.split{vcf_suffix}"
            cov_dest = dirs["cov"] / f"{sample}.merged.max_depth.per_base_coverage.tsv"
            mtcn_values, mtcn_notes = parse_mtcn(mtcn, args)
            notes.extend(mtcn_notes)
            n_hetero, n_homo, vcf_notes = count_vcf(vcf, args.low_hetero, args.low_homo, args.include_filtered)
            notes.extend(vcf_notes)
            pct100 = mad = "NA"

            status = "PASS_COLLECTION" if not missing else "FAIL_MISSING_INPUT"
            try:
                if mtcn and not missing:
                    link_or_copy(mtcn, mtcn_dest, args.copy_files)
                elif mtcn:
                    link_or_copy(mtcn, mtcn_dest, args.copy_files)
                if vcf and not missing:
                    link_or_copy(vcf, vcf_dest, args.copy_files)
                elif vcf:
                    link_or_copy(vcf, vcf_dest, args.copy_files)
                if cov_paths and (len(cov_paths) == 2 or args.allow_single_cov):
                    coverages = merge_coverage(cov_paths, cov_dest)
                    pct100, mad = coverage_metrics(coverages, args.coverage_threshold)
            except Exception as exc:  # sample-level failure should not stop all samples
                status = "FAIL_PROCESSING"
                notes.append(f"processing error: {exc}")
                logging.warning("%s: processing error: %s", sample, exc)

            if missing:
                logging.warning("%s: missing required files: %s", sample, ",".join(missing))
            for note in notes:
                logging.warning("%s: %s", sample, note)

            writer.writerow({
                "sample": sample,
                "species": species_map.get(sample, "NA"),
                **mtcn_values,
                "Percent_100": pct100,
                "MAD": mad,
                "n_hetero": n_hetero,
                "n_homo": n_homo,
                "vcf_file": str(vcf_dest if vcf_dest.exists() or vcf_dest.is_symlink() else "NA"),
                "cov_file": str(cov_dest if cov_dest.exists() else "NA"),
                "mtcn_file": str(mtcn_dest if mtcn_dest.exists() or mtcn_dest.is_symlink() else "NA"),
                "status": status,
                "missing_files": ",".join(missing) if missing else "NA",
                "notes": "; ".join(notes) if notes else "NA",
            })
    return 0


if __name__ == "__main__":
    sys.exit(main())
