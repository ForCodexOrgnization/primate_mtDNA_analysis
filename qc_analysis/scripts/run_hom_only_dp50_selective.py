#!/usr/bin/env python3
"""Isolated DP50 HOM-cohort selective rerun. Run from repository root.

Recomputes native SOURCE QC, local artifacts, and terminal QC; reuses coverage,
liftover and functional annotations after strict native-key SOURCE_QC refresh.
Does not edit the original qc_coverage_only output tree.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OLD = "results/qc_coverage_only"
NEW = "results/qc_coverage_only_dp50"
CONFIG = "config/qc_preprocessing_coverage_hom_only.yaml"
# Directories generated afresh, never linked to the historical output.
REBUILT = {"pre_liftover_variant_qc", "local_heteroplasmy_qc",
           "final_filter_hom_only", "homoplasmy_qc_only"}
# Annotated inputs to final_filter in priority order, plus their fallbacks.
VCFS = {
    "coordinate_liftover": "vcf_lifted_raw",
    "codon_match": "vcf_codon",
    "trna_match": "vcf_trna",
    "rrna_match": "vcf_rrna",
}

def run(*args):
    print("+", " ".join(str(x) for x in args), flush=True)
    subprocess.run([str(x) for x in args], cwd=ROOT, check=True)

def link_children(source, target, excluded=()):
    target.mkdir(parents=True, exist_ok=True)
    for p in source.iterdir():
        if p.name in excluded:
            continue
        q = target / p.name
        if q.exists() or q.is_symlink():
            raise RuntimeError(f"Refusing existing staged file: {q}")
        q.symlink_to(p.resolve(), target_is_directory=p.is_dir())

def initialize(old_root, new_root, source_config, staged_config):
    if new_root.exists() or new_root.is_symlink():
        raise RuntimeError(f"Destination already exists; refusing overwrite: {new_root}")
    if not old_root.is_dir():
        raise RuntimeError(f"Old result root missing: {old_root}")
    source = source_config.read_text()
    if OLD not in source:
        raise RuntimeError("Expected original coverage-only paths not found in YAML")
    if "variant_dp_min: 50" not in source:
        raise RuntimeError("Source config must explicitly specify DP >= 50")
    if staged_config.exists():
        raise RuntimeError(f"Staged config already exists: {staged_config}")
    new_root.mkdir(parents=True)
    for p in old_root.iterdir():
        q = new_root / p.name
        if p.name in REBUILT:
            q.mkdir()
        elif p.name in VCFS:
            link_children(p, q, excluded={VCFS[p.name]})
            (q / VCFS[p.name]).mkdir()
        else:
            q.symlink_to(p.resolve(), target_is_directory=p.is_dir())
    staged_config.write_text(source.replace(OLD, NEW))
    print("Created isolated tree:", new_root)
    print("Generated config:", staged_config)

def source_lookup(prelift):
    path = prelift / "reports/source_variant_qc.tsv"
    if not path.is_file():
        raise RuntimeError(f"Missing merged source QC: {path}")
    lookup = {}
    with path.open(newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            key = tuple(row[k] for k in
                        ("sample", "source_chrom", "source_pos", "source_ref", "source_alt"))
            if key in lookup:
                raise RuntimeError(f"Duplicate native key: {key}")
            lookup[key] = row
    if not lookup:
        raise RuntimeError("Empty source QC report; stopping")
    return lookup

def op(path, mode):
    return gzip.open(path, mode + "t") if path.suffix == ".gz" else path.open(mode, encoding="utf-8")

def refresh_vcf(src, dst, sample, lookup):
    count = 0
    with op(src, "r") as fin, op(dst, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            tokens = fields[7].split(";")
            info = {t.split("=", 1)[0]: t.split("=", 1)[1] for t in tokens if "=" in t}
            source_names = ("SOURCE_CHROM", "SOURCE_POS", "SOURCE_REF", "SOURCE_ALT")
            if any(not info.get(k) or info[k] == "." for k in source_names):
                raise RuntimeError(f"Unverifiable native key in {src}, line {count + 1}")
            key = (sample, *(info[k] for k in source_names))
            if key not in lookup:
                raise RuntimeError(f"No matching native variant for {key} in {src}")
            fresh = lookup[key]
            # Original DP/AF/class must agree; their biological meaning does NOT
            # change with the QC threshold. Only SOURCE_VARIANT_QC is refreshed.
            for name, column in (("SOURCE_CALL_CLASS", "source_call_class"),
                                 ("SOURCE_DP", "source_dp")):
                if info.get(name) not in (None, fresh[column]):
                    try:
                        same = float(info[name]) == float(fresh[column]) if name == "SOURCE_DP" else False
                    except (ValueError, TypeError):
                        same = False
                    if not same:
                        raise RuntimeError(f"Conflicting {name} for {key} in {src}")
            tokens = [t for t in tokens if not t.startswith("SOURCE_VARIANT_QC=")]
            tokens.append("SOURCE_VARIANT_QC=" + fresh["source_variant_qc"])
            fields[7] = ";".join(tokens)
            fout.write("\t".join(fields) + "\n")
            count += 1
    return count

def refresh_all(old_root, new_root):
    lookup = source_lookup(new_root / "pre_liftover_variant_qc")
    total_files, total_variants = 0, 0
    for stage, folder in VCFS.items():
        srcdir = old_root / stage / folder
        dstdir = new_root / stage / folder
        if not srcdir.is_dir():
            print("No VCF directory, skipping:", srcdir)
            continue
        for src in sorted(srcdir.iterdir()):
            if not src.is_file() or (src.suffix != ".vcf" and not src.name.endswith(".vcf.gz")):
                continue
            sample = src.name.split(".", 1)[0]
            dest = dstdir / src.name
            total_variants += refresh_vcf(src, dest, sample, lookup)
            total_files += 1
    if not total_files:
        raise RuntimeError("No annotated VCF files were refreshed")
    print(f"Refreshed {total_variants} VCF rows across {total_files} files")
    return total_files

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--execute", action="store_true", help="Legacy sequential execution")
    p.add_argument("--phase", choices=["prepare","refresh","final","hom"], help="Slurm stage")
    p.add_argument("--old-root", default=OLD)
    p.add_argument("--new-root", default=NEW)
    p.add_argument("--config", default=CONFIG)
    p.add_argument("--staged-config", default="config/qc_preprocessing_coverage_hom_only_dp50_selective.yaml")
    a = p.parse_args()
    def resolve(x):
        path = Path(x).expanduser()
        return path if path.is_absolute() else ROOT / path
    old, new, source, staged = map(resolve, (a.old_root, a.new_root, a.config, a.staged_config))
    # Path substitution is deliberately tied to the default internal layout.
    if a.old_root != OLD or a.new_root != NEW:
        p.error("Custom roots not supported: source YAML substitution requires the standard paths")
    print("PLAN: stage independent DP50 tree; recompute source QC and local artifacts;")
    print("      refresh immutable-key SOURCE_VARIANT_QC in annotated VCF copies;")
    print("      rerun final_filter and homoplasmy_qc; reuse all original read data")
    print("      and original reference/functional annotations.")
    if a.phase == "prepare":
        if not new.exists() and not staged.exists():
            initialize(old, new, source, staged)
        elif new.is_dir() and staged.is_file():
            print("Resuming existing workspace; preserving partial files")
        else:
            raise RuntimeError("Partial staging requires manual inspection")
        return
    if a.phase:
        if not new.is_dir() or not staged.is_file():
            raise RuntimeError("Missing staged workspace")
        if a.phase == "refresh":
            refresh_all(old, new)
        elif a.phase == "final":
            run(sys.executable, ROOT / "qc_analysis/scripts/run_final_filter_with_heteroplasmy.py", "--config", staged)
        elif a.phase == "hom":
            run(sys.executable, ROOT / "qc_analysis/scripts/run_homoplasmy_qc.py", "--config", staged)
        return
    if a.execute:
        raise RuntimeError("Sequential execution disabled; use Slurm submit script")
    print("DRY RUN ONLY. Submit via submit_hom_only_dp50_selective.sh")

if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as exc:
        print("ERROR:", exc, file=sys.stderr)
        sys.exit(2)
