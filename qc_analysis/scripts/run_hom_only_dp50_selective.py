#!/usr/bin/env python3
"""Isolated DP50 HOM-cohort selective rerun. Run from repository root.

Recomputes native SOURCE QC, local artifacts, and terminal QC; reuses coverage,
liftover and functional annotations after strict native-key SOURCE_QC refresh.
Does not edit the original qc_coverage_only output tree.
"""
from __future__ import annotations

import argparse
import csv
import math
import sqlite3
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

def source_lookup(prelift, work_root):
    """Index records by native identity, allowing multiple depth/AF observations."""
    source = prelift / 'reports/source_variant_qc.tsv'
    if not source.is_file():
        raise RuntimeError('Missing merged source QC: ' + str(source))
    database = work_root / 'dp50_source_lookup.sqlite'
    if database.exists():
        database.unlink()
    db = sqlite3.connect(str(database))
    db.execute('CREATE TABLE v (sample TEXT, chrom TEXT, pos TEXT, ref TEXT, alt TEXT, qc TEXT, dp TEXT, klass TEXT, af TEXT, source_filter TEXT)')
    rows = []
    count = 0
    with source.open(newline='') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            rows.append(tuple(row[k] for k in ('sample','source_chrom','source_pos','source_ref','source_alt','source_variant_qc','source_dp','source_call_class','source_af','source_filter')))
            if len(rows) >= 10000:
                db.executemany('INSERT INTO v VALUES (?,?,?,?,?,?,?,?,?,?)', rows)
                count += len(rows)
                rows.clear()
        if rows:
            db.executemany('INSERT INTO v VALUES (?,?,?,?,?,?,?,?,?,?)', rows)
            count += len(rows)
    db.execute('CREATE INDEX v_native_key ON v(sample,chrom,pos,ref,alt)')
    db.commit()
    if not count:
        raise RuntimeError('Empty source QC report')
    print('Indexed %d source VCF rows; duplicated native identities are retained for DP/AF matching' % count, flush=True)
    return db

def same_numeric(a, b, tol=1e-8):
    try:
        x, y = float(a), float(b)
        return math.isfinite(x) and math.isfinite(y) and abs(x - y) <= tol
    except (ValueError, TypeError):
        return False

def choose_source_record(lookup, key, info, src):
    candidates = lookup.execute(
        'SELECT qc,dp,klass,af,source_filter FROM v WHERE sample=? AND chrom=? AND pos=? AND ref=? AND alt=?',
        key,
    ).fetchall()
    if not candidates:
        raise RuntimeError('Missing native source record %r in %s' % (key,src))
    dp = info.get('SOURCE_DP')
    af = info.get('SOURCE_AF')
    klass = info.get('SOURCE_CALL_CLASS')
    source_filter = info.get('SOURCE_FILTER')
    # When two records share a position, depth and source AF identify
    # the exact native observation. Never choose the first by position alone.
    if dp not in (None, '.', ''):
        candidates = [row for row in candidates if same_numeric(row[1], dp, tol=0.0)]
    if af not in (None, '.', ''):
        candidates = [row for row in candidates if same_numeric(row[3], af)]
    if klass not in (None, '.', ''):
        candidates = [row for row in candidates if row[2] == klass]
    if source_filter not in (None, '.', ''):
        candidates = [row for row in candidates if row[4] == source_filter]
    if not candidates:
        raise RuntimeError('No source DP/AF/class match for %r in %s (DP=%s, AF=%s, class=%s)' % (key,src,dp,af,klass))
    decisions = {row[0] for row in candidates}
    if len(decisions) != 1:
        raise RuntimeError('Ambiguous source QC decisions for %r in %s: %r' % (key,src,candidates))
    # Multiple indistinguishable records carrying the same decision are safe
    # for refreshing QC only; the original VCF records are never deduplicated.
    return candidates[0][0]

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
            fresh_qc = choose_source_record(lookup, key, info, src)
            tokens = [t for t in tokens if not t.startswith("SOURCE_VARIANT_QC=")]
            tokens.append("SOURCE_VARIANT_QC=" + fresh_qc)
            fields[7] = ";".join(tokens)
            fout.write("\t".join(fields) + "\n")
            count += 1
    return count

def refresh_all(old_root, new_root):
    lookup = source_lookup(new_root / "pre_liftover_variant_qc", new_root)
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
    lookup.close()
    if not total_files:
        raise RuntimeError("No annotated VCF files were refreshed")
    print(f"Refreshed {total_variants} VCF rows across {total_files} files")
    return total_files

def reuse_annotated_vcfs(old_root, new_root):
    """Symlink original annotation VCFs, preserving all old files.

    Existing SOURCE_VARIANT_QC may still reflect DP100. Current final filter
    uses SOURCE_DP and VCF FILTER directly; use refresh if updated provenance
    tags are needed in exported VCFs.
    """
    counts = {}
    for stage, folder in VCFS.items():
        srcdir = old_root / stage / folder
        dstdir = new_root / stage / folder
        if not srcdir.is_dir() or not dstdir.is_dir():
            raise RuntimeError('Missing source/staged VCF directory: %s %s' % (srcdir, dstdir))
        created = 0
        existing = 0
        for src in sorted(srcdir.iterdir()):
            if not src.is_file() or not src.name.endswith(('.vcf', '.vcf.gz')):
                continue
            dst = dstdir / src.name
            if dst.is_symlink():
                if dst.resolve() != src.resolve():
                    raise RuntimeError('Refusing conflicting symlink: '+str(dst))
                existing += 1
            elif dst.exists():
                if dst.stat().st_size == 0:
                    raise RuntimeError('Refusing empty staged VCF: '+str(dst))
                existing += 1
            else:
                dst.symlink_to(src.resolve())
                created += 1
        counts[stage] = created + existing
        print('%s: %d VCFs available (%d links created; %d existing)' % (stage,counts[stage],created,existing),flush=True)
    if not counts.get('rrna_match'):
        raise RuntimeError('No final annotated rRNA VCFs available')
    print('NOTE: Existing SOURCE_VARIANT_QC tags are historical; final filter recomputes DP gate from SOURCE_DP.',flush=True)
    return counts

def diagnose_final_inputs(config_path):
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))
    from qc_analysis.lib.simple_yaml import read_simple_yaml
    cfg = read_simple_yaml(config_path)
    sec = cfg.get('final_filter') or {}
    qc_spec = (sec.get('sample_reports') or {}).get('sample_qc') or {}
    qc_path = ROOT / str(qc_spec.get('path',''))
    if not qc_path.is_file():
        raise RuntimeError('Missing sample QC report: ' + str(qc_path))
    with qc_path.open(newline='') as handle:
        samples = [r['sample'] for r in csv.DictReader(handle,delimiter='\t') if r.get('homoplasmy_qc_status') == 'PASS']
    sources = sec.get('vcf_sources') or {}
    sample_covered = set()
    print('HOM eligible samples:',len(samples),flush=True)
    for name,spec in sources.items():
        directory = ROOT / str(spec['dir'])
        pattern = str(spec.get('pattern','{sample}.vcf'))
        found = []
        for sample in samples:
            candidate = directory / pattern.format(sample=sample)
            if candidate.is_file() or Path(str(candidate)+'.gz').is_file():
                found.append(sample)
        sample_covered.update(found)
        print(name, len(found), '/', len(samples), 'directory=',directory,flush=True)
        print('  actual filenames:', [p.name for p in sorted(directory.iterdir())[:5]] if directory.is_dir() else 'MISSING DIRECTORY',flush=True)
    print('Eligible samples with any VCF:',len(sample_covered),'/',len(samples),flush=True)
    for sample in samples[:5]:
        print('EXAMPLE EXPECTED:',sample,[str(spec['dir'])+'/'+str(spec['pattern']).format(sample=sample) for spec in sources.values()],flush=True)
    return len(sample_covered),len(samples)

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--execute", action="store_true", help="Legacy sequential execution")
    p.add_argument("--phase", choices=["prepare","refresh","reuse","diagnose","final","hom"], help="Slurm stage")
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
        elif a.phase == "reuse":
            reuse_annotated_vcfs(old, new)
            diagnose_final_inputs(staged)
        elif a.phase == "diagnose":
            diagnose_final_inputs(staged)
        elif a.phase == "final":
            found, eligible = diagnose_final_inputs(staged)
            if eligible and found < int(eligible * 0.95):
                raise RuntimeError("Only %d/%d HOM-eligible samples have annotated VCF; refusing incomplete final QC" % (found, eligible))
            # Only alter the isolated DP50 staged YAML, never the original QC config.
            config_text = staged.read_text()
            if "final_filter:\n" not in config_text:
                raise RuntimeError("Staged YAML missing final_filter section")
            section_start = config_text.index("final_filter:\n") + len("final_filter:\n")
            existing = "  exclude_duplicate_source_variants: true\n"
            if "  exclude_duplicate_source_variants:" not in config_text:
                config_text = config_text[:section_start] + existing + config_text[section_start:]
            staged.write_text(config_text)
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
