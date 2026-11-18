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

def choose_source_record(lookup, key, info, fields, src):
    candidates = lookup.execute(
        'SELECT qc,dp,klass,af,source_filter FROM v WHERE sample=? AND chrom=? AND pos=? AND ref=? AND alt=?', key
    ).fetchall()
    if not candidates:
        raise RuntimeError('No native source record for %r in %s' % (key,src))
    fmt = dict(zip(fields[8].split(':'),fields[9].split(':'))) if len(fields)>9 else {}
    dp = info.get('SOURCE_DP') or fmt.get('DP') or info.get('DP')
    af = info.get('SOURCE_AF') or fmt.get('AF') or info.get('AF')
    klass = info.get('SOURCE_CALL_CLASS')
    source_filter = info.get('SOURCE_FILTER')
    if dp not in (None,'.',''):
        candidates = [r for r in candidates if same_numeric(r[1],dp,tol=0)]
    if af not in (None,'.',''):
        status = info.get('LIFTOVER_ALLELE_STATUS')
        if status in ('ALT_REF_FLIP', 'REF_ALT_FLIP'):
            flip = True
        elif status in ('REF_MATCH', 'NO_FLIP', None, '', '.'):
            # Verify orientation with the immutable native alleles if
            # the original VCF did not annotate the flip status.
            source_ref, source_alt = key[3], key[4]
            target_ref, target_alt = fields[3], fields[4]
            if (target_ref, target_alt) == (source_alt, source_ref):
                flip = True
            elif (target_ref, target_alt) == (source_ref, source_alt):
                flip = False
            else:
                raise RuntimeError(
                    'Unverifiable REF/ALT orientation for %r in %s: '
                    'target=%s>%s status=%s'
                    % (key, src, target_ref, target_alt, status)
                )
        else:
            raise RuntimeError(
                'Unexpected liftover allele status for %r in %s: %s'
                % (key, src, status)
            )
        def af_match(row):
            try:
                expected = 1.0-float(af) if flip else float(af)
                # Source QC AF is formatted with ~3 significant decimals
                # (e.g. 0.989), whereas the lifted FORMAT AF has more
                # precision (e.g. 0.0111397 after ALT_REF_FLIP).
                # Respect the source value's displayed precision.
                original = str(row[3])
                from decimal import Decimal
                source_value = Decimal(original)
                half_unit = float(Decimal('0.5').scaleb(source_value.as_tuple().exponent))
                return abs(float(original)-expected) <= half_unit + 1e-8
            except (ValueError,TypeError):
                return False
        candidates = [r for r in candidates if af_match(r)]
    if klass not in (None,'.',''):
        candidates = [r for r in candidates if r[2] == klass]
    if source_filter not in (None,'.',''):
        candidates = [r for r in candidates if r[4] == source_filter]
    if not candidates:
        raise RuntimeError('No matching native DP/AF evidence: %r in %s DP=%s AF=%s' % (key,src,dp,af))
    unique = {(r[0],r[1],r[2],r[3]) for r in candidates}
    if len(unique) != 1:
        raise RuntimeError('Ambiguous native source records %r in %s: %r' % (key,src,candidates))
    return candidates[0]

def op(path, mode):
    return gzip.open(path,mode+'t') if path.suffix=='.gz' else path.open(mode,encoding='utf-8')

def refresh_vcf(src, dst, sample, lookup):
    # NEVER open dst for writing: it can be a symlink to original DP100 VCF.
    import tempfile
    tags = (
        ('SOURCE_CHROM','String'), ('SOURCE_POS','Integer'),
        ('SOURCE_REF','String'), ('SOURCE_ALT','String'),
        ('SOURCE_AF','Float'), ('SOURCE_DP','Integer'),
        ('SOURCE_CALL_CLASS','String'), ('SOURCE_VARIANT_QC','String')
    )
    count = 0
    dst.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=dst.parent,prefix='.dp50_refresh_',
                                 suffix='.vcf.gz' if dst.name.endswith('.gz') else '.vcf',
                                 delete=False) as fh:
        tmp = Path(fh.name)
    try:
        declared = set()
        with op(src,'r') as fin,op(tmp,'w') as fout:
            for line in fin:
                if line.startswith('##INFO=<ID='):
                    declared.add(line.split('ID=',1)[1].split(',',1)[0])
                if line.startswith('#CHROM'):
                    for tag,typ in tags:
                        if tag not in declared:
                            fout.write('##INFO=<ID=%s,Number=1,Type=%s,Description="Native-call provenance for DP50 QC">\n' % (tag,typ))
                    fout.write(line)
                    continue
                if line.startswith('#'):
                    fout.write(line)
                    continue
                f = line.rstrip('\n').split('\t')
                if len(f)<8:
                    raise RuntimeError('Malformed VCF record: %s' % src)
                tokens = f[7].split(';')
                info = dict(t.split('=',1) for t in tokens if '=' in t)
                source_keys = ('CHROM','POS','REF','ALT')
                native = tuple(info.get('SOURCE_'+k) or info.get('SRC_'+k) for k in source_keys)
                if any(v in (None,'.','') for v in native):
                    raise RuntimeError('Missing SOURCE/SRC native identity for sample %s in %s' % (sample,src))
                key = (sample, *native)
                qc,dp,klass,af,_ = choose_source_record(lookup,key,info,f,src)
                newinfo = dict(zip(('SOURCE_CHROM','SOURCE_POS','SOURCE_REF','SOURCE_ALT'),native))
                newinfo.update(SOURCE_AF=af,SOURCE_DP=dp,SOURCE_CALL_CLASS=klass,SOURCE_VARIANT_QC=qc)
                for k,v in newinfo.items():
                    old = info.get(k)
                    if k in ('SOURCE_CHROM','SOURCE_POS','SOURCE_REF','SOURCE_ALT') and old not in (None,'.',v):
                        raise RuntimeError('Conflicting %s for %r in %s: %s vs %s' % (k,key,src,old,v))
                cleaned = [t for t in tokens if t.split('=',1)[0] not in newinfo and t!='.']
                f[7] = ';'.join(cleaned+['%s=%s'%(k,v) for k,v in newinfo.items()])
                fout.write('\t'.join(f)+'\n')
                count+=1
        os.replace(tmp,dst)
    finally:
        tmp.unlink(missing_ok=True)
    return count

def refresh_all(old_root,new_root,stage='rrna_match',samples=None):
    if stage not in VCFS:
        raise RuntimeError('Unsupported refresh stage %s'%stage)
    lookup = source_lookup(new_root/'pre_liftover_variant_qc',new_root)
    folder = VCFS[stage]
    srcdir = old_root/stage/folder
    dstdir = new_root/stage/folder
    if not srcdir.is_dir():
        raise RuntimeError('Missing original annotated VCF directory: '+str(srcdir))
    count=0
    rows=0
    try:
        for src in sorted(srcdir.iterdir()):
            if not src.is_file() or not src.name.endswith(('.vcf','.vcf.gz')):
                continue
            sample=src.name.split('.',1)[0]
            if samples is not None and sample not in samples:
                continue
            dest=dstdir/src.name
            rows+=refresh_vcf(src,dest,sample,lookup)
            count+=1
            if count%100==0:
                print('Refreshed %d VCFs, %d variants'%(count,rows),flush=True)
    finally:
        lookup.close()
    if not count:
        raise RuntimeError('No annotated VCF files refreshed')
    print('Refreshed %d %s VCFs and %d records'%(count,stage,rows),flush=True)
    return count

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
    p.add_argument("--sample", help="Refresh only one sample (validation before full cohort)")
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
            refresh_all(old, new, stage="rrna_match", samples={a.sample} if a.sample else None)
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
            for key in ("exclude_duplicate_source_variants", "require_source_call_class"):
                if "  " + key + ":" not in config_text:
                    config_text = config_text[:section_start] + "  " + key + ": true\n" + config_text[section_start:]
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
