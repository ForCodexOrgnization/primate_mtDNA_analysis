#!/usr/bin/env python3
"""Build immutable, race-safe task manifests for QC Slurm arrays."""
import argparse, csv, datetime, os, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from qc_analysis.lib.simple_yaml import read_simple_yaml

SAMPLE_STEPS = {"coordinate_liftover", "codon_match", "trna_match", "rrna_match"}
GLOBAL_STEPS = {"collect_variant_calling_results", "discover_global_anchor", "build_primate_codon_table",
 "compare_genbank_mitos2", "mitos2_prepare_tasks", "mitos2_merge", "codon_match_validate",
 "codon_match_merge", "intraspecies_contamination"}

def table_samples(path):
    path=Path(path)
    if not path.is_file(): return []
    with path.open(newline='') as f:
        rows=list(csv.reader(f, delimiter='\t'))
    if not rows:return []
    header=[x.strip().lower() for x in rows[0]]
    start=1 if any(x in header for x in ('sample','sample_id','name')) else 0
    col=next((header.index(x) for x in ('sample','sample_id','name') if x in header),0)
    return [r[col].strip() for r in rows[start:] if len(r)>col and r[col].strip()]

def valid_vcf(path, tag):
    try:
        if Path(path).stat().st_size == 0:return False
        text=Path(path).read_text(errors='ignore')
        return '#CHROM' in text and tag in text
    except OSError:return False

def paths_for(step,s,cfg):
    if step=='coordinate_liftover':
        p=cfg[step]['paths']; return '',str(Path(p['output_dir'])/'vcf_lifted_raw'/f'{s}.lifted.raw.vcf'),'MTCODON'
    sec=cfg[step];p=sec['paths']; st=sec['settings']
    if step=='codon_match': inp=Path(p['input_vcf_dir'])/st['input_vcf_pattern'].format(sample=s); tag='MTCODON'
    elif step=='trna_match':
        a=Path(p['input_vcf_dir'])/st['input_vcf_pattern'].format(sample=s); b=Path(p['fallback_input_vcf_dir'])/st['fallback_input_vcf_pattern'].format(sample=s);inp=a if a.exists() else b;tag='MTTRNA'
    else:
        choices=[Path(p['input_vcf_dir'])/st['input_vcf_pattern'].format(sample=s),Path(p['fallback_codon_vcf_dir'])/st['fallback_codon_vcf_pattern'].format(sample=s),Path(p['fallback_raw_vcf_dir'])/st['fallback_raw_vcf_pattern'].format(sample=s)];inp=next((x for x in choices if x.exists()),choices[0]);tag='MTRRNA'
    folder={'codon_match':'vcf_codon','trna_match':'vcf_trna','rrna_match':'vcf_rrna'}[step]
    suffix=st['output_suffix'] if step!='trna_match' or str(inp).startswith(str(p['input_vcf_dir'])) else '.lifted.trna.vcf'
    return str(inp),str(Path(p['output_dir'])/folder/f'{s}{suffix}'),tag

def main():
    ap=argparse.ArgumentParser();ap.add_argument('step');ap.add_argument('config');ap.add_argument('--sample');ap.add_argument('--outdir',default='results/qc/job_arrays');ap.add_argument('--force',action='store_true');ap.add_argument('--retry',action='store_true');a=ap.parse_args()
    cfg=read_simple_yaml(Path(a.config)); step=a.step
    if a.sample: candidates=[a.sample]
    elif step in GLOBAL_STEPS: candidates=[step]
    elif step=='mitos2_annotation':
        task=Path(cfg[step]['paths'].get('mitos2_reference_tasks',''))
        if task.is_file():
            with task.open() as f:candidates=['task:'+r['task_id'] for r in csv.DictReader(f,delimiter='\t') if r.get('status')!='completed']
        else:
            # This also makes a dry-run/all graph possible before its dependent
            # prepare job runs; --reference is a supported, race-safe worker key.
            sample_file=cfg[step]['paths'].get('sample_ref_file','')
            candidates=['reference:'+x for x in table_samples(sample_file)]
    elif step in SAMPLE_STEPS:
        if step=='coordinate_liftover': source=cfg[step]['paths']['sample_ref_file']
        elif step=='codon_match': source=cfg[step]['paths']['sample_reference_map']
        else: source=cfg.get('coordinate_liftover',{}).get('paths',{}).get('sample_ref_file','')
        candidates=table_samples(source)
    else: raise SystemExit(f'Unsupported array step: {step}')
    candidates=sorted(set(x for x in candidates if x))
    rows=[]; done=missing=invalid=0
    for item in candidates:
        inp=out='';tag=''
        if step in SAMPLE_STEPS: inp,out,tag=paths_for(step,item,cfg)
        if inp and not Path(inp).exists(): missing+=1;continue
        complete=bool(out and valid_vcf(out,tag))
        if out and Path(out).exists() and not complete:invalid+=1
        if complete: done+=1
        include=(not complete or a.force) if not a.retry else (not complete)
        if include:rows.append((item,inp,out,'force_rerun' if complete else 'pending'))
    if not rows: raise SystemExit(f'ERROR: no eligible tasks for {step} (candidates={len(candidates)}, completed={done}, missing_inputs={missing})')
    outdir=Path(a.outdir);outdir.mkdir(parents=True,exist_ok=True);stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    task=outdir/f'{step}.{stamp}.tasks.txt'; manifest=outdir/f'{step}.{stamp}.manifest.tsv'
    tmp=task.with_suffix(task.suffix+'.tmp');tmp.write_text(''.join(r[0]+'\n' for r in rows));os.replace(tmp,task)
    mt=manifest.with_suffix(manifest.suffix+'.tmp')
    with mt.open('w',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(('task_index','item_type','item','expected_input','expected_output','status_at_submission','config','submitted_at'))
        typ='sample' if step in SAMPLE_STEPS else 'reference' if step=='mitos2_annotation' else 'singleton'
        for i,r in enumerate(rows,1):w.writerow((i,typ,*r,a.config,stamp))
    os.replace(mt,manifest)
    pointer=outdir/f'{step}.current.tsv';pt=pointer.with_suffix('.tmp');pt.write_text(str(manifest)+'\n');os.replace(pt,pointer)
    print(f'TASK_FILE={task}');print(f'MANIFEST={manifest}');print(f'COUNT={len(rows)}')
    print(f'STATS=total candidates {len(candidates)}; already completed {done}; scheduled {len(rows)}; missing inputs {missing}; invalid existing outputs {invalid}',file=sys.stderr)
if __name__=='__main__':main()
