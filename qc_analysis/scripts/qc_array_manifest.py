#!/usr/bin/env python3
"""Build immutable, race-safe task manifests for QC Slurm arrays."""
import argparse, csv, datetime, os, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from qc_analysis.lib.simple_yaml import read_simple_yaml

SAMPLE_STEPS = {"coordinate_liftover", "codon_match", "trna_match", "rrna_match"}
GLOBAL_STEPS = {"collect_variant_calling_results", "discover_global_anchor", "build_primate_codon_table",
 "compare_genbank_mitos2", "mitos2_prepare_tasks", "mitos2_merge", "codon_match_validate",
 "codon_match_merge", "build_trna_indexes", "trna_match_merge", "trna_gene_qc", "build_primate_homo_background",
 "intraspecies_contamination", "sample_variant_filtering", "human_contamination", "final_filter"}

STEP_SECTIONS = {
    "collect_variant_calling_results": "collect_variant_calling",
    "discover_global_anchor": "global_anchor_discovery",
    "coordinate_liftover": "coordinate_liftover",
    "human_contamination": "human_contamination",
    "mitos2_prepare_tasks": "mitos2_annotation", "mitos2_annotation": "mitos2_annotation",
    "mitos2_merge": "mitos2_annotation", "build_primate_codon_table": "build_primate_codon_table",
    "compare_genbank_mitos2": "genbank_mitos2_comparison",
    "codon_match_validate": "codon_match", "codon_match": "codon_match", "codon_match_merge": "codon_match",
    "build_trna_indexes": "trna_match", "trna_match": "trna_match",
    "trna_match_merge": "trna_match", "rrna_match": "rrna_match",
    "intraspecies_contamination": "intraspecies_contamination",
    "sample_variant_filtering": "sample_variant_filtering",
    "final_filter": "final_filter",
}
FALLBACK_OUTPUTS = {
    "collect_variant_calling_results": "results/qc/variant_calling_collection",
    "discover_global_anchor": "results/qc/coordinate_liftover/global_anchor",
    "coordinate_liftover": "results/qc/coordinate_liftover",
    "human_contamination": "results/qc/human_contamination",
    "mitos2_prepare_tasks": "results/qc/mitos2_annotation", "mitos2_annotation": "results/qc/mitos2_annotation",
    "mitos2_merge": "results/qc/mitos2_annotation", "build_primate_codon_table": "results/qc/codon_table_build",
    "compare_genbank_mitos2": "results/qc/genbank_mitos2_comparison",
    "codon_match_validate": "results/qc/codon_match", "codon_match": "results/qc/codon_match",
    "codon_match_merge": "results/qc/codon_match", "build_trna_indexes": "results/qc/trna_match",
    "trna_match": "results/qc/trna_match", "trna_match_merge": "results/qc/trna_match",
    "rrna_match": "results/qc/rrna_match", "intraspecies_contamination": "results/qc/intraspecies_contamination",
    "sample_variant_filtering": "results/qc/sample_variant_filtering",
    "final_filter": "results/qc/final_filter",
}

def resolve_runtime_paths(step, cfg):
    """Resolve biological output, scheduler metadata, and log roots for a step."""
    section=cfg.get(STEP_SECTIONS.get(step,''),{}) or {}; paths=section.get('paths',{}) or {}
    explicit=paths.get('job_array_dir')
    output=paths.get('output_dir')
    if not output:
        reports=paths.get('reports_dir')
        if reports: output=str(Path(reports).parent)
    configured_output=bool(output)
    output=str(output or FALLBACK_OUTPUTS.get(step, f'.workflow/qc_preprocessing/{step}'))
    fallback_metadata=step not in FALLBACK_OUTPUTS and not configured_output
    job_array=str(explicit or (Path(output) if fallback_metadata else Path(output)/'job_arrays'))
    log_dir=paths.get('log_dir') or str(Path(output)/'logs'/'job_arrays')
    return {'output_dir':output, 'job_array_dir':job_array, 'log_dir':str(log_dir)}

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
    ap=argparse.ArgumentParser();ap.add_argument('step');ap.add_argument('config');ap.add_argument('--sample');ap.add_argument('--outdir');ap.add_argument('--force',action='store_true');ap.add_argument('--retry',action='store_true');ap.add_argument('--resolve-paths',action='store_true');a=ap.parse_args()
    cfg=read_simple_yaml(Path(a.config)); step=a.step
    runtime=resolve_runtime_paths(step,cfg)
    if a.resolve_paths:
        print(f'OUTPUT_DIR={runtime["output_dir"]}');print(f'JOB_ARRAY_DIR={runtime["job_array_dir"]}');print(f'LOG_DIR={runtime["log_dir"]}');return
    # --sample narrows only genuinely sample-level steps. In particular, do
    # not turn a singleton/reference manifest entry into a sample identifier.
    if a.sample and step in SAMPLE_STEPS:
        if step == 'trna_match':
            mapped=table_samples(cfg[step]['paths'].get('sample_reference_map',''))
            candidates=[a.sample] if a.sample in mapped else []
        else: candidates=[a.sample]
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
        elif step in {'codon_match','trna_match'}: source=cfg[step]['paths']['sample_reference_map']
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
    outdir=Path(a.outdir or runtime['job_array_dir']);outdir.mkdir(parents=True,exist_ok=True);stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    purpose='.retry' if a.retry else ''
    task=outdir/f'{step}.{stamp}{purpose}.tasks.txt'; manifest=outdir/f'{step}.{stamp}{purpose}.manifest.tsv'
    tmp=task.with_suffix(task.suffix+'.tmp');tmp.write_text(''.join(r[0]+'\n' for r in rows));os.replace(tmp,task)
    mt=manifest.with_suffix(manifest.suffix+'.tmp')
    with mt.open('w',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(('task_index','item_type','item','expected_input','expected_output','status_at_submission','config','submitted_at'))
        typ='sample' if step in SAMPLE_STEPS else 'reference' if step=='mitos2_annotation' else 'singleton'
        for i,r in enumerate(rows,1):w.writerow((i,typ,*r,a.config,stamp))
    os.replace(mt,manifest)
    pointer=outdir/f'{step}.current.tsv';pt=pointer.with_suffix('.tmp');pt.write_text(str(manifest)+'\n');os.replace(pt,pointer)
    print(f'TASK_FILE={task}');print(f'MANIFEST={manifest}');print(f'COUNT={len(rows)}')
    print(f'OUTPUT_DIR={runtime["output_dir"]}');print(f'JOB_ARRAY_DIR={outdir}');print(f'LOG_DIR={runtime["log_dir"]}')
    print(f'STATS=total candidates {len(candidates)}; already completed {done}; scheduled {len(rows)}; missing inputs {missing}; invalid existing outputs {invalid}',file=sys.stderr)
if __name__=='__main__':main()
