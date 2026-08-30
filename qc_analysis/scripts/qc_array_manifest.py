#!/usr/bin/env python3
"""Build immutable, race-safe task manifests for QC Slurm arrays."""
import argparse, csv, datetime, os, re, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from qc_analysis.lib.simple_yaml import read_simple_yaml
from qc_analysis.lib.reference_utils import normalized_fasta_sequence_sha256

SAMPLE_STEPS = {"coordinate_liftover", "codon_match", "trna_match", "rrna_match"}
DEFERRED_INPUT_STEPS = set(SAMPLE_STEPS)
WORKFLOW_REVALIDATE_STEPS = set(SAMPLE_STEPS)
GLOBAL_STEPS = {"collect_variant_calling_results", "discover_global_anchor", "build_primate_codon_table",
 "compare_genbank_mitos2", "mitos2_prepare_tasks", "mitos2_merge", "codon_match_validate",
 "codon_match_merge", "build_trna_indexes", "trna_match_merge", "trna_gene_qc", "rrna_match_merge", "build_primate_homo_background",
 "intraspecies_contamination", "sample_variant_filtering", "human_contamination", "final_filter", "interspecies_contamination"}

STEP_SECTIONS = {
    "collect_variant_calling_results": "collect_variant_calling",
    "discover_global_anchor": "global_anchor_discovery",
    "coordinate_liftover": "coordinate_liftover",
    "interspecies_contamination": "interspecies_contamination",
    "human_contamination": "human_contamination",
    "build_primate_homo_background": "primate_homo_background",
    "mitos2_prepare_tasks": "mitos2_annotation", "mitos2_annotation": "mitos2_annotation",
    "mitos2_merge": "mitos2_annotation", "build_primate_codon_table": "build_primate_codon_table",
    "compare_genbank_mitos2": "genbank_mitos2_comparison",
    "codon_match_validate": "codon_match", "codon_match": "codon_match", "codon_match_merge": "codon_match",
    "build_trna_indexes": "trna_match", "trna_match": "trna_match",
    "trna_match_merge": "trna_match", "rrna_match": "rrna_match", "rrna_match_merge": "rrna_match",
    "intraspecies_contamination": "intraspecies_contamination",
    "sample_variant_filtering": "sample_variant_filtering",
    "final_filter": "final_filter",
}
FALLBACK_OUTPUTS = {
    "collect_variant_calling_results": "results/qc/variant_calling_collection",
    "discover_global_anchor": "results/qc/coordinate_liftover/global_anchor",
    "coordinate_liftover": "results/qc/coordinate_liftover",
    "interspecies_contamination": "results/qc/interspecies_contamination",
    "human_contamination": "results/qc/human_contamination",
    "build_primate_homo_background": "results/qc/primate_homo_background",
    "mitos2_prepare_tasks": "results/qc/mitos2_annotation", "mitos2_annotation": "results/qc/mitos2_annotation",
    "mitos2_merge": "results/qc/mitos2_annotation", "build_primate_codon_table": "results/qc/codon_table_build",
    "compare_genbank_mitos2": "results/qc/genbank_mitos2_comparison",
    "codon_match_validate": "results/qc/codon_match", "codon_match": "results/qc/codon_match",
    "codon_match_merge": "results/qc/codon_match", "build_trna_indexes": "results/qc/trna_match",
    "trna_match": "results/qc/trna_match", "trna_match_merge": "results/qc/trna_match",
    "rrna_match": "results/qc/rrna_match", "rrna_match_merge": "results/qc/rrna_match",
    "intraspecies_contamination": "results/qc/intraspecies_contamination",
    "sample_variant_filtering": "results/qc/sample_variant_filtering",
    "final_filter": "results/qc/final_filter",
}

def resolve_runtime_paths(step, cfg):
    section=cfg.get(STEP_SECTIONS.get(step,''),{}) or {}; paths=section.get('paths',{}) or {}
    explicit=paths.get('job_array_dir'); output=paths.get('output_dir') or section.get('outdir')
    if not output:
        reports=paths.get('reports_dir')
        if reports: output=str(Path(reports).parent)
    configured_output=bool(output); output=str(output or FALLBACK_OUTPUTS.get(step, f'.workflow/qc_preprocessing/{step}'))
    fallback_metadata=step not in FALLBACK_OUTPUTS and not configured_output
    job_array=str(explicit or (Path(output) if fallback_metadata else Path(output)/'job_arrays'))
    log_dir=paths.get('log_dir') or str(Path(output)/'logs'/'job_arrays')
    return {'output_dir':output,'job_array_dir':job_array,'log_dir':str(log_dir)}

def table_rows(path):
    path=Path(path)
    if not path.is_file(): return []
    with path.open(newline='') as f: rows=[r for r in csv.reader(f,delimiter='\t') if any(x.strip() for x in r)]
    if not rows:return []
    header=[x.strip() for x in rows[0]]; lower=[x.lower() for x in header]
    if any(x in lower for x in ('sample','sample_id','name','target_species','reference_key')):
        return [dict(zip(header,r)) for r in rows[1:]]
    return []

def table_samples(path):
    path=Path(path)
    if not path.is_file(): return []
    with path.open(newline='') as f: rows=list(csv.reader(f,delimiter='\t'))
    if not rows:return []
    header=[x.strip().lower() for x in rows[0]]; start=1 if any(x in header for x in ('sample','sample_id','name')) else 0
    col=next((header.index(x) for x in ('sample','sample_id','name') if x in header),0)
    return [r[col].strip() for r in rows[start:] if len(r)>col and r[col].strip()]

def sample_species_rows(path):
    """Read sample/species metadata in either headered or legacy two-column form."""
    path=Path(path)
    if not path.is_file(): return []
    with path.open(newline='') as f: rows=[r for r in csv.reader(f,delimiter='\t') if any(x.strip() for x in r)]
    if not rows:return []
    header=[x.strip() for x in rows[0]]; lower=[x.lower() for x in header]
    sample_aliases=('sample','sample_id','name'); species_aliases=('target_species','species')
    sample_col=next((lower.index(x) for x in sample_aliases if x in lower),None)
    species_col=next((lower.index(x) for x in species_aliases if x in lower),None)
    if sample_col is not None and species_col is not None:
        start=1
    else:
        sample_col,species_col,start=0,1,0
    result=[]
    for row in rows[start:]:
        if len(row)<=max(sample_col,species_col): continue
        sample=row[sample_col].strip(); species=row[species_col].strip()
        if sample and species: result.append({'sample':sample,'species':species})
    return result

def species_key(value): return re.sub(r'_+','_',re.sub(r'\s+','_',str(value or '').lower())).strip('_')

def resolved_reference_inventory(cfg):
    """Return current sequence-identity MITOS2 workers and targets with usable FASTAs."""
    sec=cfg.get('mitos2_annotation') or {}; paths=sec.get('paths') or {}; manifest=table_rows(paths.get('reference_manifest',''))
    fasta_dir=Path(paths.get('final_chrM_fasta_dir',paths.get('fasta_dir','references/variant_calling/Ref_chrM')))
    references=set(); resolved_targets=set(); unresolved=[]
    for row in manifest:
        target=(row.get('target_species') or '').strip()
        if not target: continue
        no_chrm=any((row.get(k) or '').strip() in ('wg_only_no_chrM','missing_chrM_ref') for k in ('final_reference_strategy','chrM_reference_context','status'))
        if no_chrm and (row.get('chrM_selection_status') or '').strip()=='missing_chrM_ref': continue
        standardized=fasta_dir/f'{target}.fa'; raw=(row.get('chrM_expected_output_fasta') or '').strip(); manifest_fasta=Path(raw) if raw else None
        fasta=standardized if standardized.is_file() else manifest_fasta if manifest_fasta and manifest_fasta.is_file() else None
        if fasta is None: unresolved.append(target); continue
        try: sequence_sha=normalized_fasta_sequence_sha256(fasta)['sequence_sha256']
        except (OSError,ValueError): unresolved.append(target); continue
        references.add('reference:mtref_'+sequence_sha); resolved_targets.add(species_key(target))
    return sorted(references),resolved_targets,sorted(set(unresolved))

def sample_inventory(cfg):
    """Return the broad current sample inventory without consulting annotation QC."""
    mitos_paths=(cfg.get('mitos2_annotation',{}).get('paths',{}) or {})
    liftover_paths=(cfg.get('coordinate_liftover',{}).get('paths',{}) or {})
    sample_file=mitos_paths.get('sample_ref_file') or liftover_paths.get('sample_ref_file','')
    return table_samples(sample_file) if sample_file else []

def resolved_static_samples(cfg):
    mitos_paths=(cfg.get('mitos2_annotation',{}).get('paths',{}) or {})
    sample_file=mitos_paths.get('sample_ref_file') or cfg.get('coordinate_liftover',{}).get('paths',{}).get('sample_ref_file','')
    rows=sample_species_rows(sample_file); _,resolved_targets,_=resolved_reference_inventory(cfg); samples=[]
    for row in rows:
        sample=(row.get('sample') or '').strip(); target=(row.get('species') or '').strip()
        if sample and species_key(target) in resolved_targets: samples.append(sample)
    return samples

def candidate_samples(step,cfg):
    """Plan from current metadata/reference inventory, never a stale downstream map.

    Downstream maps are preferred only when no current static reference-resolved
    inventory can be derived. As a final safety-first fallback, schedule the broad
    sample inventory and let runtime eligibility perform the current map/QC gate.
    """
    if step=='coordinate_liftover':
        paths=(cfg.get(step,{}).get('paths',{}) or {}); sample_file=paths.get('sample_ref_file','')
        return table_samples(sample_file) if sample_file else []
    if step in {'codon_match','trna_match','rrna_match'}:
        static=resolved_static_samples(cfg)
        if static:return static
    configured=(cfg.get(step,{}).get('paths',{}) or {}).get('sample_reference_map','')
    if configured and Path(configured).is_file(): return table_samples(configured)
    return sample_inventory(cfg)

def mitos2_reference_candidates(cfg):
    """Always derive current MITOS2 work from stable sequence SHA identity."""
    references,_,unresolved=resolved_reference_inventory(cfg)
    if not references:
        detail=', '.join(unresolved[:10])
        raise SystemExit('ERROR: cannot derive MITOS2 reference tasks from current reference_manifest/FASTA files'+(f' (unresolved examples: {detail})' if detail else ''))
    if unresolved: print('WARNING: MITOS2 references could not be pre-resolved for: '+', '.join(unresolved[:20]),file=sys.stderr)
    return references

def valid_vcf(path,tag):
    try:
        if Path(path).stat().st_size==0:return False
        text=Path(path).read_text(errors='ignore'); return '#CHROM' in text and (not tag or tag in text)
    except OSError:return False

def paths_for(step,s,cfg,defer_input=False):
    if step=='coordinate_liftover':
        p=cfg[step]['paths']; return '',str(Path(p['output_dir'])/'vcf_lifted_raw'/f'{s}.lifted.raw.vcf'),'##INFO=<ID=SRC_POS'
    sec=cfg[step];p=sec['paths'];st=sec['settings']
    if step=='codon_match': inp=Path(p['input_vcf_dir'])/st['input_vcf_pattern'].format(sample=s);tag='MTCODON'
    elif step=='trna_match':
        a=Path(p['input_vcf_dir'])/st['input_vcf_pattern'].format(sample=s);b=Path(p['fallback_input_vcf_dir'])/st['fallback_input_vcf_pattern'].format(sample=s);inp=a if defer_input or a.exists() else b;tag='MTTRNA'
    else:
        choices=[Path(p['input_vcf_dir'])/st['input_vcf_pattern'].format(sample=s),Path(p.get('fallback_trna_vcf_dir',p['input_vcf_dir']))/st.get('fallback_trna_vcf_pattern','{sample}.lifted.trna.vcf').format(sample=s),Path(p['fallback_codon_vcf_dir'])/st['fallback_codon_vcf_pattern'].format(sample=s),Path(p['fallback_raw_vcf_dir'])/st['fallback_raw_vcf_pattern'].format(sample=s)]
        inp=choices[0] if defer_input else next((x for x in choices if x.exists()),choices[0]);tag='MTRRNA'
    folder={'codon_match':'vcf_codon','trna_match':'vcf_trna','rrna_match':'vcf_rrna'}[step]
    suffix=st['output_suffix'] if step!='trna_match' or str(inp).startswith(str(p['input_vcf_dir'])) else '.lifted.trna.vcf'
    return str(inp),str(Path(p['output_dir'])/folder/f'{s}{suffix}'),tag

def main():
    ap=argparse.ArgumentParser();ap.add_argument('step');ap.add_argument('config');ap.add_argument('--sample');ap.add_argument('--outdir');ap.add_argument('--force',action='store_true');ap.add_argument('--retry',action='store_true');ap.add_argument('--resolve-paths',action='store_true');ap.add_argument('--workflow-run',action='store_true')
    a=ap.parse_args();cfg=read_simple_yaml(Path(a.config));step=a.step;runtime=resolve_runtime_paths(step,cfg)
    if a.resolve_paths:
        print(f'OUTPUT_DIR={runtime["output_dir"]}');print(f'JOB_ARRAY_DIR={runtime["job_array_dir"]}');print(f'LOG_DIR={runtime["log_dir"]}');return
    if a.sample and step in SAMPLE_STEPS:
        mapped=candidate_samples(step,cfg)
        candidates=[a.sample] if not mapped or a.sample in mapped else []
    elif step in GLOBAL_STEPS:candidates=[step]
    elif step=='mitos2_annotation':candidates=mitos2_reference_candidates(cfg)
    elif step in SAMPLE_STEPS:candidates=candidate_samples(step,cfg)
    else:raise SystemExit(f'Unsupported array step: {step}')
    candidates=sorted(set(x for x in candidates if x))
    if step in SAMPLE_STEPS and not candidates:
        detail=f' for requested sample {a.sample!r}' if a.sample else ''
        raise SystemExit(f'ERROR: no candidate samples for {step}{detail}; no eligible tasks under current sample/reference inventory')
    rows=[];done=missing=invalid=0
    static_inventory=bool(resolved_static_samples(cfg)) if step in {'codon_match','trna_match','rrna_match'} else False
    defer_inputs=(step=='coordinate_liftover' or step=='codon_match' or (step in {'trna_match','rrna_match'} and static_inventory))
    for item in candidates:
        inp=out=tag=''
        if step in SAMPLE_STEPS:inp,out,tag=paths_for(step,item,cfg,defer_input=defer_inputs)
        input_missing=bool(inp and not Path(inp).exists())
        if input_missing:missing+=1
        if input_missing and not defer_inputs:continue
        complete=bool(out and valid_vcf(out,tag))
        if out and Path(out).exists() and not complete:invalid+=1
        if complete:done+=1
        workflow_revalidate=step in WORKFLOW_REVALIDATE_STEPS
        include=workflow_revalidate or ((not complete or a.force) if not a.retry else (not complete))
        if include:
            status='runtime_revalidate' if workflow_revalidate else 'force_rerun' if complete else 'pending_input' if input_missing else 'pending'
            rows.append((item,inp,out,status))
    outdir=Path(a.outdir or runtime['job_array_dir']);outdir.mkdir(parents=True,exist_ok=True);stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ');purpose='.retry' if a.retry else ''
    task=outdir/f'{step}.{stamp}{purpose}.tasks.txt';manifest=outdir/f'{step}.{stamp}{purpose}.manifest.tsv'
    tmp=task.with_suffix(task.suffix+'.tmp');tmp.write_text(''.join(r[0]+'\n' for r in rows));os.replace(tmp,task)
    mt=manifest.with_suffix(manifest.suffix+'.tmp')
    with mt.open('w',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(('task_index','item_type','item','expected_input','expected_output','status_at_submission','config','submitted_at'))
        typ='sample' if step in SAMPLE_STEPS else 'reference' if step=='mitos2_annotation' else 'singleton'
        for i,r in enumerate(rows,1):w.writerow((i,typ,*r,a.config,stamp))
    os.replace(mt,manifest);pointer=outdir/f'{step}.current.tsv';pt=pointer.with_suffix('.tmp');pt.write_text(str(manifest)+'\n');os.replace(pt,pointer)
    print(f'TASK_FILE={task}');print(f'MANIFEST={manifest}');print(f'COUNT={len(rows)}');print(f'OUTPUT_DIR={runtime["output_dir"]}');print(f'JOB_ARRAY_DIR={outdir}');print(f'LOG_DIR={runtime["log_dir"]}')
    print(f'STATE={"complete_noop" if not rows else "scheduled"}')
    print(f'STATS=total candidates {len(candidates)}; already completed {done}; scheduled {len(rows)}; missing inputs {missing}; invalid existing outputs {invalid}',file=sys.stderr)
if __name__=='__main__':main()
