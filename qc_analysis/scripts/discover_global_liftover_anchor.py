#!/usr/bin/env python3
"""Discover reproducible reference-level mtDNA liftover anchors from a multi-reference MSA."""
from __future__ import annotations

import argparse, csv, math, shlex, sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(REPO_ROOT))
from qc_analysis.lib.mt_anchor_utils import AMBIGUOUS_DNA, mask_ambiguity_for_alignment, derive_reference_id, rotate_sequence, rotated_to_original, sequence_sha256
from qc_analysis.lib.alignment_runner import check_aligner_environment, run_aligner
from qc_analysis.scripts.run_coordinate_liftover import read_simple_yaml, read_workflow_config, read_fasta, write_fasta, iter_sample_ref_rows, find_species_fasta, infer_anchor_with_status

@dataclass
class Ref:
    reference_id: str; species: str; family: str; species_fasta: Path; seq: str; sha: str
    samples: List[str]=field(default_factory=list)
    target_sequence: Optional[str]=None; eligible: bool=True; reason: str=""
    coarse_anchor: int=0; coarse_method: str=""; coarse_k: str=""; coarse_status: str="PENDING"


def cfg_section(data, name):
    v = data.get(name, {}) if isinstance(data, dict) else {}
    return v if isinstance(v, dict) else {}

def cfg_get(c,k,default=None): return c.get(k, default)
def cfg_int(c,k,d): return int(cfg_get(c,k,d))
def cfg_float(c,k,d): return float(cfg_get(c,k,d))
def cfg_bool(c,k,d):
    v=cfg_get(c,k,d)
    return v if isinstance(v,bool) else str(v).lower() in {"1","true","yes","on"}

def ambiguity_qc(seq: str) -> dict:
    types = sorted(set(seq) & AMBIGUOUS_DNA)
    count = sum(base in AMBIGUOUS_DNA for base in seq)
    return {
        'ambiguous_base_count': count,
        'ambiguous_base_fraction': count / len(seq) if seq else 0,
        'ambiguous_base_types': ''.join(types),
    }

def read_alignment(path: Path) -> Dict[str,str]:
    recs={}; name=None; chunks=[]
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            if name: recs[name]=''.join(chunks).upper()
            name=line[1:].split()[0]; chunks=[]
        elif line.strip(): chunks.append(line.strip())
    if name: recs[name]=''.join(chunks).upper()
    if len({len(s) for s in recs.values()}) != 1: raise ValueError('MSA records have different lengths')
    return recs

def run_msa(in_fa: Path, out_fa: Path, cfg: dict) -> None:
    aligner=str(cfg_get(cfg,'aligner','mafft')); opts=shlex.split(str(cfg_get(cfg,'aligner_options','--auto --quiet')))
    try:
        run_aligner(aligner, opts, in_fa, out_fa, cfg_int(cfg, 'threads', 1), cfg_bool(cfg, 'use_conda_env', True), str(cfg_get(cfg, 'module_load', 'miniconda/24.11.3')), str(cfg_get(cfg, 'conda_env', 'mafft_env')))
        return
    except RuntimeError:
        if not cfg_bool(cfg,'allow_simple_alignment_fallback',False):
            raise
        recs=[]; n=None; ch=[]
        for l in in_fa.read_text().splitlines():
            if l.startswith('>'):
                if n: recs.append((n,''.join(ch).upper()))
                n=l[1:].split()[0]; ch=[]
            elif l: ch.append(l.strip())
        if n: recs.append((n,''.join(ch).upper()))
        m=max(len(s) for _,s in recs)
        write_multi_fasta(out_fa, [(n, s.ljust(m, '-')) for n, s in recs])
        return

def write_multi_fasta(path: Path, records: Iterable[tuple[str,str]]) -> None:
    with path.open('w') as out:
        for name, seq in records:
            out.write(f'>{name}\n')
            for i in range(0,len(seq),80): out.write(seq[i:i+80]+'\n')

def entropy(counts):
    tot=sum(counts)
    return 0.0 if not tot else -sum((c/tot)*math.log2(c/tot) for c in counts if c)

def homopolymer(seq):
    mx=cur=0; prev=None
    for b in seq.replace('-',''):
        cur = cur+1 if b==prev else 1; prev=b; mx=max(mx,cur)
    return mx

def column_metrics(cols):
    total=len(cols); nongap=sum(1 for b in cols if b!='-'); counts=[sum(1 for b in cols if b==x) for x in 'ACGT']; called=sum(counts)
    return dict(occupancy=nongap/total, major=(max(counts)/called if called else 0), entropy=entropy(counts), gap=(total-nongap)/total)

def score_candidates(aln: Dict[str,str], human_name: str, cfg: dict):
    names=sorted(aln); L=len(next(iter(aln.values()))); w=cfg_int(cfg,'candidate_window_size',31); half=w//2; rows=[]
    for start in range(1, L-w+2):
        end=start+w-1; mets=[]
        for col in range(start-1,end): mets.append(column_metrics([aln[n][col] for n in names]))
        human_gap=sum(1 for col in range(start-1,end) if aln[human_name][col]=='-')
        elig_gap=sum(1 for n in names if n!=human_name for col in range(start-1,end) if aln[n][col]=='-')
        hp=max(homopolymer(aln[n][start-1:end]) for n in names)
        row={
          'alignment_start':start,'alignment_end':end,'center_column':start+half,'number_of_sequences':len(names),
          'mean_occupancy':sum(m['occupancy'] for m in mets)/w,'minimum_column_occupancy':min(m['occupancy'] for m in mets),
          'mean_major_allele_fraction':sum(m['major'] for m in mets)/w,'mean_shannon_entropy':sum(m['entropy'] for m in mets)/w,
          'total_gap_fraction':sum(1 for n in names for b in aln[n][start-1:end] if b=='-')/(len(names)*w),
          'maximum_homopolymer_run':hp,'human_gap_count':human_gap,'eligible_sequence_gap_count':elig_gap,
          'distance_from_alignment_edge':min(start-1,L-end),'candidate_status':'PASS'}
        if row['mean_occupancy'] < cfg_float(cfg,'min_window_mean_occupancy',0.95) or row['total_gap_fraction'] > cfg_float(cfg,'max_window_gap_fraction',0.05) or hp > cfg_int(cfg,'max_homopolymer_run',6) or human_gap:
            row['candidate_status']='FAIL'
        rows.append(row)
    rows.sort(key=lambda r:(r['candidate_status']!='PASS', -r['mean_occupancy'], -r['mean_major_allele_fraction'], r['total_gap_fraction'], r['mean_shannon_entropy'], -r['distance_from_alignment_edge'], r['center_column']))
    return rows

def select_column(aln, human_name, window, cfg):
    names=sorted(aln); best=None
    for c in range(window['alignment_start'], window['alignment_end']+1):
        vals=[aln[n][c-1] for n in names]; m=column_metrics(vals)
        if m['occupancy'] < cfg_float(cfg,'min_anchor_column_occupancy',0.98) or aln[human_name][c-1]=='-': continue
        key=(-m['occupancy'], -m['major'], m['entropy'], c)
        if best is None or key < best[0]: best=(key,c,m)
    if not best: raise RuntimeError('No eligible global anchor column')
    return best[1], best[2]

def write_exclusions(out: Path, manifest_fields: List[str], refs: List[Ref], unresolved_rows: List[dict]) -> None:
    with (out/'excluded_references.tsv').open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=manifest_fields,delimiter='\t',extrasaction='ignore')
        w.writeheader()
        for row in unresolved_rows:
            w.writerow({
                'reference_id':row.get('reference_id',''),
                'species':row.get('species') or row.get('sample',''),
                'species_fasta':row.get('species_fasta',''),
                'sequence_sha256':'',
                'sequence_length':'',
                'n_fraction':'', 'ambiguous_base_count':'', 'ambiguous_base_fraction':'', 'ambiguous_base_types':'',
                'sample_count':1 if row.get('sample') else '',
                'sample_names':row.get('sample',''),
                'discovery_eligible':False,
                'exclusion_reason':row.get('exclusion_reason',''),
                'coarse_anchor_position':'',
                'coarse_anchor_method':'',
                'coarse_anchor_kmer_length':'',
                'coarse_anchor_status':'NOT_RUN',
            })
        for r in refs:
            if not r.eligible:
                w.writerow({'reference_id':r.reference_id,'species':r.species,'species_fasta':r.species_fasta,'sequence_sha256':r.sha,'sequence_length':len(r.seq),'n_fraction':r.seq.count('N')/len(r.seq) if r.seq else 1,**ambiguity_qc(r.seq),'sample_count':len(r.samples),'sample_names':','.join(sorted(r.samples)),'discovery_eligible':r.eligible,'exclusion_reason':r.reason,'coarse_anchor_position':r.coarse_anchor,'coarse_anchor_method':r.coarse_method,'coarse_anchor_kmer_length':r.coarse_k,'coarse_anchor_status':r.coarse_status})

def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument('--config',required=True); ap.add_argument('--check-environment', action='store_true'); args=ap.parse_args(argv)
    config_path=Path(args.config)
    data=read_simple_yaml(config_path); cfg=cfg_section(data,'global_anchor_discovery')
    if args.check_environment:
        info = check_aligner_environment(str(cfg_get(cfg, 'aligner', 'mafft')), shlex.split(str(cfg_get(cfg, 'aligner_options', '--auto --quiet'))), cfg_int(cfg, 'threads', 1), cfg_bool(cfg, 'use_conda_env', True), str(cfg_get(cfg, 'module_load', 'miniconda/24.11.3')), str(cfg_get(cfg, 'conda_env', 'mafft_env')))
        for key in ('aligner', 'resolved_executable', 'version', 'threads', 'environment', 'status'):
            print(f'{key}={info[key]}')
        return 0
    liftover_cfg=read_workflow_config(config_path)
    out=Path(cfg_get(cfg,'output_dir','results/qc/global_anchor')); out.mkdir(parents=True,exist_ok=True)
    human=read_fasta(Path(cfg_get(cfg,'human_fasta','data/reference_tables/human_chrM.fa')))
    human_sha=sequence_sha256(human.seq)
    minlen,maxlen=cfg_int(cfg,'min_reference_length',14000),cfg_int(cfg,'max_reference_length',19000); maxn=cfg_float(cfg,'max_n_fraction',0.02)
    bysha: Dict[str,Ref]={}; exclusions=[]; manifest_rows=0; resolved_sample_rows=0
    for row in iter_sample_ref_rows(Path(cfg_get(cfg,'sample_ref_file','config/sample_ref_file.tsv'))):
        manifest_rows += 1
        sp=row.get('species') or row.get('sample') or ''
        try:
            fasta=Path(row.get('species_fasta')) if row.get('species_fasta') else find_species_fasta(sp, liftover_cfg)
            rec=read_fasta(fasta,row.get('target_sequence') or None); seq=rec.seq; sha=sequence_sha256(seq); rid=(row.get('reference_id') or derive_reference_id(sp,fasta,sha))
            resolved_sample_rows += 1
        except Exception as e: exclusions.append({**row,'species':sp,'species_fasta':row.get('species_fasta',''),'exclusion_reason':str(e)}); continue
        if sha not in bysha: bysha[sha]=Ref(rid,sp,row.get('family') or '',fasta,seq,sha,[row.get('sample','')],row.get('target_sequence') or None)
        else: bysha[sha].samples.append(row.get('sample',''))
    refs=sorted(bysha.values(), key=lambda r:(r.reference_id,r.sha))
    dup_ids=defaultdict(set)
    for r in refs: dup_ids[r.reference_id].add(r.sha)
    for r in refs:
        nfrac=r.seq.count('N')/len(r.seq) if r.seq else 1
        reasons=[]
        if not r.seq: reasons.append('EMPTY_SEQUENCE')
        if len(r.seq)<minlen: reasons.append('REFERENCE_TOO_SHORT')
        if len(r.seq)>maxlen: reasons.append('REFERENCE_TOO_LONG')
        if nfrac>maxn: reasons.append('HIGH_N_FRACTION')
        if len(dup_ids[r.reference_id])>1: reasons.append('ANCHOR_REFERENCE_ID_COLLISION')
        r.reason=';'.join(reasons); r.eligible=not reasons
    manifest_fields='reference_id species species_fasta sequence_sha256 sequence_length n_fraction ambiguous_base_count ambiguous_base_fraction ambiguous_base_types sample_count sample_names discovery_eligible exclusion_reason coarse_anchor_position coarse_anchor_method coarse_anchor_kmer_length coarse_anchor_status'.split()
    unresolved_sample_rows=len(exclusions)
    if not refs:
        with (out/'unique_reference_manifest.tsv').open('w',newline='') as h:
            csv.DictWriter(h,fieldnames=manifest_fields,delimiter='\t').writeheader()
        write_exclusions(out, manifest_fields, refs, exclusions)
        raise RuntimeError("No species references were resolved from sample_ref_file")
    human_anchor=1
    rotated=[]
    for r in refs:
        if r.eligible:
            if any(s.strip() for s in r.samples): pass
            r.coarse_anchor, fb = infer_anchor_with_status(r.seq,human.seq); r.coarse_method='PAIRWISE_SHARED_KMER'; r.coarse_k='31,25,21,15'; r.coarse_status='FAILED' if fb else 'PASS'
            if fb: r.eligible=False; r.reason=(r.reason+';' if r.reason else '')+'COARSE_ANCHOR_FAILED'
            else: rotated.append((r.reference_id, rotate_sequence(mask_ambiguity_for_alignment(r.seq),r.coarse_anchor)))
    with (out/'unique_reference_manifest.tsv').open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=manifest_fields,delimiter='\t'); w.writeheader()
        for r in refs: w.writerow({'reference_id':r.reference_id,'species':r.species,'species_fasta':r.species_fasta,'sequence_sha256':r.sha,'sequence_length':len(r.seq),'n_fraction':r.seq.count('N')/len(r.seq) if r.seq else 1,**ambiguity_qc(r.seq),'sample_count':len(r.samples),'sample_names':','.join(sorted(r.samples)),'discovery_eligible':r.eligible,'exclusion_reason':r.reason,'coarse_anchor_position':r.coarse_anchor,'coarse_anchor_method':r.coarse_method,'coarse_anchor_kmer_length':r.coarse_k,'coarse_anchor_status':r.coarse_status})
    write_exclusions(out, manifest_fields, refs, exclusions)
    if not rotated:
        raise RuntimeError("No eligible species references remained after QC")
    infa=out/'unique_references.coarse_rotated.fa'; write_multi_fasta(infa, [('human_chrM',rotate_sequence(mask_ambiguity_for_alignment(human.seq),human_anchor)), *rotated])
    msa_records=[('human_chrM',rotate_sequence(mask_ambiguity_for_alignment(human.seq),human_anchor)), *rotated]
    if len(msa_records) < 2:
        raise RuntimeError("Global anchor MSA requires human plus at least one species reference")
    alnfa=out/'all_references.aligned.fa'; run_msa(infa,alnfa,cfg); aln=read_alignment(alnfa)
    candidates=score_candidates(aln,'human_chrM',cfg); selected=next((r for r in candidates if r['candidate_status']=='PASS'), None)
    if not selected: raise RuntimeError('No PASS candidate windows')
    col, cm=select_column(aln,'human_chrM',selected,cfg); selected['selected_anchor_column']=col; selected['anchor_column_occupancy']=cm['occupancy']; selected['anchor_column_major_allele_fraction']=cm['major']
    cfields=list(candidates[0].keys())+['selected_anchor_column','anchor_column_occupancy','anchor_column_major_allele_fraction']
    with (out/'global_anchor_candidates.tsv').open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=cfields,delimiter='\t',extrasaction='ignore'); w.writeheader(); w.writerows(candidates)
    with (out/'global_anchor_selection.tsv').open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=cfields+['selection_criteria'],delimiter='\t'); w.writeheader(); w.writerow({**selected,'selection_criteria':'deterministic occupancy/conservation/gap/entropy/edge/column ranking'})
    pos_fields='reference_id species species_fasta sequence_sha256 sequence_length ambiguous_base_count ambiguous_base_fraction ambiguous_base_types anchor_original_position anchor_alignment_column anchor_base human_anchor_original_position anchor_method anchor_window_start anchor_window_end anchor_column_occupancy anchor_window_mean_occupancy anchor_window_gap_fraction anchor_window_conservation anchor_qc_status anchor_qc_notes'.split()
    rows=[]; fails=[]; human_rot=sum(1 for b in aln['human_chrM'][:col] if b!='-'); human_orig=rotated_to_original(human_rot,human_anchor,len(human.seq))
    for r in refs:
        if not r.eligible or r.reference_id not in aln:
            fails.append({'reference_id':r.reference_id,'anchor_qc_status':'DISCOVERY_EXCLUDED','anchor_qc_notes':r.reason}); continue
        base=aln[r.reference_id][col-1]
        if base=='-': fails.append({'reference_id':r.reference_id,'anchor_qc_status':'GLOBAL_ANCHOR_UNAVAILABLE','anchor_qc_notes':'gap_at_selected_column'}); continue
        rot=sum(1 for b in aln[r.reference_id][:col] if b!='-'); orig=rotated_to_original(rot,r.coarse_anchor,len(r.seq)); original_base=r.seq[orig-1]
        if original_base in AMBIGUOUS_DNA: fails.append({'reference_id':r.reference_id,'anchor_qc_status':'GLOBAL_ANCHOR_AMBIGUOUS_BASE','anchor_qc_notes':f'ambiguous_base={original_base}'}); continue
        rows.append({'reference_id':r.reference_id,'species':r.species,'species_fasta':r.species_fasta,'sequence_sha256':r.sha,'sequence_length':len(r.seq),**ambiguity_qc(r.seq),'anchor_original_position':orig,'anchor_alignment_column':col,'anchor_base':original_base,'human_anchor_original_position':human_orig,'anchor_method':'GLOBAL_MSA_ANCHOR','anchor_window_start':selected['alignment_start'],'anchor_window_end':selected['alignment_end'],'anchor_column_occupancy':cm['occupancy'],'anchor_window_mean_occupancy':selected['mean_occupancy'],'anchor_window_gap_fraction':selected['total_gap_fraction'],'anchor_window_conservation':selected['mean_major_allele_fraction'],'anchor_qc_status':'PASS','anchor_qc_notes':''})
    for fn in ['reference_anchor_positions.tsv','family_anchor_positions.tsv']:
        with (out/fn).open('w',newline='') as h: w=csv.DictWriter(h,fieldnames=pos_fields,delimiter='\t'); w.writeheader(); w.writerows(rows)
    with (out/'reference_anchor_failures.tsv').open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=['reference_id','anchor_qc_status','anchor_qc_notes'],delimiter='\t'); w.writeheader(); w.writerows(fails)
    (out/'family_anchor_candidates.tsv').write_text('family\tstatus\tnote\nALL\tNOT_RUN\tglobal_anchor_available_or_no_family_fallback_needed\n')
    summary={
        'manifest_rows':manifest_rows,
        'resolved_sample_rows':resolved_sample_rows,
        'unresolved_sample_rows':unresolved_sample_rows,
        'unique_references':len(refs),
        'eligible_references':len(rotated),
        'coarse_anchor_failed_references':sum(1 for r in refs if r.coarse_status == 'FAILED'),
        'msa_sequence_count':len(msa_records),
        'reference_anchors_written':len(rows),
        'reference_anchor_failures':len(fails),
        'selected_anchor_column':col,
    }
    (out/'global_anchor_summary.tsv').write_text(''.join(f'{k}\t{v}\n' for k,v in summary.items()))
    if not rows:
        raise RuntimeError("Global anchor was selected, but it could not be projected onto any species reference")
    return 0
if __name__=='__main__': sys.exit(main())
