#!/usr/bin/env python3
"""Annotate lifted variants with interval and human-reference-guided rRNA data."""
import argparse,csv,sys
from collections import Counter
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from qc_analysis.lib.match_utils import *

BASE=['STATUS','S_GENE','H_GENE','GENE_MATCH','S_LOCAL','H_LOCAL','S_LEN','H_LEN','S_FRAC','H_FRAC','FRAC_DELTA','STRAND_MATCH','REGION_MATCH']
STRUCT=['H_CLASS','H_ELEMENT','H_PAIR_POS','H_PAIR_LOCAL','H_PAIR_TYPE','H_PAIR_STATE','H_ALT_PAIR_TYPE','H_ALT_EFFECT','S_PAIR_EXPECTED_POS','S_PAIR_LIFTED_HPOS','PAIR_POS_MATCH','LOCAL_MATCH','MATCH_TIER']
FIELDS=[('MTRRNA_'+x,'rRNA region or human-reference-guided structural annotation') for x in BASE+STRUCT]

def normalize_rrna_gene(gene):
 return {'12S':'MT-RNR1','RNR1':'MT-RNR1','MT-RNR1':'MT-RNR1','16S':'MT-RNR2','RNR2':'MT-RNR2','MT-RNR2':'MT-RNR2'}.get(str(gene or '').upper(),gene)
norm=normalize_rrna_gene
def load(path,species=False): return rows(path)
def load_rrna_structure_table(path):
 required={'rrna_gene','human_pos','local_pos','struct_class'}; data=rows(path)
 if not data or not required.issubset(data[0]): raise ValueError(f'rRNA structure table {path} is missing required columns: {", ".join(sorted(required))}')
 result={}
 for row in data:
  try: result[(normalize_rrna_gene(row['rrna_gene']),int(row['human_pos']))]=row
  except (ValueError,KeyError): continue
 return result
def hit(rs,pos,chrom='',sample=''):
 for r in rs:
  if sample and r.get('sample',r.get('species','')) not in {'',sample}:continue
  if chrom and r.get('chrom','') not in {'',chrom}:continue
  try:
   if int(r['start'])<=pos<=int(r['end']):return r
  except (ValueError,KeyError):pass
def local(r,point):
 if not r:return '.','.','.'
 n=int(r['end'])-int(r['start'])+1; v=point-int(r['start'])+1 if r.get('strand','+')!='-' else int(r['end'])-point+1
 return v,n,v/n
def infer_species_pair_pos_from_human_pair_local(interval,human_paired_local_pos):
 try:
  v=int(human_paired_local_pos); return str(int(interval['end'])-v+1 if interval.get('strand','+')=='-' else int(interval['start'])+v-1)
 except (TypeError,ValueError,KeyError): return '.'
def map_for(directory,sample):
 for suffix in ('.coordinate_map.tsv','.coordinate_map.tsv.gz'):
  path=Path(directory)/f'{sample}{suffix}'
  if path.exists():return load_coordinate_map(path)
 return {}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--config',required=True);ap.add_argument('--sample');ap.add_argument('--input');ap.add_argument('--output');a=ap.parse_args();c=yaml(a.config);sec=c['rrna_match'];p,s=sec['paths'],sec['settings'];hs=load(p['human_rrna_table']);ss=load(p['species_rrna_table'])
 enabled=bool(s.get('use_rrna_structure_table',False)); spath=s.get('human_rrna_structure_table','');
 if enabled and (not spath or not Path(spath).exists()):raise SystemExit(f'rRNA structure annotation is enabled but human structure table is missing: {spath or "<unset>"}')
 structure=load_rrna_structure_table(spath) if enabled else {}
 samples=[a.sample] if a.sample else sample_names(c)
 if a.input:samples=[a.sample or Path(a.input).name.split('.')[0]]
 allrows=[]
 for sample in samples:
  choices=[(Path(p['input_vcf_dir'])/str(s['input_vcf_pattern']).format(sample=sample),True),(Path(p['fallback_codon_vcf_dir'])/str(s['fallback_codon_vcf_pattern']).format(sample=sample),False),(Path(p['fallback_raw_vcf_dir'])/str(s['fallback_raw_vcf_pattern']).format(sample=sample),False)];inp=Path(a.input) if a.input else next((x for x,_ in choices if x.exists()),None)
  if not inp:raise SystemExit(f'Missing rRNA input VCF for {sample}')
  cmap=map_for(p.get('coordinate_map_dir',''),sample);out=Path(a.output) if a.output else Path(p['output_dir'])/'vcf_rrna'/f"{sample}{s['output_suffix']}";out.parent.mkdir(parents=True,exist_ok=True);head=[];body=[];co=Counter();yes=0
  for line in open_text(inp):
   if line.startswith('#'):head.append(line);continue
   x=line.rstrip().split('\t');inf=info_parse(x[7]);sch,pos,_,_=source(inf);hp=human_pos(x,inf);sr=hit(ss,pos,sch,sample) if pos else None;hr=hit(hs,hp,x[0]) if hp else None
   status='MISSING_COORD' if not pos or not hp else 'NO_SPECIES_OR_HUMAN_RRNA' if not sr and not hr else 'NO_SPECIES_RRNA' if not sr else 'NO_HUMAN_RRNA' if not hr else 'GENE_MISMATCH' if norm(sr.get('rrna_gene'))!=norm(hr.get('rrna_gene')) else 'OK'
   sl,slen,sf=local(sr,pos) if sr else ('.','.','.');hl,hlen,hf=local(hr,hp) if hr else ('.','.','.');gm=bool(sr and hr and norm(sr.get('rrna_gene'))==norm(hr.get('rrna_gene')));strand=bool(sr and hr and sr.get('strand','+')==hr.get('strand','+'));region=gm and (not s.get('require_same_strand',False) or strand);yes+=region
   v={'MTRRNA_STATUS':status,'MTRRNA_S_GENE':norm(sr.get('rrna_gene')) if sr else '.','MTRRNA_H_GENE':norm(hr.get('rrna_gene')) if hr else '.','MTRRNA_GENE_MATCH':'yes' if gm else 'no','MTRRNA_S_LOCAL':sl,'MTRRNA_H_LOCAL':hl,'MTRRNA_S_LEN':slen,'MTRRNA_H_LEN':hlen,'MTRRNA_S_FRAC':sf,'MTRRNA_H_FRAC':hf,'MTRRNA_FRAC_DELTA':abs(sf-hf) if isinstance(sf,float) and isinstance(hf,float) else '.','MTRRNA_STRAND_MATCH':'yes' if strand else 'no','MTRRNA_REGION_MATCH':'yes' if region else 'no'}
   st=structure.get((norm(hr.get('rrna_gene')) if hr else '',hp)) if enabled and hp else None
   structural={f'MTRRNA_{k}':'.' for k in STRUCT};structural['MTRRNA_MATCH_TIER']='NA'
   if st:
    klass=st.get('struct_class','.') or '.'; ppos=st.get('paired_human_pos','.') or '.'; plocal=st.get('paired_local_pos','.') or '.'; refpt=st.get('pair_type','.') or '.'; refpt=pair_type(st.get('base'),st.get('paired_base')) if refpt=='.' else refpt
    expected=infer_species_pair_pos_from_human_pair_local(sr,plocal) if sr else '.'; lifted=lift_source_pos_to_human(expected,cmap); pmatch=compare_values(lifted,ppos); lmatch=compare_values(sl,hl)
    altpt=pair_type(x[4],st.get('paired_base')) if klass=='stem' else '.'
    structural.update({'MTRRNA_H_CLASS':klass,'MTRRNA_H_ELEMENT':st.get('struct_element','.') or '.','MTRRNA_H_PAIR_POS':ppos,'MTRRNA_H_PAIR_LOCAL':plocal,'MTRRNA_H_PAIR_TYPE':refpt,'MTRRNA_H_PAIR_STATE':st.get('pair_state','.') or pair_state(refpt),'MTRRNA_H_ALT_PAIR_TYPE':altpt,'MTRRNA_H_ALT_EFFECT':pair_effect(refpt,altpt) if klass=='stem' else '.','MTRRNA_S_PAIR_EXPECTED_POS':expected,'MTRRNA_S_PAIR_LIFTED_HPOS':lifted,'MTRRNA_PAIR_POS_MATCH':pmatch,'MTRRNA_LOCAL_MATCH':lmatch})
    if status=='OK' and gm and lmatch=='yes' and klass=='stem': structural['MTRRNA_MATCH_TIER']='HIGH_CONF_STEM' if pmatch=='yes' else 'MODERATE_CONF_STEM'
    elif status=='OK' and gm and lmatch=='yes' and klass=='loop': structural['MTRRNA_MATCH_TIER']='HIGH_CONF_LOOP'
    else: structural['MTRRNA_MATCH_TIER']='LOW_CONF'
   v.update(structural);inf.update(v);x[7]=info_format(inf);body.append('\t'.join(x)+'\n');co[status]+=1
  with out.open('w') as f:f.writelines(inject_headers(head,FIELDS,'MTRRNA'));f.writelines(body)
  row={'sample':sample,'input_vcf':str(inp),'output_vcf':str(out),'total_records':len(body),**{f'status_{q}':co[q] for q in ['OK','NO_SPECIES_RRNA','NO_HUMAN_RRNA','NO_SPECIES_OR_HUMAN_RRNA','GENE_MISMATCH','MISSING_COORD']},'rrna_region_match_yes':yes,'rrna_region_match_no':len(body)-yes,'status':'completed'};write_summary(Path(p['reports_dir'])/f'{sample}.rrna_match_summary.tsv',row);allrows.append(row)
 if allrows:
  q=Path(p['reports_dir'])/'all_samples.rrna_match_summary.tsv';q.parent.mkdir(parents=True,exist_ok=True)
  with q.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(allrows[0]),delimiter='\t');w.writeheader();w.writerows(allrows)
if __name__=='__main__':main()
