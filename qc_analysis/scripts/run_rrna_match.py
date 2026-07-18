#!/usr/bin/env python3
"""Interval-based rRNA annotation of coordinate-lifted variants; never filters records."""
import argparse,csv,sys
from collections import Counter
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from qc_analysis.lib.match_utils import *
FIELDS=[('MTRRNA_'+x,'rRNA region match annotation') for x in ['STATUS','S_GENE','H_GENE','GENE_MATCH','S_LOCAL','H_LOCAL','S_LEN','H_LEN','S_FRAC','H_FRAC','FRAC_DELTA','STRAND_MATCH','REGION_MATCH']]
def norm(x): return {'12S':'MT-RNR1','RNR1':'MT-RNR1','MT-RNR1':'MT-RNR1','16S':'MT-RNR2','RNR2':'MT-RNR2','MT-RNR2':'MT-RNR2'}.get((x or '').upper(),x)
def load(path,species=False):
 return rows(path)
def hit(rs,pos,chrom='',sample=''):
 for r in rs:
  if sample and r.get('sample',r.get('species','')) not in {'',sample}:continue
  if chrom and r.get('chrom','') not in {'',chrom}:continue
  try:
   if int(r['start'])<=pos<=int(r['end']):return r
  except (ValueError,KeyError):pass
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--config',required=True);ap.add_argument('--sample');ap.add_argument('--input');ap.add_argument('--output');a=ap.parse_args();c=yaml(a.config);sec=c['rrna_match'];p,s=sec['paths'],sec['settings']; hs=load(p['human_rrna_table']);ss=load(p['species_rrna_table']); samples=[a.sample] if a.sample else sample_names(c)
 if a.input:samples=[a.sample or Path(a.input).name.split('.')[0]]
 allrows=[]
 for sample in samples:
  choices=[(Path(p['input_vcf_dir'])/str(s['input_vcf_pattern']).format(sample=sample),True),(Path(p['fallback_codon_vcf_dir'])/str(s['fallback_codon_vcf_pattern']).format(sample=sample),False),(Path(p['fallback_raw_vcf_dir'])/str(s['fallback_raw_vcf_pattern']).format(sample=sample),False)]; inp=Path(a.input) if a.input else next((x for x,_ in choices if x.exists()),None)
  if not inp:raise SystemExit(f'Missing rRNA input VCF for {sample}')
  out=Path(a.output) if a.output else Path(p['output_dir'])/'vcf_rrna'/f"{sample}{s['output_suffix']}";out.parent.mkdir(parents=True,exist_ok=True);head=[];body=[];co=Counter();yes=0
  for line in inp.open():
   if line.startswith('#'):head.append(line);continue
   x=line.rstrip().split('\t');inf=info_parse(x[7]);sch,pos,_,_=source(inf);hp=human_pos(x,inf);sr=hit(ss,pos,sch,sample) if pos else None;hr=hit(hs,hp,x[0]) if hp else None
   if not pos or not hp: status='MISSING_COORD'
   elif not sr and not hr:status='NO_SPECIES_OR_HUMAN_RRNA'
   elif not sr:status='NO_SPECIES_RRNA'
   elif not hr:status='NO_HUMAN_RRNA'
   elif norm(sr.get('rrna_gene'))!=norm(hr.get('rrna_gene')):status='GENE_MISMATCH'
   else:status='OK'
   def vals(r,point):
    if not r:return '.','.','.'
    n=int(r['end'])-int(r['start'])+1; local=point-int(r['start'])+1 if r.get('strand','+')!='-' else int(r['end'])-point+1;return local,n,local/n
   sl,slen,sf=vals(sr,pos) if sr else ('.','.','.');hl,hlen,hf=vals(hr,hp) if hr else ('.','.','.');gm=bool(sr and hr and norm(sr.get('rrna_gene'))==norm(hr.get('rrna_gene')));strand=bool(sr and hr and sr.get('strand','+')==hr.get('strand','+')); region=gm and (not s.get('require_same_strand',False) or strand);yes+=region
   inf.update({'MTRRNA_STATUS':status,'MTRRNA_S_GENE':norm(sr.get('rrna_gene')) if sr else '.','MTRRNA_H_GENE':norm(hr.get('rrna_gene')) if hr else '.','MTRRNA_GENE_MATCH':'yes' if gm else 'no','MTRRNA_S_LOCAL':sl,'MTRRNA_H_LOCAL':hl,'MTRRNA_S_LEN':slen,'MTRRNA_H_LEN':hlen,'MTRRNA_S_FRAC':sf,'MTRRNA_H_FRAC':hf,'MTRRNA_FRAC_DELTA':abs(sf-hf) if isinstance(sf,float) and isinstance(hf,float) else '.','MTRRNA_STRAND_MATCH':'yes' if strand else 'no','MTRRNA_REGION_MATCH':'yes' if region else 'no'});x[7]=info_format(inf);body.append('\t'.join(x)+'\n');co[status]+=1
  with out.open('w') as f:f.writelines(inject_headers(head,FIELDS,'MTRRNA'));f.writelines(body)
  row={'sample':sample,'input_vcf':str(inp),'output_vcf':str(out),'total_records':len(body),**{f'status_{q}':co[q] for q in ['OK','NO_SPECIES_RRNA','NO_HUMAN_RRNA','NO_SPECIES_OR_HUMAN_RRNA','GENE_MISMATCH','MISSING_COORD']},'rrna_region_match_yes':yes,'rrna_region_match_no':len(body)-yes,'status':'completed'};write_summary(Path(p['reports_dir'])/f'{sample}.rrna_match_summary.tsv',row);allrows.append(row)
 if allrows:
  q=Path(p['reports_dir'])/'all_samples.rrna_match_summary.tsv';q.parent.mkdir(parents=True,exist_ok=True)
  with q.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(allrows[0]),delimiter='\t');w.writeheader();w.writerows(allrows)
if __name__=='__main__':main()
