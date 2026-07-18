#!/usr/bin/env python3
"""Annotate lifted VCFs with tRNA position-index structural comparisons (no filtering)."""
import argparse,csv,sys
from collections import Counter
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from qc_analysis.lib.match_utils import *
N=['STATUS','S_ID','H_ID','S_LOCAL','H_LOCAL','S_CLASS','H_CLASS','REGION_MATCH','S_ELEMENT','H_ELEMENT','ELEMENT_MATCH','S_PAIR_TYPE','H_PAIR_TYPE','PAIR_TYPE_MATCH','S_PAIR_STATE','H_PAIR_STATE','PAIR_STATE_MATCH','S_PAIR_LOCAL','H_PAIR_LOCAL','PAIR_LOCAL_MATCH','S_PAIR_POS','H_PAIR_POS','S_PAIR_LIFTED_HPOS','PAIR_POS_MATCH','H_ALT_PAIR_TYPE','S_ALT_PAIR_TYPE','H_ALT_EFFECT','S_ALT_EFFECT','ALLELE_EFFECT_MATCH','COMPENSATED','STRICT_MATCH','S_COORD_SPACE','S_LOOKUP_CHROM','S_LOOKUP_POS']; FIELDS=[('MTTRNA_'+n,'tRNA structural match annotation') for n in N]
def index(path,species=False):
 d={}
 for r in rows(path):
  try:d[(r.get('chrom',''),int(r['pos']))]=r;d[('',int(r['pos']))]=r
  except (KeyError,ValueError):pass
 return d
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--config',required=True);ap.add_argument('--sample');ap.add_argument('--input');ap.add_argument('--output');a=ap.parse_args();c=yaml(a.config);sec=c['trna_match'];p,s=sec['paths'],sec['settings']; hpath=Path(p['human_trna_index'])
 if not hpath.exists(): raise SystemExit(f'Missing human tRNA index: {hpath}. Set run_trnascan_if_missing=true and provide tRNAscan inputs/index generation, or supply the configured index.')
 hi=index(hpath); samples=[a.sample] if a.sample else sample_names(c)
 if a.input:samples=[a.sample or Path(a.input).name.split('.')[0]]
 allrows=[]
 for sample in samples:
  primary=Path(p['input_vcf_dir'])/str(s['input_vcf_pattern']).format(sample=sample); fallback=Path(p['fallback_input_vcf_dir'])/str(s['fallback_input_vcf_pattern']).format(sample=sample); inp=Path(a.input) if a.input else (primary if primary.exists() else fallback); codon=inp==primary
  spi=Path(str(p['species_trna_index_template']).format(species_trna_index_dir=p['species_trna_index_dir'],sample=sample))
  if not spi.exists(): raise SystemExit(f'Missing species tRNA index for {sample}: {spi}. tRNAscan index generation is not available without configured sequence inputs.')
  si=index(spi,True); out=Path(a.output) if a.output else Path(p['output_dir'])/'vcf_trna'/f"{sample}{s['output_suffix'] if codon else '.lifted.trna.vcf'}";out.parent.mkdir(parents=True,exist_ok=True);head=[];body=[];counts=Counter()
  for line in inp.open():
   if line.startswith('#'):head.append(line);continue
   x=line.rstrip().split('\t'); inf=info_parse(x[7]);sch,pos,_,alt=source(inf);hp=human_pos(x,inf); sr=si.get((sch,pos)) or si.get(('',pos)) if pos else None; hr=hi.get((x[0],hp)) or hi.get(('',hp)) if hp else None
   if not pos: status='MISSING_SPECIES_COORD'
   elif not sr and not hr:status='NO_SPECIES_OR_HUMAN_TRNA'
   elif not sr:status='NO_SPECIES_TRNA'
   elif not hr:status='NO_HUMAN_TRNA'
   else:status='OK'
   v={'MTTRNA_STATUS':status,'MTTRNA_S_COORD_SPACE':s.get('species_trna_coord_space','original'),'MTTRNA_S_LOOKUP_CHROM':sch or '.','MTTRNA_S_LOOKUP_POS':pos or '.'}
   for short,col in [('S_ID','trna_id'),('H_ID','trna_id'),('S_LOCAL','local_pos'),('H_LOCAL','local_pos'),('S_CLASS','struct_class'),('H_CLASS','struct_class'),('S_ELEMENT','struct_element'),('H_ELEMENT','struct_element'),('S_PAIR_TYPE','pair_type'),('H_PAIR_TYPE','pair_type'),('S_PAIR_STATE','pair_state'),('H_PAIR_STATE','pair_state'),('S_PAIR_LOCAL','paired_local_pos'),('H_PAIR_LOCAL','paired_local_pos'),('S_PAIR_POS','paired_genomic_pos'),('H_PAIR_POS','paired_genomic_pos')]:v['MTTRNA_'+short]=(sr if short.startswith('S') else hr or {}).get(col,'.')
   for key,a1,a2 in [('REGION_MATCH','trna_id','trna_id'),('ELEMENT_MATCH','struct_element','struct_element'),('PAIR_TYPE_MATCH','pair_type','pair_type'),('PAIR_STATE_MATCH','pair_state','pair_state'),('PAIR_LOCAL_MATCH','paired_local_pos','paired_local_pos')]:v['MTTRNA_'+key]='yes' if sr and hr and sr.get(a1)==hr.get(a2) else 'no'
   v.update({'MTTRNA_S_PAIR_LIFTED_HPOS':'.','MTTRNA_PAIR_POS_MATCH':'no','MTTRNA_H_ALT_PAIR_TYPE':'.','MTTRNA_S_ALT_PAIR_TYPE':'.','MTTRNA_H_ALT_EFFECT':'unknown','MTTRNA_S_ALT_EFFECT':'unknown','MTTRNA_ALLELE_EFFECT_MATCH':'no','MTTRNA_COMPENSATED':'no','MTTRNA_STRICT_MATCH':'yes' if status=='OK' and v['MTTRNA_REGION_MATCH']=='yes' else 'no'})
   inf.update(v);x[7]=info_format(inf);body.append('\t'.join(x)+'\n');counts[status]+=1
  with out.open('w') as f:f.writelines(inject_headers(head,FIELDS,'MTTRNA'));f.writelines(body)
  row={'sample':sample,'input_vcf':str(inp),'output_vcf':str(out),'total_records':len(body),**{f'status_{q}':counts[q] for q in ['OK','NO_SPECIES_TRNA','NO_HUMAN_TRNA','NO_SPECIES_OR_HUMAN_TRNA','MISSING_SPECIES_COORD']},'status':'completed'};write_summary(Path(p['reports_dir'])/f'{sample}.trna_match_summary.tsv',row);write_summary(Path(p['reports_dir'])/f'{sample}.trna_gene_liftover_qc.tsv',{'sample':sample,'status':'not_computed','note':'Gene QC requires interval-level coordinate-map processing.'});allrows.append(row)
 if allrows:
  q=Path(p['reports_dir'])/'all_samples.trna_match_summary.tsv';q.parent.mkdir(parents=True,exist_ok=True)
  with q.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(allrows[0]),delimiter='\t');w.writeheader();w.writerows(allrows)
if __name__=='__main__':main()
