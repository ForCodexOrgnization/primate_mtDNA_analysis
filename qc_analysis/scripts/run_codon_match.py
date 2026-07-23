#!/usr/bin/env python3
"""Annotate coordinate-lifted VCF records with source and human codon matches."""
import argparse,csv,sys
from collections import Counter
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from qc_analysis.lib.match_utils import yaml,rows,info_parse,info_format,source,human_pos,inject_headers,write_summary,sample_names
FIELDS=[('MTCODON_STATUS','Codon match status'),('MTCODON_MATCH','Human codon matches source reference or alternate codon'),('MTCODON_STRICT_PHASE','Strict phase matching enabled'),('MTCODON_GENE_MATCH','Mitochondrial gene match'),('MTCODON_PHASE_MATCH','Codon phase match'),('MTCODON_PRIMATE_GENE','Source gene'),('MTCODON_PRIMATE_CODON','Source reference codon'),('MTCODON_PRIMATE_ALT_CODON','Source alternate codon constructed from SRC_ALT'),('MTCODON_PRIMATE_PHASE','Source codon phase'),('MTCODON_HUMAN_GENE','Human gene'),('MTCODON_HUMAN_CODON','Human codon'),('MTCODON_HUMAN_PHASE','Human codon phase')]
def complement_base(base):
 return {'A':'T','T':'A','C':'G','G':'C'}.get(str(base).upper(),str(base).upper())
def mutate_codon(codon,phase,alt_base):
 if codon in {'',None,'.','NA'} or phase in {'',None,'.','NA'} or alt_base in {'',None,'.','NA'}: return '.'
 try: phase=int(phase)
 except (TypeError,ValueError): return '.'
 bases=list(str(codon).upper()); alt=str(alt_base).upper()
 if len(bases)!=3 or not 1<=phase<=3 or not alt: return '.'
 bases[phase-1]=alt[0]
 return ''.join(bases)
def load(path,key_column=None):
 d={}; ident=None
 for r in rows(path):
  if key_column and ident is None: ident=key_column
  try:k=((r.get(ident,'') if ident else ''),int(r['pos']))
  except (ValueError,KeyError):continue
  d.setdefault(k,r)
 return d,ident
def load_sample_reference_map(path):
 mapping={}
 for row in rows(path):
  sample=(row.get('sample') or '').strip(); reference_key=(row.get('reference_key') or '').strip()
  if sample and reference_key: mapping[sample]=reference_key
 return mapping
def main():
 a=argparse.ArgumentParser();a.add_argument('--config',required=True);a.add_argument('--sample');a.add_argument('--input');a.add_argument('--output');z=a.parse_args(); c=yaml(z.config); sec=c['codon_match'];p=sec['paths'];s=sec['settings']; strict=bool(s.get('strict_phase_match',True))
 reference_table=p.get('reference_codon_table')
 map_table=p.get('sample_reference_map')
 if reference_table and map_table:
  sp,ident=load(reference_table,'reference_key'); sample_references=load_sample_reference_map(map_table)
 else: # compatibility with historical sample-level tables
  sp,ident=load(p['all_primate_position_codon_table'],'sample'); sample_references={}
 hu,_=load(p['human_codon_table'])
 samples=[z.sample] if z.sample else sample_names(c)
 if z.input: samples=[z.sample or Path(z.input).name.split('.')[0]]
 if not samples: raise SystemExit('No samples found; supply --sample or --input.')
 allrows=[]
 for sample in samples:
  inp=Path(z.input) if z.input else Path(p['input_vcf_dir'])/str(s['input_vcf_pattern']).format(sample=sample)
  out=Path(z.output) if z.output else Path(p['output_dir'])/'vcf_codon'/f"{sample}{s['output_suffix']}"
  if not inp.exists(): raise SystemExit(f'Missing input VCF for {sample}: {inp}')
  out.parent.mkdir(parents=True,exist_ok=True); header=[]; body=[]; counts=Counter()
  with inp.open() as f:
   for line in f:
    if line.startswith('#'): header.append(line);continue
    x=line.rstrip('\n').split('\t'); inf=info_parse(x[7]); source_chrom,pos,source_ref,source_alt=source(inf); hp=human_pos(x,inf); reference_key=sample_references.get(sample, sample); sr=sp.get((reference_key,pos)) if pos else None
    hr=hu.get(('',hp)) if hp else None
    strand=(sr or {}).get('strand','+')
    alt_for_codon=complement_base(source_alt) if strand=='-' else source_alt
    vals={'MTCODON_STRICT_PHASE':'yes' if strict else 'no'}
    if not pos or not hp: status='MISSING_COORD'
    elif not sr: status='SKIPPED_NONCODING'
    elif not hr: status='NO_HUMAN_CODON'
    else:
     gene=sr.get('gene',''); phase=str(sr.get('codon_pos_in_triplet','')); ref_codon=sr.get('codon_seq','.'); alt_codon=mutate_codon(ref_codon,phase,alt_for_codon); gm=gene==hr.get('gene',''); pm=phase==str(hr.get('codon_pos_in_triplet','')); match=hr.get('codon_seq','') in {ref_codon,alt_codon}
     vals.update(MTCODON_GENE_MATCH='yes' if gm else 'no',MTCODON_PHASE_MATCH='yes' if pm else 'no',MTCODON_MATCH='yes' if match else 'no')
     status='GENE_MISMATCH' if strict and not gm else 'PHASE_MISMATCH' if strict and not pm else 'PASS' if match else 'MISMATCH'
    vals.update(MTCODON_STATUS=status,MTCODON_PRIMATE_GENE=(sr or {}).get('gene','.'),MTCODON_PRIMATE_CODON=(sr or {}).get('codon_seq','.'),MTCODON_PRIMATE_ALT_CODON=mutate_codon((sr or {}).get('codon_seq'),(sr or {}).get('codon_pos_in_triplet'),alt_for_codon),MTCODON_PRIMATE_PHASE=(sr or {}).get('codon_pos_in_triplet','.'),MTCODON_HUMAN_GENE=(hr or {}).get('gene','.'),MTCODON_HUMAN_CODON=(hr or {}).get('codon_seq','.'),MTCODON_HUMAN_PHASE=(hr or {}).get('codon_pos_in_triplet','.'))
    for k in ('MTCODON_MATCH','MTCODON_GENE_MATCH','MTCODON_PHASE_MATCH'): vals.setdefault(k,'no')
    inf.update(vals);x[7]=info_format(inf);body.append('\t'.join(x)+'\n');counts[status]+=1
  with out.open('w') as f:f.writelines(inject_headers(header,FIELDS,'MTCODON'));f.writelines(body)
  row={'sample':sample,'input_vcf':str(inp),'output_vcf':str(out),'total_records':len(body),**{f'status_{q}':counts[q] for q in ['PASS','SKIPPED_NONCODING','NO_HUMAN_CODON','GENE_MISMATCH','PHASE_MISMATCH','MISMATCH','MISSING_COORD']},'strict_phase_match':strict,'status':'completed'};write_summary(Path(p['reports_dir'])/f'{sample}.codon_match_summary.tsv',row);allrows.append(row)
 if allrows:
  path=Path(p['reports_dir'])/'all_samples.codon_match_summary.tsv';path.parent.mkdir(parents=True,exist_ok=True)
  with path.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(allrows[0]),delimiter='\t');w.writeheader();w.writerows(allrows)
if __name__=='__main__':main()
