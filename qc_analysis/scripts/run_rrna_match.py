#!/usr/bin/env python3
"""Annotate lifted variants with rRNA interval and two-sided structure data."""
import argparse,csv,sys
from collections import Counter
from pathlib import Path
import re

sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from qc_analysis.lib.match_utils import *

BASE=['STATUS','S_GENE','H_GENE','GENE_MATCH','S_LOCAL','H_LOCAL','S_LEN','H_LEN','S_FRAC','H_FRAC','FRAC_DELTA','STRAND_MATCH','REGION_MATCH']
STRUCT=[
 'H_CLASS','S_CLASS','H_ELEMENT','S_ELEMENT','H_PAIR_POS','S_PAIR_POS',
 'H_PAIR_LOCAL','S_PAIR_LOCAL','H_PAIR_TYPE','S_PAIR_TYPE','H_PAIR_STATE',
 'S_PAIR_STATE','H_ALT_PAIR_TYPE','H_ALT_EFFECT','S_PAIR_LIFTED_HPOS',
 'STRUCTURE_MATCH','PAIR_RELATION_MATCH','LOCAL_MATCH','MATCH_TIER',
 # Deprecated compatibility fields.  They are no longer human-projected.
 'S_PAIR_EXPECTED_POS','PAIR_POS_MATCH',
]

def field_specs():
 specs=[]
 for name in BASE:
  specs.append((f'MTRRNA_{name}','rRNA interval annotation'))
 for name in STRUCT:
  desc='independent human/reference rRNA secondary-structure annotation'
  if name=='S_PAIR_EXPECTED_POS':
   desc='Deprecated; human-projected species pair coordinate is no longer populated'
  elif name=='PAIR_POS_MATCH':
   desc='Deprecated alias of MTRRNA_PAIR_RELATION_MATCH'
  specs.append((f'MTRRNA_{name}',desc))
 return specs

FIELDS=field_specs()
MISSING={'',None,'.','NA','None'}

def normalize_rrna_gene(gene):
 key=re.sub(r'[^A-Z0-9]','',str(gene or '').upper())
 if key.startswith('MT'): key=key[2:]
 return {'12S':'MT-RNR1','RNR1':'MT-RNR1','RRNS':'MT-RNR1','16S':'MT-RNR2','RNR2':'MT-RNR2','RRNL':'MT-RNR2'}.get(key,gene)
norm=normalize_rrna_gene
def load(path,species=False): return rows(path)

def load_species_rrna_regions(path):
 data=rows(path)
 if not data: raise ValueError(f'MITOS2 reference rRNA region table {path} is empty')
 required={'reference_key','rrna_gene','start','end','strand'}
 missing=required-set(data[0])
 if missing: raise ValueError(f'MITOS2 reference rRNA region table {path} is missing required columns: {", ".join(sorted(missing))}')
 return data

def first_value(row,names,default='.'):
 for name in names:
  value=row.get(name)
  if value not in MISSING:return value
 return default

def int_or_none(value):
 try:return int(value)
 except (TypeError,ValueError):return None

def has_pos(value): return int_or_none(value) is not None

def normalize_struct_class(row,pair_pos):
 klass=str(row.get('struct_class','') or '').strip().lower()
 if has_pos(pair_pos): return 'stem'
 if klass in {'stem','paired'}: return 'unknown'
 if klass in {'loop','unpaired'}: return 'loop'
 return 'unknown'

def normalize_structure_row(row,coord_field,pair_pos_fields):
 pair_pos=first_value(row,pair_pos_fields)
 klass=normalize_struct_class(row,pair_pos)
 pair_kind=rrna_pair_type(row.get('base'),row.get('paired_base')) if klass=='stem' else '.'
 pair_state_value=rrna_pair_state(pair_kind,klass)
 return {
  'rrna_gene':norm(row.get('rrna_gene')),
  'pos':str(row.get(coord_field,'.') or '.'),
  'local_pos':first_value(row,('local_pos','model_position')),
  'base':row.get('base','.') or '.',
  'struct_class':klass,
  'struct_element':row.get('struct_element','.') or '.',
  'paired_pos':str(pair_pos) if pair_pos not in MISSING else '.',
  'paired_local_pos':first_value(row,('paired_local_pos','pair_local_pos')),
  'paired_base':row.get('paired_base','.') or '.',
  'pair_type':pair_kind,
  'pair_state':pair_state_value,
  'structure_source':row.get('structure_source','.') or '.',
  'coordinate_reference_sequence_sha256':row.get('coordinate_reference_sequence_sha256','') or '',
 }

def load_structure_index(path,key_mode):
 data=rows(path)
 if not data: raise ValueError(f'rRNA structure table {path} is empty')
 headers=set(data[0])
 coord_options=('human_pos','genomic_pos','pos') if key_mode=='human' else ('genomic_pos','pos','human_pos')
 coord_field=next((name for name in coord_options if name in headers),None)
 required={'rrna_gene','struct_class'}
 if key_mode=='species': required.add('reference_key')
 if not coord_field: required.add('genomic_pos' if key_mode=='species' else 'human_pos')
 missing=required-headers
 if missing: raise ValueError(f'rRNA structure table {path} is missing required columns: {", ".join(sorted(missing))}')
 result={};reference_sha={}
 pair_fields=('paired_human_pos','paired_genomic_pos','paired_pos') if key_mode=='human' else ('paired_genomic_pos','paired_human_pos','paired_pos')
 for row in data:
  pos=int_or_none(row.get(coord_field))
  if pos is None: continue
  normalized=normalize_structure_row(row,coord_field,pair_fields)
  gene=norm(row.get('rrna_gene'))
  if key_mode=='human':
   result[(gene,pos)]=normalized
  else:
   ref=row.get('reference_key','')
   if not ref: continue
   result[(ref,gene,pos)]=normalized
   sha=row.get('coordinate_reference_sequence_sha256','') or ''
   if sha:
    previous=reference_sha.get(ref)
    reference_sha[ref]='CONFLICT' if previous and previous!=sha else sha
 return result,reference_sha

def load_rrna_structure_table(path):
 """Backward-compatible human-table loader used by tests and older callers."""
 result,_sha=load_structure_index(path,'human')
 return result

def load_species_rrna_structure_table(path):
 return load_structure_index(path,'species')

def load_sample_reference_map(path):
 if not path or not Path(path).exists():return {}
 result={}
 for row in rows(path):
  sample=row.get('sample','')
  if sample and sample not in result: result[sample]=row
 return result

def hit(rs,pos,chrom='',sample=''):
 for r in rs:
  sample_columns=('sample','species','species_key','accession','reference_id')
  identifiers=[r.get(column,'') for column in sample_columns if column in r and r.get(column,'')]
  if sample and identifiers and sample not in identifiers:continue
  if chrom and r.get('chrom','') not in {'',chrom}:continue
  try:
   if int(r['start'])<=pos<=int(r['end']):return r
  except (ValueError,KeyError):pass

def species_region_for(sample_ref,regions,pos,chrom=''):
 """Find a source-coordinate rRNA only within the sample's exact reference."""
 if not sample_ref or not pos:return None
 reference_key=sample_ref.get('reference_key','')
 if not reference_key:return None
 expected_sha=sample_ref.get('coordinate_reference_sequence_sha256','') or ''
 for region in regions:
  if region.get('reference_key','') != reference_key:continue
  observed_sha=region.get('coordinate_reference_sequence_sha256','') or ''
  if expected_sha and observed_sha and expected_sha != observed_sha:continue
  if chrom and region.get('chrom','') not in {'',chrom}:continue
  try:
   if int(region['start'])<=pos<=int(region['end']):return region
  except (ValueError,KeyError):pass

def local(r,point):
 if not r:return '.','.','.'
 n=int(r['end'])-int(r['start'])+1; v=point-int(r['start'])+1 if r.get('strand','+')!='-' else int(r['end'])-point+1
 return v,n,v/n

def infer_species_pair_pos_from_human_pair_local(interval,human_paired_local_pos):
 """Deprecated compatibility helper; rRNA matching no longer uses this path."""
 try:
  v=int(human_paired_local_pos); return str(int(interval['end'])-v+1 if interval.get('strand','+')=='-' else int(interval['start'])+v-1)
 except (TypeError,ValueError,KeyError): return '.'

def map_for(directory,sample):
 for suffix in ('.coordinate_map.tsv','.coordinate_map.tsv.gz'):
  path=Path(directory)/f'{sample}{suffix}'
  if path.exists():return load_coordinate_map(path)
 return {}

def structure_match(hclass,sclass):
 h=str(hclass or '').lower();s=str(sclass or '').lower()
 if h not in {'stem','loop'} or s not in {'stem','loop'}:return 'UNKNOWN'
 if h=='stem' and s=='stem':return 'STEM_STEM'
 if h=='loop' and s=='loop':return 'LOOP_LOOP'
 if h=='stem' and s=='loop':return 'STEM_LOOP'
 return 'LOOP_STEM'

def match_tier(status,gm,smatch,pair_relation):
 if smatch=='UNKNOWN':return 'STRUCTURE_UNKNOWN'
 if smatch in {'STEM_LOOP','LOOP_STEM'}:return 'STRUCTURE_DISCORDANT'
 if status=='OK' and gm and smatch=='STEM_STEM' and pair_relation=='yes':return 'HIGH_CONF_STEM'
 if status=='OK' and gm and smatch=='LOOP_LOOP':return 'HIGH_CONF_LOOP'
 return 'LOW_CONF'

def sample_reference_valid(sample_ref,reference_sha,reference_key):
 if not reference_key:return False
 expected=sample_ref.get('coordinate_reference_sequence_sha256','') if sample_ref else ''
 observed=reference_sha.get(reference_key,'')
 if observed=='CONFLICT':return False
 return not (expected and observed and expected!=observed)

def species_structure_for(sample_ref,reference_sha,species_structure,sr,pos):
 if not sample_ref or not sr or not pos:return None
 reference_key=sample_ref.get('reference_key','')
 if not sample_reference_valid(sample_ref,reference_sha,reference_key):return None
 return species_structure.get((reference_key,norm(sr.get('rrna_gene')),pos))

def empty_structural(enabled):
 if not enabled:
  result={f'MTRRNA_{key}':'.' for key in STRUCT}
  result['MTRRNA_MATCH_TIER']='NA'
  return result
 result={f'MTRRNA_{key}':'.' for key in STRUCT}
 result.update({
  'MTRRNA_H_CLASS':'unknown','MTRRNA_S_CLASS':'unknown',
  'MTRRNA_STRUCTURE_MATCH':'UNKNOWN','MTRRNA_PAIR_RELATION_MATCH':'NA',
  'MTRRNA_PAIR_POS_MATCH':'NA','MTRRNA_MATCH_TIER':'STRUCTURE_UNKNOWN',
 })
 return result

def apply_structure(structural,hst,sst,cmap,alt,status,gm):
 hclass=hst.get('struct_class','unknown') if hst else 'unknown'
 sclass=sst.get('struct_class','unknown') if sst else 'unknown'
 smatch=structure_match(hclass,sclass)
 spair=sst.get('paired_pos','.') if sst else '.'
 hpair=hst.get('paired_pos','.') if hst else '.'
 lifted=lift_source_pos_to_human(spair,cmap) if spair not in MISSING else '.'
 pair_relation='NA'
 if smatch=='STEM_STEM' and lifted not in MISSING and hpair not in MISSING:
  pair_relation='yes' if str(lifted)==str(hpair) else 'no'
 h_pair_type=hst.get('pair_type','.') if hst else '.'
 alt_pair_type=rrna_pair_type(alt,hst.get('paired_base')) if hst and hclass=='stem' else '.'
 structural.update({
  'MTRRNA_H_CLASS':hclass,'MTRRNA_S_CLASS':sclass,
  'MTRRNA_H_ELEMENT':hst.get('struct_element','.') if hst else '.',
  'MTRRNA_S_ELEMENT':sst.get('struct_element','.') if sst else '.',
  'MTRRNA_H_PAIR_POS':hpair,'MTRRNA_S_PAIR_POS':spair,
  'MTRRNA_H_PAIR_LOCAL':hst.get('paired_local_pos','.') if hst else '.',
  'MTRRNA_S_PAIR_LOCAL':sst.get('paired_local_pos','.') if sst else '.',
  'MTRRNA_H_PAIR_TYPE':h_pair_type,
  'MTRRNA_S_PAIR_TYPE':sst.get('pair_type','.') if sst else '.',
  'MTRRNA_H_PAIR_STATE':hst.get('pair_state','unknown') if hst else 'unknown',
  'MTRRNA_S_PAIR_STATE':sst.get('pair_state','unknown') if sst else 'unknown',
  'MTRRNA_H_ALT_PAIR_TYPE':alt_pair_type,
  'MTRRNA_H_ALT_EFFECT':rrna_pair_effect(h_pair_type,alt_pair_type) if hclass=='stem' else '.',
  'MTRRNA_S_PAIR_LIFTED_HPOS':lifted,
  'MTRRNA_STRUCTURE_MATCH':smatch,
  'MTRRNA_PAIR_RELATION_MATCH':pair_relation,
  'MTRRNA_PAIR_POS_MATCH':pair_relation,
  'MTRRNA_S_PAIR_EXPECTED_POS':'.',
  'MTRRNA_MATCH_TIER':match_tier(status,gm,smatch,pair_relation),
 })

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--config',required=True);ap.add_argument('--sample');ap.add_argument('--input');ap.add_argument('--output');a=ap.parse_args();c=yaml(a.config);sec=c['rrna_match'];p,s=sec['paths'],sec['settings'];hs=load(p['human_rrna_table']);ss=load_species_rrna_regions(p['species_rrna_table'])
 enabled=bool(s.get('use_rrna_structure_table',False)); hpath=p.get('human_rrna_structure_table',s.get('human_rrna_structure_table','')); spath=p.get('species_rrna_structure_table',s.get('species_rrna_structure_table',''))
 if enabled and (not hpath or not Path(hpath).exists()):raise SystemExit(f'rRNA structure annotation is enabled but human structure table is missing: {hpath or "<unset>"}')
 if enabled and (not spath or not Path(spath).exists()):raise SystemExit(f'rRNA structure annotation is enabled but species/reference structure table is missing: {spath or "<unset>"}')
 human_structure=load_rrna_structure_table(hpath) if enabled else {}; species_structure,reference_sha=load_species_rrna_structure_table(spath) if enabled else ({},{})
 sample_refs=load_sample_reference_map(p.get('sample_reference_map',''))
 samples=[a.sample] if a.sample else sample_names(c)
 if a.input:samples=[a.sample or Path(a.input).name.split('.')[0]]
 allrows=[]
 for sample in samples:
  choices=[(Path(p['input_vcf_dir'])/str(s['input_vcf_pattern']).format(sample=sample),True),(Path(p['fallback_codon_vcf_dir'])/str(s['fallback_codon_vcf_pattern']).format(sample=sample),False),(Path(p['fallback_raw_vcf_dir'])/str(s['fallback_raw_vcf_pattern']).format(sample=sample),False)];inp=Path(a.input) if a.input else next((x for x,_ in choices if x.exists()),None)
  if not inp:raise SystemExit(f'Missing rRNA input VCF for {sample}')
  cmap=map_for(p.get('coordinate_map_dir',''),sample);sample_ref=sample_refs.get(sample,{})
  out=Path(a.output) if a.output else Path(p['output_dir'])/'vcf_rrna'/f"{sample}{s['output_suffix']}";out.parent.mkdir(parents=True,exist_ok=True);head=[];body=[];co=Counter();yes=0;sc=Counter();pc=Counter();tc=Counter();annotated=0;unknown=0
  for line in open_text(inp):
   if line.startswith('#'):head.append(line);continue
   x=line.rstrip().split('\t');inf=info_parse(x[7]);sch,pos,_,_=source(inf);hp=human_pos(x,inf);sr=species_region_for(sample_ref,ss,pos,sch) if pos else None;hr=hit(hs,hp,x[0]) if hp else None
   status='MISSING_COORD' if not pos or not hp else 'NO_SPECIES_OR_HUMAN_RRNA' if not sr and not hr else 'NO_SPECIES_RRNA' if not sr else 'NO_HUMAN_RRNA' if not hr else 'GENE_MISMATCH' if norm(sr.get('rrna_gene'))!=norm(hr.get('rrna_gene')) else 'OK'
   sl,slen,sf=local(sr,pos) if sr else ('.','.','.');hl,hlen,hf=local(hr,hp) if hr else ('.','.','.');gm=bool(sr and hr and norm(sr.get('rrna_gene'))==norm(hr.get('rrna_gene')));strand=bool(sr and hr and sr.get('strand','+')==hr.get('strand','+'));region=gm and (not s.get('require_same_strand',False) or strand);yes+=region
   v={'MTRRNA_STATUS':status,'MTRRNA_S_GENE':norm(sr.get('rrna_gene')) if sr else '.','MTRRNA_H_GENE':norm(hr.get('rrna_gene')) if hr else '.','MTRRNA_GENE_MATCH':'yes' if gm else 'no','MTRRNA_S_LOCAL':sl,'MTRRNA_H_LOCAL':hl,'MTRRNA_S_LEN':slen,'MTRRNA_H_LEN':hlen,'MTRRNA_S_FRAC':sf,'MTRRNA_H_FRAC':hf,'MTRRNA_FRAC_DELTA':abs(sf-hf) if isinstance(sf,float) and isinstance(hf,float) else '.','MTRRNA_STRAND_MATCH':'yes' if strand else 'no','MTRRNA_REGION_MATCH':'yes' if region else 'no'}
   structural=empty_structural(enabled);structural['MTRRNA_LOCAL_MATCH']=compare_values(sl,hl)
   if enabled:
    hst=human_structure.get((norm(hr.get('rrna_gene')) if hr else '',hp)) if hp else None
    sst=species_structure_for(sample_ref,reference_sha,species_structure,sr,pos)
    apply_structure(structural,hst,sst,cmap,x[4],status,gm)
    smatch=structural['MTRRNA_STRUCTURE_MATCH'];tier=structural['MTRRNA_MATCH_TIER'];relation=structural['MTRRNA_PAIR_RELATION_MATCH']
    sc[smatch]+=1;tc[tier]+=1
    if relation=='yes':pc['yes']+=1
    elif relation=='no':pc['no']+=1
    if smatch=='UNKNOWN':unknown+=1
    else:annotated+=1
   v.update(structural);inf.update(v);x[7]=info_format(inf);body.append('\t'.join(x)+'\n');co[status]+=1
  with out.open('w') as f:f.writelines(inject_headers(head,FIELDS,'MTRRNA'));f.writelines(body)
  row={'sample':sample,'input_vcf':str(inp),'output_vcf':str(out),'total_records':len(body),**{f'status_{q}':co[q] for q in ['OK','NO_SPECIES_RRNA','NO_HUMAN_RRNA','NO_SPECIES_OR_HUMAN_RRNA','GENE_MISMATCH','MISSING_COORD']},'rrna_region_match_yes':yes,'rrna_region_match_no':len(body)-yes,'rrna_structure_annotated':annotated,'rrna_structure_unknown':unknown,'n_stem_stem':sc['STEM_STEM'],'n_loop_loop':sc['LOOP_LOOP'],'n_stem_loop':sc['STEM_LOOP'],'n_loop_stem':sc['LOOP_STEM'],'n_pair_relation_match':pc['yes'],'n_pair_relation_mismatch':pc['no'],'n_high_conf_stem':tc['HIGH_CONF_STEM'],'n_high_conf_loop':tc['HIGH_CONF_LOOP'],'n_structure_discordant':tc['STRUCTURE_DISCORDANT'],'n_structure_unknown':tc['STRUCTURE_UNKNOWN'],'status':'completed'};write_summary(Path(p['reports_dir'])/f'{sample}.rrna_match_summary.tsv',row);allrows.append(row)
 if allrows and not a.sample and not a.input:
  q=Path(p['reports_dir'])/'all_samples.rrna_match_summary.tsv';q.parent.mkdir(parents=True,exist_ok=True)
  with q.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(allrows[0]),delimiter='\t');w.writeheader();w.writerows(allrows)
if __name__=='__main__':main()
