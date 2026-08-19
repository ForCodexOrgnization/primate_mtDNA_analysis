#!/usr/bin/env python3
"""Consolidate orthology annotations and build a pre-human-QC homo background."""
from __future__ import annotations
import argparse, csv, gzip, json, math, sys
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from qc_analysis.lib.simple_yaml import read_simple_yaml

ORTHO_COLUMNS="sample species human_chrom human_pos human_ref human_alt source_chrom source_pos source_ref source_alt AF DP variant_class region_type orthology_match_status orthology_fail_reason".split()
HOMO_COLUMNS="sample species genus family human_pos human_ref human_alt AF DP region_type orthology_match_status".split()

def resolve(v):
 p=Path(str(v)).expanduser();return p if p.is_absolute() else ROOT/p
def parse_info(s):return {x.split('=',1)[0]:x.split('=',1)[1] if '=' in x else True for x in s.split(';') if x and x!='.'}
def number(v):
 try:
  x=float(v);return x if math.isfinite(x) else None
 except (TypeError,ValueError):return None
def metadata(path):
 with path.open(newline='',encoding='utf-8-sig') as h:
  rows=[(line_number,row) for line_number,row in enumerate(csv.reader(h,delimiter='\t'),1)
        if row and any(value.strip() for value in row) and not row[0].lstrip().startswith('#')]
 if not rows:return {}
 _first_line,first=rows[0];header=[value.strip().lower() for value in first]
 headered=bool(set(header)&{'sample','species','species_fasta','genus','family'})
 if headered and 'sample' not in header:
  raise ValueError(f"Metadata file {path} has a header but no sample column")
 if headered and len(set(header))!=len(header):
  raise ValueError(f"Metadata file {path} has duplicate column names after case normalization")
 data=rows[1:] if headered else rows;result={};sample_lines={}
 for line_number,fields in data:
  if headered:
   row={name:(fields[index].strip() if index<len(fields) else '') for index,name in enumerate(header)}
  else:
   if len(fields)<2:
    raise ValueError(f"Headerless metadata file {path} line {line_number} must contain at least sample and species columns")
   row={'sample':fields[0].strip(),'species':fields[1].strip()}
  sample=row.get('sample','').strip()
  if not sample:continue
  if sample in result:
   raise ValueError(f"Duplicate sample ID {sample!r} in metadata file {path} at lines {sample_lines[sample]} and {line_number}")
  result[sample]=row;sample_lines[sample]=line_number
 return result
def status(info):
 """Expose existing match decisions without changing their biological rules."""
 cs=str(info.get('MTCODON_STATUS',''))
 ts=str(info.get('MTTRNA_STATUS',''))
 rs=str(info.get('MTRRNA_STATUS',''))
 if cs and cs!='SKIPPED_NONCODING':
  if cs=='PASS':return 'CDS','PASS',''
  if 'AMBIGUOUS' in cs:return 'CDS','AMBIGUOUS',cs
  return 'CDS','FAIL',cs
 if ts and ts!='NO_SPECIES_OR_HUMAN_TRNA':
  if ts=='OK' and info.get('MTTRNA_STRICT_MATCH')=='yes':return 'tRNA','PASS',''
  if 'AMBIGUOUS' in ts:return 'tRNA','AMBIGUOUS',ts
  return 'tRNA','FAIL',ts if ts!='OK' else 'STRICT_MATCH_NO'
 if rs and rs!='NO_SPECIES_OR_HUMAN_RRNA':
  if rs=='OK' and info.get('MTRRNA_REGION_MATCH')=='yes':return 'rRNA','PASS',''
  return 'rRNA','FAIL',rs if rs!='OK' else 'REGION_MATCH_NO'
 return 'noncoding','NOT_APPLICABLE','NOT_APPLICABLE_NONCODING'
def variants(path):
 op=gzip.open if path.suffix=='.gz' else open
 with op(path,'rt',encoding='utf-8') as h:
  for line in h:
   if line.startswith('#'):continue
   f=line.rstrip().split('\t');info=parse_info(f[7]);fmt=f[8].split(':') if len(f)>8 else []; vals=f[9].split(':') if len(f)>9 else []; sample=dict(zip(fmt,vals))
   af=number(sample.get('AF','').split(',')[0]);dp=number(sample.get('DP','')) or number(info.get('DP'))
   if af is None:
    ad=[number(x) for x in sample.get('AD','').split(',')]
    if len(ad)==2 and None not in ad and sum(ad)>0:af=ad[1]/sum(ad)
   yield f,info,af,dp,sample
def marker_alleles(path):
 if not path.is_file():return set()
 with path.open(newline='',encoding='utf-8') as h:
  out=set()
  for r in csv.DictReader(h,delimiter='\t'):
   try:out.add((int(r.get('pos') or r.get('human_pos')),str(r.get('ref') or r.get('human_ref') or '').upper(),str(r.get('alt') or r.get('human_alt') or '').upper()))
   except (TypeError,ValueError):pass
  return out
def write(path,columns,rows):
 path.parent.mkdir(parents=True,exist_ok=True)
 with path.open('w',newline='',encoding='utf-8') as h:w=csv.DictWriter(h,fieldnames=columns,delimiter='\t');w.writeheader();w.writerows(rows)
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--config',required=True);a=ap.parse_args();cfg=read_simple_yaml(a.config);sec=cfg.get('primate_homo_background') or {}; paths=sec.get('paths') or {}; settings=sec.get('settings') or {}
 if sec.get('enabled',True) is False:return 0
 indir=resolve(paths.get('input_vcf_dir','results/qc/rrna_match/vcf_rrna'));out=resolve(paths.get('output_dir','results/qc/primate_homo_background'));orthodir=resolve(paths.get('orthology_reports_dir','results/qc/orthology_match/reports'));meta=metadata(resolve(paths.get('sample_ref_file','config/sample_ref_file.tsv')))
 pattern=paths.get('input_vcf_pattern','{sample}.lifted.codon.trna.rrna.vcf'); homo=float(settings.get('homoplasmy_af_min',.95));dpmin=float(settings.get('dp_min',100));pass_only=settings.get('pass_only',True) is not False;snv_only=settings.get('snv_only',True) is not False;accepted=set(str(settings.get('accepted_orthology_statuses','PASS')).split(','))
 markers=marker_alleles(resolve(paths.get('human_marker_table','data/reference_tables/human_phylotree_rcrs_v17.1_snv.tsv'))); ortho=[];homos=[];eligible_by_species=defaultdict(set);processed_samples=set()
 for sample,m in sorted(meta.items()):
  path=indir/pattern.format(sample=sample)
  if not path.is_file() and Path(str(path)+'.gz').is_file():path=Path(str(path)+'.gz')
  if not path.is_file():continue
  processed_samples.add(sample);eligible_by_species[m.get('species','')].add(sample)
  for f,info,af,dp,_ in variants(path):
   region,match,reason=status(info);srcpos=info.get('SRC_POS',info.get('MTLIFT_ORIG_POS',''));row=dict(sample=sample,species=m.get('species',''),human_chrom=f[0],human_pos=f[1],human_ref=f[3],human_alt=f[4],source_chrom=info.get('SRC_CHROM',info.get('MTLIFT_ORIG_CHROM','')),source_pos=srcpos,source_ref=info.get('SRC_REF',info.get('MTLIFT_ORIG_REF','')),source_alt=info.get('SRC_ALT',info.get('MTLIFT_ORIG_ALT','')),AF=af if af is not None else 'NA',DP=dp if dp is not None else 'NA',variant_class='SNV' if len(f[3])==len(f[4])==1 and ',' not in f[4] else 'OTHER',region_type=region,orthology_match_status=match,orthology_fail_reason=reason);ortho.append(row)
   valid=len(f[3])==len(f[4])==1 and f[3] in 'ACGT' and f[4] in 'ACGT' and ',' not in f[4]
   if (not snv_only or valid) and (not pass_only or f[6]=='PASS') and af is not None and af>=homo and dp is not None and dp>=dpmin and match in accepted:
    homos.append({k:row[k] for k in ('sample','species','human_pos','human_ref','human_alt','AF','DP','region_type','orthology_match_status')}|{'genus':m.get('genus',''),'family':m.get('family','')})
 write(orthodir/'orthology_match_report.tsv',ORTHO_COLUMNS,ortho);write(out/'primate_homo_background.tsv',HOMO_COLUMNS,homos)
 grouped=defaultdict(list)
 for r in homos:grouped[(r['human_pos'],r['human_ref'],r['human_alt'])].append(r)
 collapsed=[]
 for key,rs in sorted(grouped.items(),key=lambda x:(int(x[0][0]),x[0][1:])):
  allele=(int(key[0]),key[1],key[2]);carriers=len({r['sample'] for r in rs});collapsed.append(dict(human_pos=key[0],human_ref=key[1],human_alt=key[2],n_primate_samples=carriers,n_primate_species=len({r['species'] for r in rs if r['species']}),n_primate_genera=len({r['genus'] for r in rs if r['genus']}),n_primate_families=len({r['family'] for r in rs if r['family']}),species_list=';'.join(sorted({r['species'] for r in rs if r['species']})),background_frequency=carriers/len(processed_samples) if processed_samples else 'NA',is_human_phylotree_marker=str(allele in markers).lower()))
 write(out/'primate_homo_marker_background.tsv','human_pos human_ref human_alt n_primate_samples n_primate_species n_primate_genera n_primate_families species_list background_frequency is_human_phylotree_marker'.split(),collapsed)
 species_rows=[]
 for (pos,ref,alt),rs in grouped.items():
  for species in sorted({r['species'] for r in rs}):
   carriers=len({r['sample'] for r in rs if r['species']==species});n=len(eligible_by_species[species]);species_rows.append(dict(species=species,human_pos=pos,human_ref=ref,human_alt=alt,n_eligible_samples=n,n_homo_carriers=carriers,background_frequency=carriers/n if n else 'NA'))
 write(out/'species_marker_background.tsv','species human_pos human_ref human_alt n_eligible_samples n_homo_carriers background_frequency'.split(),species_rows)
 provenance={'homoplasmy_af_threshold':homo,'dp_threshold':dpmin,'vcf_filter_requirement':'PASS' if pass_only else 'ANY','snv_only':snv_only,'orthology_acceptance_criteria':sorted(accepted),'input_vcf_directory':str(indir),'human_marker_reference_version':settings.get('human_marker_reference_version','rcrs-v17.1'),'number_of_samples':len(processed_samples),'number_of_species':len({meta[s].get('species','') for s in processed_samples}),'number_of_distinct_homo_alleles':len(grouped)}
 (out/'primate_homo_background_metadata.json').write_text(json.dumps(provenance,indent=2)+'\n');print(f'[primate_homo_background] orthology={len(ortho)} homo={len(homos)} alleles={len(grouped)}')
if __name__=='__main__':main()
