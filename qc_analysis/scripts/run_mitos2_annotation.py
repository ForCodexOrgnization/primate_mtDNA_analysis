#!/usr/bin/env python3
"""Run MITOS2 once per final chrM FASTA and build sample-level interval/codon tables.

The materialized chrM FASTA is always the coordinate authority; MITOS2 supplies
features only.  MITOS2 command-line interfaces vary, so the invocation is kept
in one deliberately small, inspectable command builder.
"""
import argparse, csv, re, shlex, shutil, subprocess, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from qc_analysis.lib.match_utils import yaml
try:
    from Bio import SeqIO
except ImportError:
    SeqIO = None

FEATURE_FIELDS = 'reference_key reference_species coordinate_reference_accession coordinate_reference_fasta feature_type gene gene_raw start end strand score source_file annotation_source'.split()
CODON_FIELDS = 'file_name seq_name sample species species_key accession accession_version reference_id family pos ref_base_genome gene gene_raw product protein_id strand codon_index codon_pos_in_triplet codon_seq codon_pos1_genomic codon_pos2_genomic codon_pos3_genomic codon_start_qualifier transl_table cds_tail_incomplete_bases annotation_source coordinate_reference_fasta coordinate_reference_accession'.split()
SUMMARY_FIELDS = 'reference_key reference_species coordinate_reference_accession coordinate_reference_fasta status mitos2_command raw_dir n_features n_cds_features n_coding_position_rows note'.split()
GENES = {'ND1':'MT-ND1','ND2':'MT-ND2','ND3':'MT-ND3','ND4':'MT-ND4','ND4L':'MT-ND4L','ND5':'MT-ND5','ND6':'MT-ND6','COX1':'MT-CO1','COI':'MT-CO1','COX2':'MT-CO2','COII':'MT-CO2','COX3':'MT-CO3','COIII':'MT-CO3','CYTB':'MT-CYB','ATP6':'MT-ATP6','ATP8':'MT-ATP8'}
def val(r,k): return (r.get(k) or '').strip()
def sk(s): return re.sub(r'_+','_',re.sub(r'\s+','_',s.lower())).strip('_')
def norm(g):
    k=re.sub(r'[^A-Z0-9]','',g.upper().replace('MT','',1)); return GENES.get(k,g)
def read(path):
    if not Path(path).exists(): return []
    with open(path,newline='') as h:
        rows=list(csv.reader(h,delimiter='\t'))
    if not rows:return []
    return [dict(zip(rows[0],r)) for r in rows[1:]] if 'sample' in rows[0] else [{'sample':r[0],'species':r[1] if len(r)>1 else ''} for r in rows]
def write(path, fields, rows):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    with open(path,'w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
def attrs(text):
    result={}
    for x in text.split(';'):
        if '=' in x: k,v=x.split('=',1);result[k.lower()]=v
        elif ' ' in x: k,v=x.split(' ',1);result[k.lower()]=v.strip('"')
    return result
def typ(x):
    x=x.upper()
    return 'CDS' if x=='CDS' else 'tRNA' if x in ('TRNA','T-RNA') else 'rRNA' if x in ('RRNA','R-RNA') else ''
def parse_outputs(directory, ref):
    rows=[]; seen=set()
    for p in sorted(Path(directory).rglob('*')):
        if not p.is_file() or p.suffix.lower() not in ('.gff','.gff3','.bed','.tsv','.txt','.tbl'): continue
        for line in p.read_text(errors='replace').splitlines():
            if not line or line.startswith('#'):continue
            c=line.split('\t'); ft=''; start=end=''; strand='+'; score=''; raw=''
            if len(c)>=9: # GFF
                ft=typ(c[2]); start,end,strand,score,raw=c[3],c[4],c[6],c[5],attrs(c[8]).get('gene') or attrs(c[8]).get('name') or attrs(c[8]).get('product') or c[2]
            elif len(c)>=3:
                # BED is zero based; other simple tables are normally one based. MITOS BED uses BED semantics.
                ft=typ(c[3] if len(c)>3 else ''); start=str(int(c[1])+1) if p.suffix.lower()=='.bed' else c[1]; end=c[2]; raw=c[3] if len(c)>3 else ft; strand=c[5] if len(c)>5 else '+'
            if not ft: continue
            key=(ft,start,end,strand,raw)
            if key in seen:continue
            seen.add(key); rows.append({**ref,'feature_type':ft,'gene':norm(raw) if ft=='CDS' else raw,'gene_raw':raw,'start':start,'end':end,'strand':strand or '+','score':score,'source_file':str(p),'annotation_source':'MITOS2'})
    return rows
def detect(settings):
    candidates=[x.strip() for x in str(settings.get('mitos2_command_candidates','mitos2,mitos,runmitos.py')).split(',')]
    activation=f"module load {shlex.quote(str(settings.get('conda_module','miniconda/24.11.3')))} && source \"$(conda info --base)/etc/profile.d/conda.sh\" && conda activate {shlex.quote(str(settings.get('conda_env','mitos2')))}"
    probe='; '.join(f'command -v {shlex.quote(x)} && exit 0' for x in candidates) + '; exit 1'
    x=subprocess.run(['bash','-lc',activation+' && '+probe],text=True,capture_output=True)
    if x.returncode: raise RuntimeError(f"ERROR: MITOS2 was not found after activating conda env {settings.get('conda_env','mitos2')}.")
    return x.stdout.strip().splitlines()[-1], activation
def codons(features, fasta, ref, samples, code):
    record=next(SeqIO.parse(str(fasta),'fasta')); seq=str(record.seq).upper(); result=[]
    for f in features:
        if f['feature_type']!='CDS':continue
        a,b=int(f['start'])-1,int(f['end']); coords=list(range(a,b)); strand=f['strand']
        if strand=='-':
            coords.reverse(); dna=''.join(seq[x] for x in coords).translate(str.maketrans('ACGTN','TGCAN'))
        else:dna=''.join(seq[x] for x in coords)
        usable=len(dna)//3*3
        for i in range(0,usable,3):
            trip=coords[i:i+3]; cod=dna[i:i+3]
            for phase,pos in enumerate(trip,1):
                base=seq[pos]
                result.append({'file_name':Path(fasta).name,'seq_name':record.id,'sample':'','species':'','species_key':'','accession':ref['coordinate_reference_accession'],'accession_version':ref['coordinate_reference_accession'],'reference_id':ref['coordinate_reference_accession'],'family':'','pos':pos+1,'ref_base_genome':base,'gene':f['gene'],'gene_raw':f['gene_raw'],'product':f['gene_raw'],'protein_id':'','strand':strand,'codon_index':i//3+1,'codon_pos_in_triplet':phase,'codon_seq':cod,'codon_pos1_genomic':trip[0]+1,'codon_pos2_genomic':trip[1]+1,'codon_pos3_genomic':trip[2]+1,'codon_start_qualifier':'1','transl_table':code,'cds_tail_incomplete_bases':len(dna)-usable,'annotation_source':'MITOS2','coordinate_reference_fasta':str(fasta),'coordinate_reference_accession':ref['coordinate_reference_accession']})
    expanded=[]
    for s in samples:
        for r in result: expanded.append({**r,'sample':s['sample'],'species':s['species'],'species_key':sk(s['species']),'accession':s.get('accession',r['accession']),'reference_id':s.get('accession',r['reference_id'])})
    return expanded
def main():
 p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--sample');p.add_argument('--reference');p.add_argument('--force',action='store_true');p.add_argument('--dry-run',action='store_true');a=p.parse_args()
 sec=yaml(a.config).get('mitos2_annotation');
 if not sec: raise SystemExit('Missing mitos2_annotation section in config.')
 paths,settings=sec['paths'],sec.get('settings',{});
 if not settings.get('enabled',True):print('MITOS2 annotation disabled.');return
 manifest=read(paths['reference_manifest']); samples=read(paths['sample_ref_file']);
 refs={}
 for m in manifest:
  species=val(m,'final_chrM_species') or val(m,'target_species'); fasta=val(m,'chrM_expected_output_fasta') or str(Path(paths['fasta_dir'])/(species+'.fa')); accession=val(m,'final_chrM_accession') or val(m,'final_chrM_genbank_accn') or val(m,'final_chrM_refseq_accn'); key=re.sub(r'[^A-Za-z0-9_.-]+','_',accession or Path(fasta).stem)
  if a.reference and a.reference not in (species,accession,key):continue
  refs.setdefault(str(Path(fasta)),{'reference_key':key,'reference_species':species,'coordinate_reference_accession':accession,'coordinate_reference_fasta':str(Path(fasta))})
 allf=[];allc=[];summary=[]
 command=activation=''
 for fasta,ref in refs.items():
  linked=[{'sample':val(s,'sample'),'species':val(s,'species'),'accession':ref['coordinate_reference_accession']} for s in samples if sk(val(s,'species')) in (sk(ref['reference_species']), sk(next((val(m,'target_species') for m in manifest if val(m,'final_chrM_species')==ref['reference_species']),'')))]
  if a.sample:linked=[x for x in linked if x['sample']==a.sample]
  if not linked and a.sample:continue
  raw=Path(paths['mitos2_raw_dir'])/ref['reference_key']; status='completed';note=''
  try:
   if not Path(fasta).exists(): raise FileNotFoundError(f'Final chrM FASTA is missing: {fasta}')
   if a.dry_run: command='dry-run'; feats=[];status='dry_run';note='Would run MITOS2.'
   else:
    command,activation=detect(settings)
    if a.force or settings.get('overwrite_existing',False) or not raw.exists() or not any(raw.iterdir()):
     raw.mkdir(parents=True,exist_ok=True); run=f'{activation} && {shlex.quote(command)} -i {shlex.quote(fasta)} -o {shlex.quote(str(raw))} --code {settings.get("genetic_code",2)} --threads {settings.get("threads",4)}'
     x=subprocess.run(['bash','-lc',run],text=True,capture_output=True)
     if x.returncode: raise RuntimeError(x.stderr.strip() or 'MITOS2 failed')
    feats=parse_outputs(raw,ref)
   allf.extend(feats); allc.extend(codons(feats,fasta,ref,linked,str(settings.get('genetic_code',2))) if feats and not a.dry_run else [])
  except Exception as e: status='failed';note=str(e)
  summary.append({**ref,'status':status,'mitos2_command':command,'raw_dir':str(raw),'n_features':len([x for x in allf if x['reference_key']==ref['reference_key']]),'n_cds_features':len([x for x in allf if x['reference_key']==ref['reference_key'] and x['feature_type']=='CDS']),'n_coding_position_rows':len([x for x in allc if x['reference_id']==ref['coordinate_reference_accession']]),'note':note})
 write(paths['mitos2_feature_table'],FEATURE_FIELDS,allf);write(paths['mitos2_cds_table'],CODON_FIELDS,allc);write(paths['mitos2_summary_table'],SUMMARY_FIELDS,summary)
 if any(x['status']=='failed' for x in summary): raise SystemExit('One or more MITOS2 references failed; see summary table.')
 print(f'Wrote {len(allf)} features and {len(allc)} sample-level coding rows.')
if __name__=='__main__':main()
