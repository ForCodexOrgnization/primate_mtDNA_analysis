#!/usr/bin/env python3
"""Run MITOS2 per final chrM FASTA and materialize diagnostic annotation tables."""
import argparse, csv, re, shlex, subprocess, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from qc_analysis.lib.match_utils import yaml
try:
 from Bio import SeqIO
except ImportError: SeqIO = None
FEATURE_FIELDS='reference_key reference_species coordinate_reference_accession coordinate_reference_fasta feature_type gene gene_raw start end strand score source_file annotation_source'.split()
CODON_FIELDS='file_name seq_name sample species species_key accession accession_version reference_id family pos ref_base_genome gene gene_raw product protein_id strand codon_index codon_pos_in_triplet codon_seq codon_pos1_genomic codon_pos2_genomic codon_pos3_genomic codon_start_qualifier transl_table cds_tail_incomplete_bases annotation_source coordinate_reference_fasta coordinate_reference_accession'.split()
SUMMARY_FIELDS='reference_key reference_species coordinate_reference_accession coordinate_reference_fasta status command_mode mitos2_command attempted_commands return_code stdout_log stderr_log help_log raw_dir n_features n_cds_features n_coding_position_rows n_output_files_scanned n_parseable_files parser_status note'.split()
DIAG_FIELDS='reference_key file suffix n_lines n_candidate_feature_lines parser_used n_features_parsed'.split()
GENES={'ND1':'MT-ND1','ND2':'MT-ND2','ND3':'MT-ND3','ND4':'MT-ND4','ND4L':'MT-ND4L','ND5':'MT-ND5','ND6':'MT-ND6','COX1':'MT-CO1','COI':'MT-CO1','COX2':'MT-CO2','COII':'MT-CO2','COX3':'MT-CO3','COIII':'MT-CO3','CYTB':'MT-CYB','ATP6':'MT-ATP6','ATP8':'MT-ATP8'}
CODING=set(GENES)
def val(r,k): return (r.get(k) or '').strip()
def sk(s): return re.sub(r'_+','_',re.sub(r'\s+','_',s.lower())).strip('_')
def norm(g):
 k=re.sub(r'[^A-Z0-9]','',g.upper().replace('MT','',1)); return GENES.get(k,g)
def write(p,fields,rows):
 Path(p).parent.mkdir(parents=True,exist_ok=True)
 with open(p,'w',newline='') as h:
  w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
def read(p):
 p=Path(p)
 if not p.exists():return []
 with p.open(newline='') as h: rows=[x for x in csv.reader(h,delimiter='\t') if any(y.strip() for y in x)]
 if not rows:return []
 headers={'sample','target_species','final_chrM_species','final_chrM_accession','chrM_expected_output_fasta'}
 return [dict(zip(rows[0],x)) for x in rows[1:]] if headers.intersection(rows[0]) else [{'sample':x[0].strip(),'species':x[1].strip() if len(x)>1 else ''} for x in rows]
def attrs(s):
 d={}
 for x in s.split(';'):
  if '=' in x:k,v=x.split('=',1);d[k.lower()]=v.strip('"')
  elif ' ' in x:k,v=x.split(' ',1);d[k.lower()]=v.strip(' "')
 return d
def infer(raw,declared=''):
 t=declared.lower(); name=re.sub(r'[^a-z0-9]','',raw.lower())
 if t.lower() in ('cds','trna','rrna'): return {'cds':'CDS','trna':'tRNA','rrna':'rRNA'}[t]
 if t=='gene':
  if name in {x.lower() for x in CODING}:return 'CDS'
  if name.startswith(('trn','trna')):return 'tRNA'
  if any(x in name for x in ('rrns','rrnl','12s','16s','rrna')):return 'rRNA'
 return ''
def text_file(p):
 if p.suffix.lower() in ('.fa','.fasta','.fna','.gz','.bam','.png','.pdf'):return False
 try: p.read_text(errors='strict');return True
 except (UnicodeDecodeError,OSError):return False
def parse_outputs(raw,ref):
 features=[]; diagnostics=[]; seen=set(); allowed={'.gff','.gff3','.bed','.tbl','.tsv','.txt','.out','.result','.mitos',''}
 for p in sorted(Path(raw).rglob('*')):
  if not p.is_file() or p.name.startswith('mitos2.') or p.suffix.lower() not in allowed or not text_file(p):continue
  lines=p.read_text(errors='replace').splitlines(); parsed=[]; cand=0; parser='none'
  for line in lines:
   if not line or line.startswith('#'):continue
   c=line.split('\t'); ft=rawgene=''; start=end=strand=score=''
   if len(c)>=9 and c[3].isdigit() and c[4].isdigit():
    parser='gff'; at=attrs(c[8]); rawgene=at.get('gene') or at.get('name') or at.get('product') or at.get('id') or ''; ft=infer(rawgene,c[2]); start,end,strand,score=c[3],c[4],c[6] or '+',c[5]; cand+=bool(ft)
   elif len(c)>=3 and c[1].isdigit() and c[2].isdigit():
    parser='tabular'; rawgene=c[3] if len(c)>3 else ''; ft=infer(rawgene,rawgene); start=str(int(c[1])+1) if p.suffix.lower()=='.bed' else c[1];end=c[2];strand=c[5] if len(c)>5 else '+';cand+=bool(ft)
   if not ft:continue
   key=(ft,start,end,strand,rawgene)
   if key in seen:continue
   seen.add(key); parsed.append({**ref,'feature_type':ft,'gene':norm(rawgene) if ft=='CDS' else rawgene,'gene_raw':rawgene,'start':start,'end':end,'strand':strand,'score':score,'source_file':str(p),'annotation_source':'MITOS2'})
  features+=parsed;diagnostics.append({'reference_key':ref['reference_key'],'file':str(p),'suffix':p.suffix,'n_lines':len(lines),'n_candidate_feature_lines':cand,'parser_used':parser,'n_features_parsed':len(parsed)})
 return features,diagnostics
def activate(settings):
 # mitos2 is the conda environment, mitos is the package, and runmitos is the CLI.
 return f"module load {shlex.quote(str(settings.get('conda_module', 'miniconda')))} && source \"$(conda info --base)/etc/profile.d/conda.sh\" && conda activate {shlex.quote(str(settings.get('conda_env', 'mitos2')))}"
def command(settings):
 validation = activate(settings) + '''
if ! command -v runmitos >/dev/null 2>&1; then
    echo "ERROR: runmitos was not found after activating conda env mitos2." >&2
    echo "CONDA_PREFIX=${CONDA_PREFIX:-not_set}" >&2
    echo "PATH=$PATH" >&2
    exit 1
fi

echo "CONDA_PREFIX=$CONDA_PREFIX"
echo "MITOS executable=$(command -v runmitos || true)"
echo "Using MITOS2 executable: $(command -v runmitos)"
'''
 x=subprocess.run(['bash','-lc',validation],text=True,capture_output=True)
 if x.returncode != 0:
  raise RuntimeError(x.stderr.strip() or 'ERROR: runmitos validation failed after conda activation.')
 return 'runmitos',x.stdout
def templates(exe,fasta,out,settings):
 q=shlex.quote
 exe=str(exe);fasta=str(fasta);out=str(out)
 code=str(settings.get('genetic_code',2))
 refseqver=str(settings.get('refseqver','refseq81m'))
 refdir=str(settings.get('refdir','') or '')
 common=f'-c {q(code)} -o {q(out)} -r {q(refseqver)} --best --noplots'
 if refdir:common+=f' -R {q(refdir)}'
 return [f'{q(exe)} --fasta {q(fasta)} {common}',f'{q(exe)} -i {q(fasta)} {common}']
def codons(features,fasta,ref,samples,code):
 if SeqIO is None:raise RuntimeError('Biopython is required to create MITOS2 codon rows.')
 rec=next(SeqIO.parse(str(fasta),'fasta'));seq=str(rec.seq).upper();base=[]
 for f in features:
  if f['feature_type']!='CDS':continue
  coords=list(range(int(f['start'])-1,int(f['end'])));strand=f['strand']; dna=''.join(seq[i] for i in coords)
  if strand=='-':coords.reverse();dna=dna.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
  usable=len(dna)//3*3
  for i in range(0,usable,3):
   trip=coords[i:i+3]
   for phase,pos in enumerate(trip,1):base.append({'file_name':Path(fasta).name,'seq_name':rec.id,'sample':'','species':'','species_key':'','accession':ref['coordinate_reference_accession'],'accession_version':ref['coordinate_reference_accession'],'reference_id':ref['coordinate_reference_accession'],'family':'','pos':pos+1,'ref_base_genome':seq[pos],'gene':f['gene'],'gene_raw':f['gene_raw'],'product':f['gene_raw'],'protein_id':'','strand':strand,'codon_index':i//3+1,'codon_pos_in_triplet':phase,'codon_seq':dna[i:i+3],'codon_pos1_genomic':trip[0]+1,'codon_pos2_genomic':trip[1]+1,'codon_pos3_genomic':trip[2]+1,'codon_start_qualifier':'1','transl_table':code,'cds_tail_incomplete_bases':len(dna)-usable,'annotation_source':'MITOS2','coordinate_reference_fasta':str(fasta),'coordinate_reference_accession':ref['coordinate_reference_accession']})
 return [{**r,'sample':s['sample'],'species':s['species'],'species_key':sk(s['species'])} for s in samples for r in base]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--config',required=True);ap.add_argument('--sample');ap.add_argument('--reference');ap.add_argument('--force',action='store_true');ap.add_argument('--dry-run',action='store_true');a=ap.parse_args();sec=yaml(a.config).get('mitos2_annotation');
 if not sec:raise SystemExit('Missing mitos2_annotation section in config.')
 paths,settings=sec['paths'],sec.get('settings',{});manifest=read(paths['reference_manifest']);samples=read(paths['sample_ref_file']);
 if a.sample:samples=[s for s in samples if val(s,'sample')==a.sample]
 refs={}
 for m in manifest:
  target=val(m,'target_species'); species=val(m,'final_chrM_species') or target; fasta=val(m,'chrM_expected_output_fasta') or str(Path(paths['fasta_dir'])/(species+'.fa')); acc=val(m,'final_chrM_accession') or val(m,'final_chrM_genbank_accn') or val(m,'final_chrM_refseq_accn');key=re.sub(r'[^A-Za-z0-9_.-]+','_',acc or Path(fasta).stem)
  if a.reference and a.reference not in (target,species,acc,key):continue
  r=refs.setdefault(str(Path(fasta)),{'reference_key':key,'reference_species':species,'coordinate_reference_accession':acc,'coordinate_reference_fasta':str(Path(fasta)),'targets':set()});r['targets'].add(sk(target))
 allf=[];allc=[];summ=[]
 for fasta,ref in refs.items():
  linked=[{'sample':val(s,'sample'),'species':val(s,'species')} for s in samples if sk(val(s,'species')) in ref['targets'] or sk(val(s,'species'))==sk(ref['reference_species'])]
  if a.sample and not linked:continue
  raw=Path(paths['mitos2_raw_dir'])/ref['reference_key'];raw.mkdir(parents=True,exist_ok=True);logs={x:str(raw/f'mitos2.{x}.txt') for x in ('command','stdout','stderr','returncode','help')};attempted=[];rc='';note='';status='failed';exe='';mode=str(settings.get('mitos2_command_mode','auto'))
  # Create every diagnostic log even if environment activation or command discovery fails.
  for path in logs.values(): Path(path).write_text('')
  try:
   if not Path(fasta).exists():raise FileNotFoundError(f'Final chrM FASTA is missing: {fasta}')
   exe,validation=command(settings);mode='runmitos'; helpx=subprocess.run(['bash','-lc',activate(settings)+f' && {shlex.quote(exe)} --help'],text=True,capture_output=True);Path(logs['help']).write_text(validation+helpx.stdout+'\n'+helpx.stderr)
   marker=raw/'mitos2.completed.ok'
   if a.dry_run:Path(logs['command']).write_text('dry-run\n');Path(logs['stdout']).write_text('');Path(logs['stderr']).write_text('');Path(logs['returncode']).write_text('0\n');status='dry_run';features=[];diag=[]
   elif marker.exists() and not a.force:
    # Only a validated successful marker permits reuse; incomplete raw directories rerun.
    rc='0'; attempted=['reused successful output']; features,diag=parse_outputs(raw,ref);write(raw/'parsed_output_files.tsv',DIAG_FIELDS,diag)
    cds=[f for f in features if f['feature_type']=='CDS']; rows=codons(features,fasta,ref,linked,str(settings.get('genetic_code',2))) if cds else []
    if not features:status='failed_parse';note='Successful marker existed but no parseable features were found'
    elif not cds:status='failed_no_cds';note='Successful marker existed but no CDS were found'
    elif not rows:status='failed_no_coding_rows';note='Successful marker existed but CDS produced no coding rows'
    else:status='completed'
   else:
    for cmd in templates(exe,fasta,raw,settings):
     attempted.append(cmd);x=subprocess.run(['bash','-lc',activate(settings)+' && '+cmd],text=True,capture_output=True);Path(logs['stdout']).write_text((Path(logs['stdout']).read_text() if Path(logs['stdout']).exists() else '')+x.stdout);Path(logs['stderr']).write_text((Path(logs['stderr']).read_text() if Path(logs['stderr']).exists() else '')+x.stderr);rc=str(x.returncode)
     if x.returncode==0:break
    Path(logs['command']).write_text('\n'.join(attempted)+'\n');Path(logs['returncode']).write_text(rc+'\n')
    if rc != '0':raise RuntimeError(f'runmitos failed. Check mitos2_annotation.settings.refseqver and refdir. See raw/{ref["reference_key"]}/mitos2.stderr.txt')
    features,diag=parse_outputs(raw,ref);write(raw/'parsed_output_files.tsv',DIAG_FIELDS,diag)
    cds=[f for f in features if f['feature_type']=='CDS']; rows=codons(features,fasta,ref,linked,str(settings.get('genetic_code',2))) if cds else []
    if not features:status='failed_parse';note='MITOS2 ran but no parseable features were found'
    elif not cds:status='failed_no_cds';note='MITOS2 output had features but no CDS'
    elif not rows:status='failed_no_coding_rows';note='MITOS2 CDS features produced no coding rows'
    else:status='completed';marker.write_text('completed\n')
   if a.dry_run:rows=[]
  except Exception as e:features=[];diag=[];rows=[];note=str(e);rc=rc or 'exception';Path(logs['returncode']).write_text(rc+'\n');Path(logs['stderr']).write_text((Path(logs['stderr']).read_text() if Path(logs['stderr']).exists() else '')+note+'\n')
  allf+=features;allc+=rows;summ.append({**ref,'status':status,'command_mode':mode,'mitos2_command':exe,'attempted_commands':' | '.join(attempted),'return_code':rc,'stdout_log':logs['stdout'],'stderr_log':logs['stderr'],'help_log':logs['help'],'raw_dir':str(raw),'n_features':len(features),'n_cds_features':len([f for f in features if f['feature_type']=='CDS']),'n_coding_position_rows':len(rows),'n_output_files_scanned':len(diag),'n_parseable_files':sum(bool(d['n_features_parsed']) for d in diag),'parser_status':status,'note':(note+'; stderr_log='+logs['stderr']).strip('; ')})
 write(paths['mitos2_feature_table'],FEATURE_FIELDS,allf);write(paths['mitos2_cds_table'],CODON_FIELDS,allc);write(paths['mitos2_summary_table'],SUMMARY_FIELDS,summ)
 failed=[x for x in summ if x['status'] not in ('completed','dry_run')]
 if failed:
  for x in failed:print(f"ERROR: MITOS2 failed for {x['reference_key']}: {x['note']}",file=sys.stderr)
  raise SystemExit(f'One or more MITOS2 references failed; see {paths["mitos2_summary_table"]}.')
 print(f'Wrote {len(allf)} features and {len(allc)} sample-level coding rows.')
if __name__=='__main__':main()
