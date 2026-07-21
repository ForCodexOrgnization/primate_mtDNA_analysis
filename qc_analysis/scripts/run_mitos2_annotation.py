#!/usr/bin/env python3
"""Run MITOS2 per target-species variant-calling chrM FASTA."""
import argparse, csv, re, shlex, shutil, subprocess, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from qc_analysis.lib.match_utils import yaml
try:
 from Bio import SeqIO
except ImportError: SeqIO = None
REFERENCE_METADATA_FIELDS='target_species final_chrM_species final_chrM_accession coordinate_reference_fasta_from_manifest mitos2_input_fasta'.split()
FEATURE_FIELDS=('reference_key reference_species coordinate_reference_accession coordinate_reference_fasta '+ ' '.join(REFERENCE_METADATA_FIELDS) +' feature_type gene gene_raw start end strand score source_file annotation_source').split()
CODON_FIELDS=('file_name seq_name sample species species_key accession accession_version reference_id family pos ref_base_genome gene gene_raw product protein_id strand codon_index codon_pos_in_triplet codon_seq codon_pos1_genomic codon_pos2_genomic codon_pos3_genomic codon_start_qualifier transl_table cds_tail_incomplete_bases annotation_source coordinate_reference_fasta coordinate_reference_accession '+ ' '.join(REFERENCE_METADATA_FIELDS)).split()
DEBUG_FIELDS='gff_seqid fasta_record_id fasta_length cds_length usable_cds_length n_codons n_position_rows status error gene gene_raw start end strand'.split()
TASK_FIELDS=('task_id reference_key reference_species coordinate_reference_accession coordinate_reference_fasta '+ ' '.join(REFERENCE_METADATA_FIELDS) +' n_samples_using_reference status').split()
SUMMARY_FIELDS=('reference_key reference_species coordinate_reference_accession coordinate_reference_fasta '+ ' '.join(REFERENCE_METADATA_FIELDS) +' status command_mode mitos2_command attempted_commands return_code stdout_log stderr_log help_log raw_dir n_features n_cds_features n_linked_samples n_reference_coding_position_rows n_sample_level_coding_position_rows n_coding_position_rows n_output_files_scanned n_parseable_files result_gff_exists n_gff_gene_rows n_gff_cds_like_gene_rows n_gff_trna_rows n_gff_rrna_rows parser_status note').split()
DIAG_FIELDS='reference_key file suffix n_lines n_candidate_feature_lines parser_used n_features_parsed'.split()
GENES = {
 'ND1':'MT-ND1', 'NAD1':'MT-ND1', 'ND2':'MT-ND2', 'NAD2':'MT-ND2',
 'ND3':'MT-ND3', 'NAD3':'MT-ND3', 'ND4':'MT-ND4', 'NAD4':'MT-ND4',
 'ND4L':'MT-ND4L', 'NAD4L':'MT-ND4L', 'ND5':'MT-ND5', 'NAD5':'MT-ND5',
 'ND6':'MT-ND6', 'NAD6':'MT-ND6', 'COX1':'MT-CO1', 'COI':'MT-CO1',
 'COX2':'MT-CO2', 'COII':'MT-CO2', 'COX3':'MT-CO3', 'COIII':'MT-CO3',
 'COB':'MT-CYB', 'CYTB':'MT-CYB', 'ATP6':'MT-ATP6', 'ATP8':'MT-ATP8',
 'RRNS':'MT-RNR1', 'RRNL':'MT-RNR2',
}
CODING = {key for key, value in GENES.items() if value not in ('MT-RNR1', 'MT-RNR2')}
def val(r,k): return (r.get(k) or '').strip()
def sk(s): return re.sub(r'_+','_',re.sub(r'\s+','_',s.lower())).strip('_')
def write(p,fields,rows):
 Path(p).parent.mkdir(parents=True,exist_ok=True)
 with open(p,'w',newline='') as h:
  w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
def read(p):
 p=Path(p)
 if not p.exists():return []
 with p.open(newline='') as h: rows=[x for x in csv.reader(h,delimiter='\t') if any(y.strip() for y in x)]
 if not rows:return []
 headers={'sample','target_species','final_chrM_species','final_chrM_accession','chrM_expected_output_fasta','reference_key','gff_seqid','status'}
 return [dict(zip(rows[0],x)) for x in rows[1:]] if headers.intersection(rows[0]) else [{'sample':x[0].strip(),'species':x[1].strip() if len(x)>1 else ''} for x in rows]
def attrs(s):
 d={}
 for x in s.split(';'):
  if '=' in x: k,v=x.split('=',1); d[k.lower()]=v.strip('"')
  elif ' ' in x: k,v=x.split(' ',1); d[k.lower()]=v.strip(' "')
 return d
def cleanraw(raw):
 raw = (raw or '').strip()
 raw = re.sub(r'^(?:gene|transcript)_', '', raw, flags=re.I)
 return re.sub(r'\([^)]*\)$', '', raw).strip()
def norm(g):
 raw = cleanraw(g)
 key = re.sub(r'[^A-Z0-9]', '', raw.upper())
 if key.startswith('MT'): key = key[2:]
 return GENES.get(key, raw)
def infer(raw, declared=''):
 feature_type = (declared or '').lower()
 name = re.sub(r'[^a-z0-9]', '', cleanraw(raw).lower())
 if feature_type in ('cds', 'trna', 'rrna'):
  return {'cds':'CDS', 'trna':'tRNA', 'rrna':'rRNA'}[feature_type]
 if feature_type == 'gene':
  if name.upper() in CODING: return 'CDS'
  if name.startswith(('trn', 'trna')): return 'tRNA'
  if name in ('rrns', 'rrnl', '12s', '16s') or 'rrna' in name: return 'rRNA'
 return ''
def gff_diagnostics(raw):
 """Summarize the authoritative MITOS2 GFF, including gene-like CDS evidence."""
 p=Path(raw)/'result.gff'; result={'result_gff_exists':p.is_file(),'n_gff_gene_rows':0,
  'n_gff_cds_like_gene_rows':0,'n_gff_trna_rows':0,'n_gff_rrna_rows':0}
 if not result['result_gff_exists'] or not text_file(p): return result
 for line in p.read_text(errors='replace').splitlines():
  if not line or line.startswith('#'): continue
  c=line.split('\t')
  if len(c)<9: continue
  declared=c[2].lower(); at=attrs(c[8])
  rawgene=at.get('name') or at.get('gene') or at.get('gene_id') or at.get('id') or ''
  if declared=='gene':
   result['n_gff_gene_rows']+=1
   if infer(rawgene, declared)=='CDS': result['n_gff_cds_like_gene_rows']+=1
  elif declared=='trna': result['n_gff_trna_rows']+=1
  elif declared=='rrna': result['n_gff_rrna_rows']+=1
 return result
def parser_failure_status(features, gff):
 if gff['n_gff_cds_like_gene_rows'] and not any(f['feature_type']=='CDS' for f in features):
  return 'failed_parser_cds_gene_detection'
 return 'failed_parse' if not features else 'failed_no_cds'
def text_file(p):
 if p.suffix.lower() in ('.fa','.fasta','.fna','.gz','.bam','.png','.pdf'):return False
 try: p.read_text(errors='strict');return True
 except (UnicodeDecodeError,OSError):return False
def parse_file(p, ref):
 parsed=[]; diagnostics=[]; seen=set(); lines=p.read_text(errors='replace').splitlines(); cand=0; parser='none'
 for line in lines:
  if not line or line.startswith('#'): continue
  c=line.split('\t'); ft=rawgene=''; start=end=strand=score=''
  if len(c) >= 9 and c[3].isdigit() and c[4].isdigit():
   parser='gff'; declared=c[2].lower()
   # MITOS2 GFF: genes represent protein coding intervals; transcript rows represent RNA intervals.
   if declared in ('region', 'exon', 'ncrna_gene'): continue
   at=attrs(c[8]); rawgene=at.get('name') or at.get('gene') or at.get('gene_id') or at.get('id') or ''
   ft=infer(rawgene, declared); start,end,strand,score=c[3],c[4],c[6] or '+',c[5]; cand += bool(ft)
  elif len(c) >= 3 and c[1].isdigit() and c[2].isdigit():
   parser='tabular'; rawgene=c[3] if len(c)>3 else ''; ft=infer(rawgene, rawgene)
   start=str(int(c[1])+1) if p.suffix.lower()=='.bed' else c[1]; end=c[2]; strand=c[5] if len(c)>5 else '+'; cand += bool(ft)
  if not ft: continue
  rawgene=cleanraw(rawgene); key=(ft,start,end,strand,rawgene)
  if key in seen: continue
  seen.add(key)
  parsed.append({**ref,'gff_seqid':c[0] if len(c) >= 9 else '', 'feature_type':ft,'gene':norm(rawgene) if ft in ('CDS', 'rRNA') else rawgene,
                 'gene_raw':rawgene,'start':start,'end':end,'strand':strand,'score':score,
                 'source_file':str(p),'annotation_source':'MITOS2'})
 diagnostics.append({'reference_key':ref['reference_key'],'file':str(p),'suffix':p.suffix,'n_lines':len(lines),
                     'n_candidate_feature_lines':cand,'parser_used':parser,'n_features_parsed':len(parsed)})
 return parsed, diagnostics
def parse_outputs(raw,ref):
 raw=Path(raw); explicit=[raw/'result.gff', raw/'result.bed', raw/'result.mitos']; diagnostics=[]
 # Prefer the authoritative GFF and do not duplicate it with BED/MITOS output.
 for p in explicit:
  if p.is_file() and text_file(p):
   features, diag=parse_file(p, ref); diagnostics += diag
   if features: return features, diagnostics
 allowed={'.gff','.gff3','.bed','.tbl','.tsv','.txt','.out','.result','.mitos',''}
 explicit_set=set(explicit)
 for p in sorted(raw.rglob('*')):
  if p in explicit_set or not p.is_file() or p.name.startswith('mitos2.') or p.suffix.lower() not in allowed or not text_file(p): continue
  features, diag=parse_file(p, ref); diagnostics += diag
  if features: return features, diagnostics
 return [], diagnostics
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
 common=f'-c {q(code)} -o {q(out)} -r {q(refseqver)}'
 if refdir: common += f' -R {q(refdir)}'
 common += ' --best --noplots'
 return [f'{q(exe)} -i {q(fasta)} {common}', f'{q(exe)} --input {q(fasta)} {common}']
def build_reference_codon_rows(features,fasta,ref,code):
 """Build coding-position rows once per reference, before sample expansion."""
 if SeqIO is None: raise RuntimeError('Biopython is required to create MITOS2 codon rows.')
 records=list(SeqIO.parse(str(fasta),'fasta'))
 if not records: raise RuntimeError(f'No FASTA records in {fasta}')
 base=[]; debug=[]
 for f in features:
  if f['feature_type']!='CDS': continue
  gff_seqid=f.get('gff_seqid',''); rec=next((r for r in records if r.id == gff_seqid), records[0] if len(records)==1 else None)
  d={'gff_seqid':gff_seqid,'fasta_record_id':rec.id if rec else '','fasta_length':len(rec) if rec else '','cds_length':'','usable_cds_length':'','n_codons':0,'n_position_rows':0,'status':'','error':'','gene':f['gene'],'gene_raw':f['gene_raw'],'start':f['start'],'end':f['end'],'strand':f['strand']}
  try:
   if rec is None: raise ValueError(f'No FASTA record matches GFF seqid {gff_seqid!r}')
   seq=str(rec.seq).upper(); start,end=int(f['start']),int(f['end'])
   if start < 1 or end < start or end > len(seq): raise ValueError(f'CDS coordinates {start}..{end} outside FASTA length {len(seq)}')
   coords=list(range(start-1,end));strand=f['strand'];dna=''.join(seq[i] for i in coords)
   if strand=='-': coords.reverse();dna=dna.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
   usable=len(dna)//3*3; d.update(cds_length=len(dna),usable_cds_length=usable,n_codons=usable//3)
   for i in range(0,usable,3):
    trip=coords[i:i+3]
    for phase,pos in enumerate(trip,1): base.append({'file_name':Path(fasta).name,'seq_name':rec.id,'sample':'','species':'','species_key':'','accession':ref['coordinate_reference_accession'],'accession_version':ref['coordinate_reference_accession'],'reference_id':ref['coordinate_reference_accession'],'family':'','pos':pos+1,'ref_base_genome':seq[pos],'gene':f['gene'],'gene_raw':f['gene_raw'],'product':f['gene_raw'],'protein_id':'','strand':strand,'codon_index':i//3+1,'codon_pos_in_triplet':phase,'codon_seq':dna[i:i+3],'codon_pos1_genomic':trip[0]+1,'codon_pos2_genomic':trip[1]+1,'codon_pos3_genomic':trip[2]+1,'codon_start_qualifier':'1','transl_table':code,'cds_tail_incomplete_bases':len(dna)-usable,'annotation_source':'MITOS2','coordinate_reference_fasta':str(fasta),'coordinate_reference_accession':ref['coordinate_reference_accession'],**{k:ref.get(k,'') for k in REFERENCE_METADATA_FIELDS}})
   d.update(n_position_rows=usable,status='completed')
  except Exception as exc: d.update(status='failed',error=f'{type(exc).__name__}: {exc}')
  debug.append(d)
 raw_dir=ref.get('raw_dir','')
 if raw_dir: write(Path(raw_dir)/'mitos2_reference_codon_debug.tsv',DEBUG_FIELDS,debug)
 return base
def expand_reference_codon_rows(reference_codon_rows,linked_samples):
 return [{**r,'sample':sample['sample'],'species':sample['species'],'species_key':sk(sample['species'])} for sample in linked_samples for r in reference_codon_rows]
def collect_reference(ref,linked,paths,settings):
 """Return a complete result object for one reference; never abort a batch."""
 raw=Path(paths['mitos2_raw_dir'])/ref['reference_key']; ref={**ref,'raw_dir':str(raw)}; logs={x:str(raw/f'mitos2.{x}.txt') for x in ('command','stdout','stderr','returncode','help')}
 gff=gff_diagnostics(raw); features=[];diag=[];reference_codon_rows=[];sample_codon_rows=[];note=''
 recorded_status=(raw/'mitos2.status.txt').read_text().strip() if (raw/'mitos2.status.txt').exists() else ''
 marker=raw/'mitos2.completed.ok'
 if marker.exists() or recorded_status or gff['result_gff_exists']:
  features,diag=parse_outputs(raw,ref)
  if features: write(raw/'parsed_output_files.tsv',DIAG_FIELDS,diag)
 n_cds=sum(f['feature_type']=='CDS' for f in features)
 if n_cds:
  try: reference_codon_rows=build_reference_codon_rows(features,ref['coordinate_reference_fasta'],ref,str(settings.get('genetic_code',2)))
  except Exception as exc: note=f'Reference codon construction failed: {exc}'
 if reference_codon_rows:
  try: sample_codon_rows=expand_reference_codon_rows(reference_codon_rows,linked)
  except Exception as exc: note=f'Sample codon expansion failed: {exc}'
 debug_rows=read(raw/'mitos2_reference_codon_debug.tsv')
 gene_warnings=any(val(row,'status') != 'completed' for row in debug_rows)
 if not n_cds:
  status=parser_failure_status(features,gff) if (marker.exists() or recorded_status or gff['result_gff_exists']) else 'pending'
 elif not reference_codon_rows: status='failed_reference_codon_construction'
 elif not linked: status='completed_reference_no_linked_samples'
 elif sample_codon_rows: status='completed_with_gene_warnings' if gene_warnings else 'completed'
 else: status='failed_sample_expansion'
 if recorded_status and not gff['result_gff_exists'] and not features: status=recorded_status
 rc=Path(logs['returncode']).read_text().strip() if Path(logs['returncode']).exists() else ''
 command_text=Path(logs['command']).read_text().strip() if Path(logs['command']).exists() else ''
 summary={**ref,'status':status,'command_mode':'runmitos','mitos2_command':'runmitos','attempted_commands':command_text,'return_code':rc,'stdout_log':logs['stdout'],'stderr_log':logs['stderr'],'help_log':logs['help'],'raw_dir':str(raw),'n_features':len(features),'n_cds_features':n_cds,'n_linked_samples':len(linked),'n_reference_coding_position_rows':len(reference_codon_rows),'n_sample_level_coding_position_rows':len(sample_codon_rows),'n_coding_position_rows':len(sample_codon_rows),'n_output_files_scanned':len(diag),'n_parseable_files':sum(bool(d['n_features_parsed']) for d in diag),**gff,'parser_status':status,'note':note}
 return {'features':features,'reference_codon_rows':reference_codon_rows,'sample_codon_rows':sample_codon_rows,'summary_row':summary,'status':status,'note':note}
def sanitized_fallback_fasta(source, target, paths):
 """Copy a manifest fallback with MITOS2's stable ``>chrM`` record identifier."""
 source=Path(source); destination=Path(paths['mitos2_raw_dir'])/'input_fastas'/(re.sub(r'[^A-Za-z0-9_.-]+','_',target)+'.fa')
 try:
  lines=source.read_text().splitlines()
  if not lines or not lines[0].startswith('>'): raise ValueError('FASTA has no header')
  destination.parent.mkdir(parents=True,exist_ok=True)
  destination.write_text('>chrM\n'+'\n'.join(line.strip() for line in lines[1:] if line.strip())+'\n')
  return str(destination)
 except (OSError, ValueError): return str(source)
def references(paths, sample_filter=None):
 manifest=read(paths['reference_manifest']); samples=read(paths['sample_ref_file'])
 if sample_filter: samples=[s for s in samples if val(s,'sample')==sample_filter]
 refs={}
 for m in manifest:
  target=val(m,'target_species'); species=val(m,'final_chrM_species') or target
  fasta_dir=paths.get('final_chrM_fasta_dir',paths.get('fasta_dir','references/variant_calling/Ref_chrM'))
  standardized=Path(fasta_dir)/(target+'.fa')
  manifest_fasta=val(m,'chrM_expected_output_fasta')
  no_chrm=any(val(m,k) in ('wg_only_no_chrM','missing_chrM_ref') for k in ('final_reference_strategy','chrM_reference_context','status'))
  if standardized.is_file(): fasta=str(standardized); status='pending'
  elif manifest_fasta: fasta=sanitized_fallback_fasta(manifest_fasta,target,paths); status='pending'
  else: fasta=str(standardized); status='skipped_no_chrM_reference' if no_chrm else 'pending'
  acc=val(m,'final_chrM_accession') or val(m,'final_chrM_refseq_accn') or val(m,'final_chrM_genbank_accn')
  # A target species is the task identity: cross-species targets can share an accession but not a FASTA.
  key=re.sub(r'[^A-Za-z0-9_.-]+','_',target or Path(fasta).stem)
  refs[key]={'reference_key':key,'reference_species':species,'coordinate_reference_accession':acc,
             'coordinate_reference_fasta':str(Path(fasta)),'coordinate_reference_fasta_from_manifest':manifest_fasta,
             'mitos2_input_fasta':str(Path(fasta)),'target_species':target,'final_chrM_species':val(m,'final_chrM_species'),
             'final_chrM_accession':val(m,'final_chrM_accession'),'targets':{sk(target)},'initial_status':status}
 result=[]
 for ref in refs.values():
  linked=[{'sample':val(s,'sample'),'species':val(s,'species')} for s in samples if sk(val(s,'species')) in ref['targets'] or sk(val(s,'species'))==sk(ref['reference_species'])]
  result.append((ref,linked))
 return sorted(result,key=lambda pair:(pair[0]['reference_key'], pair[0]['coordinate_reference_fasta']))
def task_rows(refs, paths):
 rows=[]
 for task_id,(ref,linked) in enumerate(refs,1):
  marker=Path(paths['mitos2_raw_dir'])/ref['reference_key']/'mitos2.completed.ok'
  rows.append({'task_id':task_id,**{k:ref[k] for k in TASK_FIELDS if k in ref},'n_samples_using_reference':len(linked),'status':'completed' if marker.exists() else ref.get('initial_status','pending')})
 return rows
def merge(paths,settings,refs):
 results=[collect_reference(ref,linked,paths,settings) for ref,linked in refs]
 allf=[row for result in results for row in result['features']]
 allc=[row for result in results for row in result['sample_codon_rows']]
 allref=[row for result in results for row in result['reference_codon_rows']]
 summ=[result['summary_row'] for result in results]
 reference_table=paths.get('mitos2_reference_cds_table',str(Path(paths['output_dir'])/'all_mitos2_reference_position_codon_table.tsv'))
 write(paths['mitos2_feature_table'],FEATURE_FIELDS,allf);write(paths['mitos2_cds_table'],CODON_FIELDS,allc);write(reference_table,CODON_FIELDS,allref);write(paths['mitos2_summary_table'],SUMMARY_FIELDS,summ)
 print(f'Wrote {len(allf)} features and {len(allc)} sample-level coding rows.')
def run_reference(ref,linked,paths,settings,a):
 """Execute one MITOS2 reference and always return its materialized result."""
 fasta=ref['mitos2_input_fasta'];raw=Path(paths['mitos2_raw_dir'])/ref['reference_key'];raw.mkdir(parents=True,exist_ok=True)
 logs={x:str(raw/f'mitos2.{x}.txt') for x in ('command','stdout','stderr','returncode','help')}; marker=raw/'mitos2.completed.ok';status_file=raw/'mitos2.status.txt'
 if ref.get('initial_status') == 'skipped_no_chrM_reference':
  status_file.write_text('skipped_no_chrM_reference\n'); return collect_reference(ref,linked,paths,settings)
 if marker.exists() and not a.force:
  print(f'Skipping completed MITOS2 reference: {ref["reference_key"]}')
  return collect_reference(ref,linked,paths,settings)
 if a.force:
  for p in raw.glob('result.*'): p.unlink(missing_ok=True)
  for name in ('ignored.mitos','stst.dat','mitos2.completed.ok','mitos2.status.txt','parsed_output_files.tsv','mitos2_reference_codon_debug.tsv'):
   (raw/name).unlink(missing_ok=True)
  for name in ('blast','mitfi-global'):
   shutil.rmtree(raw/name,ignore_errors=True)
 for path in logs.values(): Path(path).write_text('')
 if not Path(fasta).exists():
  Path(logs['returncode']).write_text('exception\n');Path(logs['stderr']).write_text(f'MITOS2 input FASTA is missing: {fasta}\n');status_file.write_text('failed_missing_fasta\n');return collect_reference(ref,linked,paths,settings)
 if a.dry_run:
  Path(logs['command']).write_text('dry-run\n');Path(logs['returncode']).write_text('0\n');status_file.write_text('dry_run\n');return collect_reference(ref,linked,paths,settings)
 attempted=[];rc='';success=False;failure_status='failed_mitos2_execution'
 try:
  exe,validation=command(settings);helpx=subprocess.run(['bash','-lc',activate(settings)+f' && {shlex.quote(exe)} --help'],text=True,capture_output=True);Path(logs['help']).write_text(validation+helpx.stdout+'\n'+helpx.stderr)
  for cmd in templates(exe,fasta,raw,settings):
   attempted.append(cmd);x=subprocess.run(['bash','-lc',activate(settings)+' && '+cmd],text=True,capture_output=True);Path(logs['stdout']).write_text(Path(logs['stdout']).read_text()+x.stdout);Path(logs['stderr']).write_text(Path(logs['stderr']).read_text()+x.stderr);rc=str(x.returncode)
   if x.returncode: continue
   if not all((raw/name).is_file() for name in ('result.gff','result.bed','result.mitos')):
    failure_status='failed_mitos2_execution'; attempted.append('runmitos_missing_required_raw_output'); continue
   result=collect_reference(ref,linked,paths,settings)
   if result['status'] in ('completed','completed_with_gene_warnings','completed_reference_no_linked_samples'): success=True;break
   failure_status=result['status'];attempted.append('template_returned_zero_but_invalid_materialized_output')
  Path(logs['command']).write_text('\n'.join(attempted)+'\n');Path(logs['returncode']).write_text(rc+'\n')
 except Exception as exc:
  Path(logs['returncode']).write_text((rc or 'exception')+'\n');Path(logs['stderr']).write_text(Path(logs['stderr']).read_text()+str(exc)+'\n');failure_status='failed_mitos2_execution'
 if success: status_file.unlink(missing_ok=True);marker.write_text('completed\n')
 else:
  status_file.write_text(failure_status+'\n');print(f'MITOS2 reference failed: {ref["reference_key"]} ({failure_status}); continuing.')
 return collect_reference(ref,linked,paths,settings)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--config',required=True);ap.add_argument('--sample');ap.add_argument('--prepare-tasks',action='store_true');ap.add_argument('--task-id');ap.add_argument('--reference');ap.add_argument('--merge-only',action='store_true');ap.add_argument('--force',action='store_true');ap.add_argument('--dry-run',action='store_true');a=ap.parse_args()
 if sum(bool(x) for x in (a.prepare_tasks,a.merge_only)) and (a.task_id or a.reference): ap.error('--prepare-tasks/--merge-only cannot be combined with --task-id or --reference')
 sec=yaml(a.config).get('mitos2_annotation');
 if not sec:raise SystemExit('Missing mitos2_annotation section in config.')
 paths,settings=sec['paths'],sec.get('settings',{}); refs=references(paths,a.sample)
 task_path=paths.get('mitos2_reference_tasks',str(Path(paths['output_dir'])/'mitos2_reference_tasks.tsv'))
 if a.prepare_tasks:
  write(task_path,TASK_FIELDS,task_rows(refs,paths))
  print(f'Wrote {len(refs)} MITOS2 reference tasks to {task_path}.');return
 if a.merge_only:
  merge(paths,settings,refs);write(task_path,TASK_FIELDS,task_rows(refs,paths));return
 if a.task_id:
  selected=[pair for task,pair in zip(task_rows(refs,paths),refs) if str(task['task_id'])==str(a.task_id)]
  if not selected: raise SystemExit(f'No MITOS2 task found with task_id {a.task_id}.')
 elif a.reference: selected=[pair for pair in refs if a.reference in (pair[0]['reference_key'],pair[0]['reference_species'],pair[0]['coordinate_reference_accession'])]
 else: selected=refs
 if not selected: raise SystemExit('No MITOS2 references selected.')
 for ref,linked in selected:
  run_reference(ref,linked,paths,settings,a)
 # Array workers must not concurrently rewrite combined output tables.
 if not (a.task_id or a.reference):
  merge(paths,settings,refs);write(task_path,TASK_FIELDS,task_rows(refs,paths))
if __name__=='__main__':main()
