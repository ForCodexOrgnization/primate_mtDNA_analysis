#!/usr/bin/env python3
"""Run MITOS2 per target-species variant-calling chrM FASTA."""
import argparse, csv, hashlib, re, shlex, shutil, subprocess, sys
from contextlib import ExitStack
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from qc_analysis.lib.match_utils import yaml
from qc_analysis.lib.reference_utils import normalized_fasta_sequence_sha256
try:
 from Bio import SeqIO
except ImportError: SeqIO = None
REFERENCE_METADATA_FIELDS='target_species final_chrM_species final_chrM_accession coordinate_reference_fasta_from_manifest mitos2_input_fasta'.split()
FEATURE_FIELDS=('reference_key reference_species coordinate_reference_accession coordinate_reference_fasta '+ ' '.join(REFERENCE_METADATA_FIELDS) +' feature_type gene gene_raw start end strand score source_file annotation_source').split()
CODON_FIELDS=('file_name seq_name sample species species_key accession accession_version reference_id family pos ref_base_genome gene gene_raw product protein_id strand codon_index codon_pos_in_triplet codon_seq codon_pos1_genomic codon_pos2_genomic codon_pos3_genomic codon_start_qualifier transl_table cds_tail_incomplete_bases annotation_source coordinate_reference_fasta coordinate_reference_accession coordinate_reference_sequence_sha256 mitos2_input_sequence_sha256 mitos2_input_sequence_length '+ ' '.join(REFERENCE_METADATA_FIELDS)).split()
REFERENCE_CODON_FIELDS=('reference_key reference_species coordinate_reference_fasta coordinate_reference_accession coordinate_reference_sequence_sha256 mitos2_input_sequence_sha256 pos ref_base_genome gene gene_raw strand codon_index codon_pos_in_triplet codon_seq codon_pos1_genomic codon_pos2_genomic codon_pos3_genomic transl_table annotation_source file_name seq_name accession accession_version reference_id product protein_id codon_start_qualifier cds_tail_incomplete_bases mitos2_input_sequence_length '+ ' '.join(REFERENCE_METADATA_FIELDS)).split()
SAMPLE_REFERENCE_FIELDS='sample species species_key reference_key coordinate_reference_fasta coordinate_reference_accession coordinate_reference_sequence_sha256'.split()
DEBUG_FIELDS='gff_seqid fasta_record_id fasta_length sequence_length original_gff_start original_gff_end canonical_start canonical_end circular_wrap_used wrapped_segment_count cds_length usable_cds_length n_codons n_position_rows status error gene gene_raw start end strand'.split()
TASK_FIELDS=('task_id task_key reference_key reference_species coordinate_reference_accession coordinate_reference_fasta coordinate_reference_sequence_sha256 mitos2_input_sequence_sha256 mitos2_input_sequence_length '+ ' '.join(REFERENCE_METADATA_FIELDS) +' n_samples_using_reference status').split()
SUMMARY_FIELDS=('task_key reference_key reference_species coordinate_reference_accession coordinate_reference_fasta coordinate_reference_sequence_sha256 mitos2_input_sequence_sha256 mitos2_input_sequence_length '+ ' '.join(REFERENCE_METADATA_FIELDS) +' status production_qc_status production_qc_reasons command_mode mitos2_command attempted_commands return_code stdout_log stderr_log help_log raw_dir n_features n_cds_features n_linked_samples n_reference_coding_position_rows n_sample_level_coding_position_rows n_coding_position_rows n_output_files_scanned n_parseable_files result_gff_exists n_gff_gene_rows n_gff_cds_like_gene_rows n_gff_trna_rows n_gff_rrna_rows parser_status note').split()
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
EXPECTED_GENES={f'MT-ND{i}' for i in (1,2,3,4,5,6)}|{'MT-ND4L','MT-CO1','MT-CO2','MT-CO3','MT-CYB','MT-ATP6','MT-ATP8'}
IUPAC=set('ACGTRYSWKMBDHVN')
def val(r,k): return (r.get(k) or '').strip()
def sk(s): return re.sub(r'_+','_',re.sub(r'\s+','_',s.lower())).strip('_')
def biological_reference_key(sequence_sha256):
 if not re.fullmatch(r'[0-9a-f]{64}', sequence_sha256 or ''): return ''
 return 'mtref_'+sequence_sha256
def task_key(sequence_sha256): return 'seq_'+sequence_sha256
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
def circular_interval_coordinates(start: int, end: int, sequence_length: int) -> list[int]:
 """Return ordered zero-based coordinates for a circular GFF interval.

 GFF coordinates are 1-based and inclusive; returned coordinates are zero-based
 and retain traversal order across the circular genome origin.
 """
 if sequence_length <= 0:
  raise ValueError(f'sequence_length must be positive, observed {sequence_length}')
 if start < 1:
  raise ValueError(f'start must be >= 1, observed {start}')
 if end < start:
  raise ValueError(f'end must be >= start, observed start={start}, end={end}')
 feature_length=end-start+1
 if feature_length > sequence_length:
  raise ValueError('Circular feature length exceeds one full circular genome: '
                   f'feature_length={feature_length}, sequence_length={sequence_length}')
 return [(extended_pos-1) % sequence_length for extended_pos in range(start,end+1)]
def build_reference_codon_rows(features,fasta,ref,code):
 """Build coding-position rows once per reference, before sample expansion."""
 if SeqIO is None: raise RuntimeError('Biopython is required to create MITOS2 codon rows.')
 records=list(SeqIO.parse(str(fasta),'fasta'))
 if not records: raise RuntimeError(f'No FASTA records in {fasta}')
 base=[]; debug=[]
 for f in features:
  if f['feature_type']!='CDS': continue
  gff_seqid=f.get('gff_seqid',''); rec=next((r for r in records if r.id == gff_seqid), records[0] if len(records)==1 else None)
  d={'gff_seqid':gff_seqid,'fasta_record_id':rec.id if rec else '','fasta_length':len(rec) if rec else '',
     'sequence_length':len(rec) if rec else '','original_gff_start':f['start'],'original_gff_end':f['end'],
     'canonical_start':'','canonical_end':'','circular_wrap_used':'','wrapped_segment_count':'',
     'cds_length':'','usable_cds_length':'','n_codons':0,'n_position_rows':0,'status':'','error':'',
     'gene':f['gene'],'gene_raw':f['gene_raw'],'start':f['start'],'end':f['end'],'strand':f['strand']}
  try:
   if rec is None: raise ValueError(f'No FASTA record matches GFF seqid {gff_seqid!r}')
   seq=str(rec.seq).upper(); start,end=int(f['start']),int(f['end'])
   sequence_length=len(seq); cds_length=end-start+1
   d.update(sequence_length=sequence_length,original_gff_start=start,original_gff_end=end,
            canonical_start=(start-1) % sequence_length+1 if sequence_length else '',
            canonical_end=(end-1) % sequence_length+1 if sequence_length else '',
            circular_wrap_used='yes' if end > sequence_length else 'no',
            wrapped_segment_count=2 if end > sequence_length else 1)
   try: coords=circular_interval_coordinates(start,end,sequence_length)
   except ValueError as exc:
    raise ValueError(f"Reference {ref.get('reference_key','')}, gene {f['gene']}: {exc}") from exc
   strand=f['strand'];genomic_dna=''.join(seq[position] for position in coords)
   if strand=='-': coords=list(reversed(coords));dna=genomic_dna.translate(str.maketrans('ACGTRYSWKMBDHVN','TGCAYRSWMKVHDBN'))[::-1]
   else: dna=genomic_dna
   usable=len(dna)//3*3; d.update(cds_length=len(dna),usable_cds_length=usable,n_codons=usable//3)
   for i in range(0,usable,3):
    trip=coords[i:i+3]
    for phase,pos in enumerate(trip,1): base.append({'reference_key':ref['reference_key'],'reference_species':ref.get('reference_species',''),'file_name':Path(fasta).name,'seq_name':rec.id,'sample':'','species':'','species_key':'','accession':ref['coordinate_reference_accession'],'accession_version':ref['coordinate_reference_accession'],'reference_id':ref['coordinate_reference_accession'],'family':'','pos':pos+1,'ref_base_genome':seq[pos],'gene':f['gene'],'gene_raw':f['gene_raw'],'product':f['gene_raw'],'protein_id':'','strand':strand,'codon_index':i//3+1,'codon_pos_in_triplet':phase,'codon_seq':dna[i:i+3],'codon_pos1_genomic':trip[0]+1,'codon_pos2_genomic':trip[1]+1,'codon_pos3_genomic':trip[2]+1,'codon_start_qualifier':'1','transl_table':code,'cds_tail_incomplete_bases':len(dna)-usable,'annotation_source':'MITOS2','coordinate_reference_fasta':str(fasta),'coordinate_reference_accession':ref['coordinate_reference_accession'],'coordinate_reference_sequence_sha256':ref.get('coordinate_reference_sequence_sha256',''),'mitos2_input_sequence_sha256':ref.get('mitos2_input_sequence_sha256',''),'mitos2_input_sequence_length':ref.get('mitos2_input_sequence_length',''),**{k:ref.get(k,'') for k in REFERENCE_METADATA_FIELDS}})
   d.update(n_position_rows=usable,status='completed')
  except Exception as exc: d.update(status='failed',error=f'{type(exc).__name__}: {exc}')
  debug.append(d)
 raw_dir=ref.get('raw_dir','')
 if raw_dir: write(Path(raw_dir)/'mitos2_reference_codon_debug.tsv',DEBUG_FIELDS,debug)
 return base
def production_qc(ref, rows, settings):
 """Return strict production disposition and all reasons; never hide a failure."""
 failures=[]; warnings=[]; fasta=Path(ref.get('coordinate_reference_fasta',''))
 try:
  if not fasta.is_file(): raise OSError('not a readable file')
  actual=normalized_fasta_sequence_sha256(fasta)
  if actual['sequence_sha256'] != ref.get('coordinate_reference_sequence_sha256'): failures.append('coordinate_fasta_hash_changed')
 except (OSError,ValueError) as exc: failures.append(f'coordinate_fasta_unreadable:{exc}')
 if not ref.get('coordinate_reference_sequence_sha256') or ref.get('mitos2_input_sequence_sha256') != ref.get('coordinate_reference_sequence_sha256'): failures.append('mitos2_input_hash_mismatch')
 genes={r.get('gene') for r in rows}
 if genes != EXPECTED_GENES: failures.append('expected_13_genes_missing_or_extra:'+','.join(sorted(EXPECTED_GENES-genes)))
 if any(r.get('strand') not in ('+','-') for r in rows): failures.append('invalid_strand')
 if any(str(r.get('codon_pos_in_triplet')) not in ('1','2','3') for r in rows): failures.append('invalid_codon_position')
 if any(len(r.get('codon_seq','')) != 3 or not set(r.get('codon_seq','')).issubset(IUPAC) for r in rows): failures.append('invalid_iupac_codon')
 try: seq=''.join(str(rec.seq).upper() for rec in SeqIO.parse(str(fasta),'fasta'))
 except Exception: seq=''
 if any(not str(r.get('pos','')).isdigit() or int(r['pos'])<1 or int(r['pos'])>len(seq) or r.get('ref_base_genome') != seq[int(r['pos'])-1] for r in rows): failures.append('reference_base_disagrees_with_fasta')
 resolved={}
 for r in rows: resolved.setdefault((r.get('reference_key'),r.get('pos')),set()).add(r.get('ref_base_genome'))
 if any(len(v)>1 for v in resolved.values()): failures.append('conflicting_reference_bases')
 low=int(settings.get('min_production_coding_position_rows',9000)); high=int(settings.get('max_production_coding_position_rows',13000))
 if not low <= len(rows) <= high: failures.append(f'implausible_coding_position_rows:{len(rows)}')
 if not rows: failures.append('codon_construction_failed')
 return ('FAIL_PRODUCTION' if failures else 'WARN_PRODUCTION' if warnings else 'PASS_PRODUCTION'), ';'.join(failures+warnings)
def collect_reference(ref,linked,paths,settings):
 """Return a complete result object for one reference; never abort a batch."""
 raw=Path(paths['mitos2_raw_dir'])/ref.get('task_key',ref['reference_key']); ref={**ref,'raw_dir':str(raw)}; logs={x:str(raw/f'mitos2.{x}.txt') for x in ('command','stdout','stderr','returncode','help')}
 if ref.get('status') == 'skipped_no_chrM_reference':
  summary={**ref,'status':'skipped_no_chrM_reference','production_qc_status':'FAIL_PRODUCTION','production_qc_reasons':'coordinate_fasta_unavailable','command_mode':'not_run','mitos2_command':'','attempted_commands':'','return_code':'','stdout_log':'','stderr_log':'','help_log':'','n_features':0,'n_cds_features':0,'n_linked_samples':len(linked),'n_reference_coding_position_rows':0,'n_sample_level_coding_position_rows':0,'n_coding_position_rows':0,'n_output_files_scanned':0,'n_parseable_files':0,'result_gff_exists':False,'n_gff_gene_rows':0,'n_gff_cds_like_gene_rows':0,'n_gff_trna_rows':0,'n_gff_rrna_rows':0,'parser_status':'skipped_no_chrM_reference','note':'Manifest row has no materialized chrM FASTA; accession is not a production fallback.'}
  return {'features':[],'reference_codon_rows':[],'summary_row':summary,'status':'skipped_no_chrM_reference','note':summary['note']}
 gff=gff_diagnostics(raw); features=[];diag=[];reference_codon_rows=[];note=''
 recorded_status=(raw/'mitos2.status.txt').read_text().strip() if (raw/'mitos2.status.txt').exists() else ''
 marker=raw/'mitos2.completed.ok'
 if marker.exists() or recorded_status or gff['result_gff_exists']:
  features,diag=parse_outputs(raw,ref)
  if features: write(raw/'parsed_output_files.tsv',DIAG_FIELDS,diag)
 n_cds=sum(f['feature_type']=='CDS' for f in features)
 if n_cds:
  try: reference_codon_rows=build_reference_codon_rows(features,ref['coordinate_reference_fasta'],ref,str(settings.get('genetic_code',2)))
  except Exception as exc: note=f'Reference codon construction failed: {exc}'
 debug_rows=read(raw/'mitos2_reference_codon_debug.tsv')
 gene_warnings=any(val(row,'status') != 'completed' for row in debug_rows)
 if not n_cds:
  status=parser_failure_status(features,gff) if (marker.exists() or recorded_status or gff['result_gff_exists']) else 'pending'
 elif not reference_codon_rows: status='failed_reference_codon_construction'
 elif not linked: status='completed_reference_no_linked_samples'
 else: status='completed_with_gene_warnings' if gene_warnings else 'completed'
 if recorded_status and not gff['result_gff_exists'] and not features: status=recorded_status
 rc=Path(logs['returncode']).read_text().strip() if Path(logs['returncode']).exists() else ''
 command_text=Path(logs['command']).read_text().strip() if Path(logs['command']).exists() else ''
 qc,qc_reasons=production_qc(ref,reference_codon_rows,settings)
 summary={**ref,'status':status,'production_qc_status':qc,'production_qc_reasons':qc_reasons,'command_mode':'runmitos','mitos2_command':'runmitos','attempted_commands':command_text,'return_code':rc,'stdout_log':logs['stdout'],'stderr_log':logs['stderr'],'help_log':logs['help'],'raw_dir':str(raw),'n_features':len(features),'n_cds_features':n_cds,'n_linked_samples':len(linked),'n_reference_coding_position_rows':len(reference_codon_rows),'n_sample_level_coding_position_rows':0,'n_coding_position_rows':len(reference_codon_rows),'n_output_files_scanned':len(diag),'n_parseable_files':sum(bool(d['n_features_parsed']) for d in diag),**gff,'parser_status':status,'note':note}
 return {'features':features,'reference_codon_rows':reference_codon_rows,'summary_row':summary,'status':status,'note':note}
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
 # This is the authoritative coordinate-assignment index.  In particular,
 # reference_species is provenance and must never participate in this index.
 target_references={}
 for m in manifest:
  target=val(m,'target_species'); species=val(m,'final_chrM_species') or target
  fasta_dir=paths.get('final_chrM_fasta_dir',paths.get('fasta_dir','references/variant_calling/Ref_chrM'))
  standardized=Path(fasta_dir)/(target+'.fa')
  manifest_fasta=val(m,'chrM_expected_output_fasta')
  no_chrm=any(val(m,k) in ('wg_only_no_chrM','missing_chrM_ref') for k in ('final_reference_strategy','chrM_reference_context','status'))
  if no_chrm and val(m,'chrM_selection_status') == 'missing_chrM_ref': fasta=''; status='skipped_no_chrM_reference'
  elif standardized.is_file(): fasta=str(standardized); status='pending'
  elif manifest_fasta: fasta=sanitized_fallback_fasta(manifest_fasta,target,paths); status='pending'
  else: fasta=str(standardized); status='skipped_no_chrM_reference' if no_chrm else 'pending'
  acc=val(m,'final_chrM_accession') or val(m,'final_chrM_refseq_accn') or val(m,'final_chrM_genbank_accn')
  try:
   hash_info=normalized_fasta_sequence_sha256(fasta) if fasta and Path(fasta).is_file() else {}
  except (OSError, ValueError) as exc:
   hash_info={}
  sequence_sha=hash_info.get('sequence_sha256','')
  # Execution is deduplicated by sequence.  Species and accession are provenance,
  # never biological identity and never a reason to collapse distinct sequences.
  key=biological_reference_key(sequence_sha) or 'missing_'+hashlib.sha256(str(Path(fasta).resolve()).encode()).hexdigest()
  tkey=task_key(sequence_sha) if sequence_sha else 'missing_'+re.sub(r'[^A-Za-z0-9_.-]+','_',target or Path(fasta).stem)
  candidate={'task_key':tkey,'reference_key':key,'reference_species':species,'coordinate_reference_accession':acc,
             'coordinate_reference_fasta':str(Path(fasta).resolve()) if fasta else '',
             'coordinate_reference_sequence_sha256':sequence_sha,
             'mitos2_input_sequence_sha256':sequence_sha,
             'mitos2_input_sequence_length':hash_info.get('sequence_length',''),'coordinate_reference_fasta_from_manifest':manifest_fasta,
             'mitos2_input_fasta':str(Path(fasta).resolve()) if fasta else '','target_species':target,'final_chrM_species':val(m,'final_chrM_species'),
             'final_chrM_accession':val(m,'final_chrM_accession'),'chrM_selection_status':val(m,'chrM_selection_status'),
             'final_reference_strategy':val(m,'final_reference_strategy'),'reference_pairing_status':val(m,'reference_pairing_status'),
             'targets':{sk(target)},'target_records':{},'initial_status':status,
             **({'status':status} if val(m,'chrM_selection_status') == 'missing_chrM_ref' else {})}
  candidate['target_records'][sk(target)]={'coordinate_reference_fasta':str(Path(fasta).resolve()) if fasta else '',
   'coordinate_reference_accession':acc,'coordinate_reference_sequence_sha256':sequence_sha}
  target_key=sk(target)
  assignment={'reference_key':key,**candidate['target_records'][target_key]}
  previous=target_references.get(target_key)
  if previous and previous != assignment:
   raise SystemExit('Conflicting coordinate references for target species '
    f'{target!r}: reference_keys={[previous["reference_key"], key]}; '
    f'FASTA paths={[previous["coordinate_reference_fasta"], assignment["coordinate_reference_fasta"]]}; '
    f'accessions={[previous["coordinate_reference_accession"], acc]}; '
    f'SHA256 values={[previous["coordinate_reference_sequence_sha256"], sequence_sha]}')
  target_references[target_key]=assignment
  if key in refs:
   refs[key]['targets'].update(candidate['targets']); refs[key]['target_records'].update(candidate['target_records'])
  else: refs[key]=candidate
 result=[]
 for ref in refs.values():
  linked=[]
  for s in samples:
   sample_species=val(s,'target_species') or val(s,'species')
   species_key=sk(sample_species)
   assignment=target_references.get(species_key)
   if assignment and assignment['reference_key']==ref['reference_key']:
    provenance={k:assignment[k] for k in ('coordinate_reference_fasta','coordinate_reference_accession','coordinate_reference_sequence_sha256')}
    linked.append({'sample':val(s,'sample'),'species':val(s,'species'),**provenance})
  result.append((ref,linked))
 unresolved=[f"{val(s,'sample')} ({val(s,'target_species') or val(s,'species')})" for s in samples
             if sk(val(s,'target_species') or val(s,'species')) not in target_references]
 if unresolved:
  print('Samples with no resolved target coordinate reference: '+', '.join(unresolved),file=sys.stderr)
 return sorted(result,key=lambda pair:(pair[0]['reference_key'], pair[0]['coordinate_reference_fasta']))
def task_rows(refs, paths):
 rows=[]
 for ref,linked in refs:
  if ref.get('status') == 'skipped_no_chrM_reference': continue
  task_id=len(rows)+1
  marker=Path(paths['mitos2_raw_dir'])/ref['task_key']/'mitos2.completed.ok'
  rows.append({'task_id':task_id,**{k:ref[k] for k in TASK_FIELDS if k in ref},'n_samples_using_reference':len(linked),'status':'completed' if marker.exists() else ref.get('initial_status','pending')})
 return rows
def sample_reference_rows(refs):
 """Describe variant-calling coordinates without consulting annotation QC."""
 rows=[{'sample':sample['sample'],'species':sample['species'],'species_key':sk(sample['species']),
          'reference_key':ref['reference_key'],'coordinate_reference_fasta':sample['coordinate_reference_fasta'],
          'coordinate_reference_accession':sample['coordinate_reference_accession'],
          'coordinate_reference_sequence_sha256':sample['coordinate_reference_sequence_sha256']}
         for ref,linked in refs for sample in linked]
 return validate_sample_reference_rows(rows,'sample_coordinate_reference_map.tsv')
def validate_sample_reference_rows(rows, table_name):
 """Collapse exact duplicates and reject distinct coordinate references per sample."""
 unique=[]; seen=set(); by_sample={}
 for row in rows:
  signature=tuple((field,str(row.get(field,''))) for field in SAMPLE_REFERENCE_FIELDS)
  if signature in seen: continue
  seen.add(signature); unique.append(row)
  by_sample.setdefault(row.get('sample',''),[]).append(row)
 conflicts={sample: entries for sample,entries in by_sample.items()
            if len({entry.get('reference_key','') for entry in entries}) > 1}
 if conflicts:
  details=[]
  for sample,entries in sorted(conflicts.items()):
   details.append(f"sample={sample!r}; species={sorted({e.get('species','') for e in entries})}; "
    f"reference_keys={sorted({e.get('reference_key','') for e in entries})}; "
    f"FASTA paths={sorted({e.get('coordinate_reference_fasta','') for e in entries})}; "
    f"accessions={sorted({e.get('coordinate_reference_accession','') for e in entries})}; "
    f"SHA256 values={sorted({e.get('coordinate_reference_sequence_sha256','') for e in entries})}")
  raise SystemExit(f'Conflicting reference keys in {table_name}: '+' | '.join(details))
 return unique
def merge(paths,settings,refs):
 """Merge one reference at a time without constructing sample-expanded rows."""
 reference_table=paths.get('mitos2_reference_cds_table',str(Path(paths['output_dir'])/'all_mitos2_reference_position_codon_table.tsv'))
 mapping_table=paths.get('codon_sample_reference_map',str(Path(paths['output_dir'])/'codon_sample_reference_map.tsv'))
 outputs=[(reference_table,REFERENCE_CODON_FIELDS),(mapping_table,SAMPLE_REFERENCE_FIELDS),
          (paths['mitos2_summary_table'],SUMMARY_FIELDS)]
 if paths.get('mitos2_feature_table'): outputs.append((paths['mitos2_feature_table'],FEATURE_FIELDS))
 # Validate the complete assignment set before opening either output.  The
 # production mapping is a subset of these rows, so it cannot introduce a
 # distinct-key conflict; it is validated again as rows are selected below.
 validated_assignments=validate_sample_reference_rows(sample_reference_rows(refs),Path(mapping_table).name)
 assignment_by_sample={row['sample']:row for row in validated_assignments}
 counts={'references':0,'codons':0,'mappings':0,'features':0}; codon_mapping_rows=[]
 with ExitStack() as stack:
  writers={}
  for path,fields in outputs:
   Path(path).parent.mkdir(parents=True,exist_ok=True)
   handle=stack.enter_context(open(path,'w',newline=''))
   writers[path]=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',extrasaction='ignore');writers[path].writeheader()
  for ref,linked in refs:
   result=collect_reference(ref,linked,paths,settings)
   writers[paths['mitos2_summary_table']].writerow(result['summary_row']);counts['references']+=1
   if paths.get('mitos2_feature_table'):
    writers[paths['mitos2_feature_table']].writerows(result['features']);counts['features']+=len(result['features'])
   if result['summary_row'].get('production_qc_status')=='PASS_PRODUCTION':
    writers[reference_table].writerows(result['reference_codon_rows']);counts['codons']+=len(result['reference_codon_rows'])
    codon_mapping_rows.extend(assignment_by_sample[sample['sample']] for sample in linked)
   del result
  codon_mapping_rows=validate_sample_reference_rows(codon_mapping_rows,Path(mapping_table).name)
  counts['mappings']=len(codon_mapping_rows)
  writers[mapping_table].writerows(codon_mapping_rows)
 print(f"Wrote {counts['codons']} reference coding rows and {counts['mappings']} sample/reference mappings from {counts['references']} references.")
def run_reference(ref,linked,paths,settings,a):
 """Execute one MITOS2 reference and always return its materialized result."""
 fasta=ref['mitos2_input_fasta'];raw=Path(paths['mitos2_raw_dir'])/ref.get('task_key',ref['reference_key']);raw.mkdir(parents=True,exist_ok=True)
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
  tasks=task_rows(refs,paths);write(task_path,TASK_FIELDS,tasks)
  write(paths['sample_coordinate_reference_map'],SAMPLE_REFERENCE_FIELDS,sample_reference_rows(refs))
  print(f'Wrote {len(tasks)} MITOS2 reference tasks to {task_path}.');return
 if a.merge_only:
  merge(paths,settings,refs);write(task_path,TASK_FIELDS,task_rows(refs,paths));return
 if a.task_id:
  runnable=[pair for pair in refs if pair[0].get('status') != 'skipped_no_chrM_reference']
  selected=[pair for task,pair in zip(task_rows(refs,paths),runnable) if str(task['task_id'])==str(a.task_id)]
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
