#!/usr/bin/env python3
"""Validate codon tables and annotate coordinate-lifted VCF records."""
from __future__ import annotations
import argparse, csv, json, os, sys, tempfile, warnings
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from qc_analysis.lib.match_utils import (human_pos, info_format, info_parse, inject_headers,
    open_text, sample_names, source, write_summary, yaml)

CODON_MATCH_VERSION = "2.2"
RESOLVED_DNA = set("ACGT")
DNA = RESOLVED_DNA  # Backwards-compatible name used by codon helpers.
IUPAC_DNA = set("ACGTRYSWKMBDHVN")
STRING_FIELDS = ['MTCODON_STATUS','MTCODON_SUPPORTED_SNV','MTCODON_MATCH','MTCODON_STRICT_PHASE',
'MTCODON_GENE_MATCH','MTCODON_PHASE_MATCH','MTCODON_PRIMATE_GENE','MTCODON_PRIMATE_CODON',
'MTCODON_PRIMATE_ALT_CODON','MTCODON_PRIMATE_PHASE','MTCODON_HUMAN_GENE','MTCODON_HUMAN_CODON',
'MTCODON_HUMAN_PHASE','MTCODON_OVERLAPPING_CDS','MTCODON_AMBIGUOUS_BEST_MATCH',
'MTCODON_SOURCE_REF_MATCH','MTCODON_DUPLICATE_ANNOTATIONS','MTCODON_SOURCE_CODON_RESOLVED',
'MTCODON_HUMAN_CODON_RESOLVED','MTCODON_ANY_RESOLVED_PAIR','MTCODON_SOURCE_REF_RESOLVED',
'MTCODON_ANY_RESOLVED_SOURCE_REF']
INTEGER_FIELDS = ['MTCODON_N_PRIMATE_ANNOTATIONS','MTCODON_N_HUMAN_ANNOTATIONS','MTCODON_N_PAIR_CANDIDATES']
MULTI_FIELDS = ['MTCODON_PRIMATE_GENES','MTCODON_HUMAN_GENES','MTCODON_MATCHING_GENES']
DESCRIPTIONS = {
'MTCODON_STATUS':'Codon match status','MTCODON_SUPPORTED_SNV':'Whether source alleles form a simple SNV',
'MTCODON_MATCH':'Whether a gene/phase-matched human codon matches source reference or alternate codon',
'MTCODON_STRICT_PHASE':'Whether strict gene/phase failure-status categorization is enabled',
'MTCODON_SOURCE_REF_MATCH':'Whether SRC_REF matches the codon-table genomic reference base',
'MTCODON_DUPLICATE_ANNOTATIONS':'Whether duplicate source or human annotations were removed at this position'}
FIELDS = ([(x,'1','String',DESCRIPTIONS.get(x,x.replace('MTCODON_','').replace('_',' ').title())) for x in STRING_FIELDS]
          +[(x,'1','Integer',x.replace('MTCODON_N_','Number of ').replace('_',' ').lower()) for x in INTEGER_FIELDS]
          +[(x,'.','String','Sorted unique '+x.replace('MTCODON_','').replace('_',' ').lower()) for x in MULTI_FIELDS])
REQ_REFERENCE={'reference_key','pos','gene','strand','codon_pos_in_triplet','codon_seq','ref_base_genome'}
REQ_HISTORICAL={'sample','pos','gene','strand','codon_pos_in_triplet','codon_seq','ref_base_genome'}
REQ_HUMAN={'pos','gene','strand','codon_pos_in_triplet','codon_seq'}
REQ_MAP={'sample','reference_key'}

@dataclass(frozen=True, slots=True)
class CodonAnnotation:
    gene:str; strand:str; codon_pos_in_triplet:str; codon_seq:str; ref_base_genome:str=''
    def get(self,name,default=None): return getattr(self,name,default)
    def __getitem__(self,name): return getattr(self,name)

class CodonIndex(defaultdict):
    def __init__(self):
        super().__init__(list); self.duplicates=Counter(); self.duplicate_details=[]; self.ambiguous_details=[]
        self.ambiguous_ref_details=[]; self.resolved_ref_rows=0; self.conflicting_resolved_ref_positions=0; self.loaded=0

def normalize_codon(codon): return str(codon or '').strip().upper()
def normalize_base(base): return str(base or '').strip().upper()
def is_valid_iupac_base(base):
    base=normalize_base(base); return len(base)==1 and base in IUPAC_DNA
def is_resolved_base(base):
    base=normalize_base(base); return len(base)==1 and base in RESOLVED_DNA
def is_resolved_codon(codon):
    codon=normalize_codon(codon); return len(codon)==3 and set(codon)<=DNA
def is_valid_iupac_codon(codon):
    codon=normalize_codon(codon); return len(codon)==3 and set(codon)<=IUPAC_DNA

def _atomic_text(path, writer):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    fd,tmp=tempfile.mkstemp(prefix=f'.{path.name}.',suffix='.tmp',dir=path.parent,text=True)
    try:
        with os.fdopen(fd,'w',newline='') as h: writer(h); h.flush(); os.fsync(h.fileno())
        os.replace(tmp,path)
    except BaseException:
        try: os.unlink(tmp)
        except FileNotFoundError: pass
        raise

def _schema_error(path,kind,required,observed):
    missing=sorted(required-set(observed or []))
    raise SystemExit(f"Invalid {kind} table {path}: required columns={sorted(required)}; missing columns={missing}; observed columns={observed or []}")

def _validate_annotation(row,path,kind,line,key_column,strict=True):
    if not strict:return
    def bad(msg): raise SystemExit(f"Invalid {kind} table {path}, row {line}: {msg}; values={row}")
    try: pos=int(row['pos'])
    except (ValueError,TypeError): bad('pos must be a positive integer')
    if pos <= 0: bad('pos must be a positive integer')
    if row['strand'].strip() not in {'+','-'}: bad('strand must be + or -')
    if row['codon_pos_in_triplet'].strip() not in {'1','2','3'}: bad('codon_pos_in_triplet must be 1, 2, or 3')
    codon=normalize_codon(row['codon_seq'])
    if not is_valid_iupac_codon(codon): bad('codon_seq must contain exactly three valid IUPAC DNA symbols (ACGTRYSWKMBDHVN)')
    base=normalize_base(row.get('ref_base_genome'))
    if 'ref_base_genome' in required_for_kind(kind) and not base: bad('ref_base_genome must be non-empty')
    if base and not is_valid_iupac_base(base): bad('ref_base_genome must be exactly one valid IUPAC DNA symbol (ACGTRYSWKMBDHVN)')
    if key_column and not (row.get(key_column) or '').strip(): bad(f'{key_column} must be non-empty')

def required_for_kind(kind):
    return REQ_HUMAN if kind == 'human codon' else REQ_REFERENCE if kind == 'reference codon' else REQ_HISTORICAL

def load_codon_index(path,key_column=None,table_type=None,strict=True):
    kind=table_type or ('reference codon' if key_column=='reference_key' else 'historical sample codon' if key_column=='sample' else 'human codon')
    required=REQ_REFERENCE if key_column=='reference_key' else REQ_HISTORICAL if key_column=='sample' else REQ_HUMAN
    index=CodonIndex(); seen=defaultdict(dict); bases=defaultdict(set)
    with open_text(path) as handle:
        reader=csv.DictReader(handle,delimiter='\t'); observed=reader.fieldnames or []
        if not required.issubset(observed): _schema_error(path,kind,required,observed)
        for line,row in enumerate(reader,2):
            _validate_annotation(row,path,kind,line,key_column,strict)
            try: pos=int(row['pos'])
            except (ValueError,TypeError):
                if strict: raise
                continue
            key=(row.get(key_column,'').strip() if key_column else '',pos)
            values=(row.get('gene','').strip(),row.get('strand','').strip(),row.get('codon_pos_in_triplet','').strip(),
                    row.get('codon_seq','').strip().upper(),row.get('ref_base_genome','').strip().upper())
            optional=tuple((n,(row.get(n) or '').strip()) for n in ('codon_index','codon_pos1','codon_pos2','codon_pos3') if n in row)
            signature=values+optional
            if signature in seen[key]:
                index.duplicates[key]+=1; seen[key][signature]+=1; continue
            seen[key][signature]=1; index[key].append(CodonAnnotation(*values)); index.loaded+=1
            if not is_resolved_codon(values[3]):
                index.ambiguous_details.append({'table':kind,'reference_key':key[0],'position':pos,
                    'gene':values[0],'codon_seq':values[3],'annotation_source':(row.get('annotation_source') or '').strip(),
                    'file_name':str(path)})
            if values[4]: bases[key].add(values[4])
            if is_resolved_base(values[4]): index.resolved_ref_rows+=1
            elif is_valid_iupac_base(values[4]):
                index.ambiguous_ref_details.append({'table':kind,'reference_key':key[0],'position':pos,
                    'gene':values[0],'ref_base_genome':values[4],'codon_seq':values[3],
                    'annotation_source':(row.get('annotation_source') or '').strip(),'file_name':str(path)})
    inconsistent={k:{b for b in v if is_resolved_base(b)} for k,v in bases.items()
                  if len({b for b in v if is_resolved_base(b)})>1}
    index.conflicting_resolved_ref_positions=len(inconsistent)
    if inconsistent:
        example=next(iter(inconsistent.items()))
        raise SystemExit(f'Inconsistent ref_base_genome values in {kind} table {path}: position {example[0]} has {sorted(example[1])}; inconsistent positions={len(inconsistent)}')
    for key,sigs in seen.items():
        for sig,count in sigs.items():
            if count>1:index.duplicate_details.append((kind,key,sig[0],count-1,repr(sig)))
    return index

def load(path,key_column=None): return load_codon_index(path,key_column),key_column

def _genes(rows): return sorted({_get(r,'gene','').strip() for r in rows if _get(r,'gene','').strip()})
def find_overlapping_annotations(path,table,key_column=None):
    idx=load_codon_index(path,key_column)
    return sorted([{'table':table,'reference_key':k,'position':p,'genes':','.join(_genes(a)),
                    'number_of_annotations':len(a)} for (k,p),a in idx.items() if len(_genes(a))>1],
                  key=lambda r:(r['table'],r['reference_key'],r['position']))

def load_sample_reference_map(path,strict=True,stats=None):
    mapping={}; conflicts=0
    with open_text(path) as h:
        reader=csv.DictReader(h,delimiter='\t'); observed=reader.fieldnames or []
        if not REQ_MAP.issubset(observed): _schema_error(path,'sample-reference map',REQ_MAP,observed)
        for line,row in enumerate(reader,2):
            sample=(row.get('sample') or '').strip(); ref=(row.get('reference_key') or '').strip()
            if strict and (not sample or not ref): raise SystemExit(f'Invalid sample-reference map {path}, row {line}: sample and reference_key must be non-empty; values={row}')
            if not sample or not ref: continue
            if sample in mapping and mapping[sample]!=ref:
                conflicts+=1; raise SystemExit(f"Conflicting reference keys for sample {sample!r}: {mapping[sample]!r} versus {ref!r}")
            mapping[sample]=ref
    if stats is not None: stats['conflicting_mappings']=conflicts
    return mapping

def _get(row,name,default='.'): return row.get(name,default) if row else default
def _row_tie_break(row): return tuple(str(_get(row,n,'')) for n in ('gene','codon_seq','codon_pos_in_triplet','strand'))
def complement_base(base): return {'A':'T','T':'A','C':'G','G':'C'}.get(str(base).upper(),str(base).upper())
def is_supported_snv(ref,alt):
    ref=str(ref or '').strip().upper(); alt=str(alt or '').strip().upper()
    return len(ref)==len(alt)==1 and ref in DNA and alt in DNA

def mutate_codon(codon,phase,alt_base):
    try: phase=int(phase)
    except (TypeError,ValueError): return '.'
    bases=list(str(codon or '').upper()); alt=str(alt_base or '').strip().upper()
    if not is_resolved_codon(codon) or phase not in (1,2,3) or alt not in DNA:return '.'
    bases[phase-1]=alt; return ''.join(bases)

def evaluate_candidates(source_rows,human_rows,source_alt,supported_snv=True,source_ref=''):
    result=[]
    for s in source_rows:
        strand_alt=complement_base(source_alt) if supported_snv and _get(s,'strand','+')=='-' else source_alt
        sg=_get(s,'gene',''); sp=str(_get(s,'codon_pos_in_triplet','')); sc=_get(s,'codon_seq','.')
        sb=normalize_base(_get(s,'ref_base_genome','')); sr=normalize_base(source_ref)
        ref_resolved=is_resolved_base(sb); ref_match=ref_resolved and is_resolved_base(sr) and sb==sr
        ac=mutate_codon(sc,sp,strand_alt) if supported_snv and ref_match else '.'
        for h in human_rows:
            hg=_get(h,'gene',''); hp=str(_get(h,'codon_pos_in_triplet','')); hc=_get(h,'codon_seq','')
            result.append({'source':s,'human':h,'source_gene':sg,'human_gene':hg,'source_phase':sp,'human_phase':hp,
            'source_codon':sc,'human_codon':hc,'alternate_codon':ac,'gene_match':sg==hg,'phase_match':sp==hp,
            'source_resolved':is_resolved_codon(sc),'human_resolved':is_resolved_codon(hc),
            'source_ref_base':sb,'source_ref_resolved':ref_resolved,'source_ref_match':ref_match,
            'codon_match':supported_snv and ref_match and is_resolved_codon(sc) and is_resolved_codon(hc) and hc in {sc,ac}})
    return result

def _score(c): return (int(c['gene_match']),int(c['gene_match'] and c['phase_match']),int(c['gene_match'] and c['phase_match'] and c['codon_match']))
def _tie(c): return tuple(str(c[n]) for n in ('source_gene','human_gene','source_codon','human_codon','source_phase','human_phase','alternate_codon'))

def annotate(source_rows,human_rows,source_alt,strict,source_ref='A',source_duplicate=False,human_duplicate=False):
    supported=is_supported_snv(source_ref,source_alt)
    candidates=evaluate_candidates(source_rows,human_rows,source_alt,supported,source_ref)
    sg,hg=_genes(source_rows),_genes(human_rows); valid=[c for c in candidates if c['gene_match'] and c['phase_match'] and c['codon_match']]
    any_gene=any(c['gene_match'] for c in candidates); any_phase=any(c['gene_match'] and c['phase_match'] for c in candidates)
    vals={'MTCODON_SUPPORTED_SNV':'yes' if supported else 'no','MTCODON_STRICT_PHASE':'yes' if strict else 'no',
    'MTCODON_DUPLICATE_ANNOTATIONS':'yes' if source_duplicate or human_duplicate else 'no',
    'MTCODON_N_PRIMATE_ANNOTATIONS':str(len(source_rows)),'MTCODON_N_HUMAN_ANNOTATIONS':str(len(human_rows)),
    'MTCODON_N_PAIR_CANDIDATES':str(len(candidates)),'MTCODON_OVERLAPPING_CDS':'yes' if len(sg)>1 or len(hg)>1 else 'no',
    'MTCODON_PRIMATE_GENES':','.join(sg) or '.','MTCODON_HUMAN_GENES':','.join(hg) or '.',
    'MTCODON_GENE_MATCH':'yes' if any_gene else 'no','MTCODON_PHASE_MATCH':'yes' if any_phase else 'no',
    'MTCODON_MATCH':'yes' if valid else 'no','MTCODON_MATCHING_GENES':','.join(sorted({c['source_gene'] for c in valid})) or '.',
    'MTCODON_ANY_RESOLVED_PAIR':'yes' if any(c['source_resolved'] and c['human_resolved'] for c in candidates) else 'no',
    'MTCODON_ANY_RESOLVED_SOURCE_REF':'yes' if any(is_resolved_base(_get(r,'ref_base_genome','')) for r in source_rows) else 'no'}
    best=None; ambiguous=False
    if candidates:
        score=max(map(_score,candidates)); tied=[c for c in candidates if _score(c)==score]; ambiguous=len({_tie(c) for c in tied})>1; best=min(tied,key=_tie)
    vals['MTCODON_AMBIGUOUS_BEST_MATCH']='yes' if ambiguous else 'no'
    ss=best['source'] if best else min(source_rows,key=_row_tie_break,default=None); hh=best['human'] if best else min(human_rows,key=_row_tie_break,default=None)
    alt='.'
    selected_base=normalize_base(_get(ss,'ref_base_genome','')) if ss else ''
    selected_ref_resolved=is_resolved_base(selected_base)
    selected_ref_match=selected_ref_resolved and is_resolved_base(source_ref) and selected_base==normalize_base(source_ref)
    if ss and supported and selected_ref_match:
        a=complement_base(source_alt) if _get(ss,'strand','+')=='-' else source_alt; alt=mutate_codon(_get(ss,'codon_seq'),_get(ss,'codon_pos_in_triplet'),a)
    vals.update(MTCODON_PRIMATE_GENE=_get(ss,'gene'),MTCODON_PRIMATE_CODON=_get(ss,'codon_seq'),MTCODON_PRIMATE_ALT_CODON=alt,
    MTCODON_PRIMATE_PHASE=_get(ss,'codon_pos_in_triplet'),MTCODON_HUMAN_GENE=_get(hh,'gene'),MTCODON_HUMAN_CODON=_get(hh,'codon_seq'),MTCODON_HUMAN_PHASE=_get(hh,'codon_pos_in_triplet'))
    vals['MTCODON_SOURCE_CODON_RESOLVED']='NA' if not ss else ('yes' if is_resolved_codon(_get(ss,'codon_seq')) else 'no')
    vals['MTCODON_HUMAN_CODON_RESOLVED']='NA' if not hh else ('yes' if is_resolved_codon(_get(hh,'codon_seq')) else 'no')
    vals['MTCODON_SOURCE_REF_RESOLVED']='NA' if not ss or not selected_base else ('yes' if selected_ref_resolved else 'no')
    vals['MTCODON_SOURCE_REF_MATCH']='NA' if not ss or not selected_base or not selected_ref_resolved else ('yes' if selected_ref_match else 'no')
    return vals,candidates

def determine_status(source_pos,human_position,source_rows,human_rows,values,candidates,strict=True):
    if not source_pos or not human_position:return 'MISSING_COORD'
    if not source_rows:return 'SKIPPED_NONCODING'
    if not human_rows:return 'NO_HUMAN_CODON'
    compatible=[c for c in candidates if c['gene_match'] and c['phase_match']]
    if compatible and any(c['source_ref_resolved'] and not c['source_ref_match'] for c in compatible) and not any(c['source_ref_match'] for c in compatible):return 'SOURCE_REF_MISMATCH'
    if values['MTCODON_SUPPORTED_SNV']=='no':return 'UNSUPPORTED_NON_SNV'
    if values['MTCODON_MATCH']=='yes':return 'PASS'
    if compatible and all(is_valid_iupac_base(c['source_ref_base']) and not c['source_ref_resolved'] for c in compatible):return 'AMBIGUOUS_SOURCE_REF'
    if compatible and not any(c['source_resolved'] and c['human_resolved'] for c in compatible):return 'AMBIGUOUS_CODON'
    if strict and values['MTCODON_GENE_MATCH']=='no':return 'GENE_MISMATCH'
    if strict and values['MTCODON_PHASE_MATCH']=='no':return 'PHASE_MISMATCH'
    return 'MISMATCH'

def _strict_setting(settings):
    new=settings.get('strict_gene_phase_status'); old=settings.get('strict_phase_match')
    if new is not None and old is not None and bool(new)!=bool(old): warnings.warn('strict_gene_phase_status conflicts with legacy strict_phase_match; using strict_gene_phase_status',UserWarning)
    return bool(new if new is not None else old if old is not None else True)

def _write_rows(path,fieldnames,data):
    def writer(h):
        w=csv.DictWriter(h,fieldnames=fieldnames,delimiter='\t'); w.writeheader(); w.writerows(data)
    _atomic_text(path,writer)

def merge_summaries(reports_dir):
    reports=Path(reports_dir); files=sorted(p for p in reports.glob('*.codon_match_summary.tsv') if p.name!='all_samples.codon_match_summary.tsv')
    schema=None; by_sample={}
    for path in files:
        with path.open() as h:
            reader=csv.DictReader(h,delimiter='\t')
            if schema is None:schema=reader.fieldnames
            elif reader.fieldnames!=schema:raise SystemExit(f'Incompatible summary schema: {path}: {reader.fieldnames}; expected {schema}')
            for row in reader:
                sample=row.get('sample','')
                if sample in by_sample and by_sample[sample]!=row:raise SystemExit(f'Conflicting duplicate summary rows for sample {sample!r}')
                by_sample[sample]=row
    if not schema: raise SystemExit(f'No per-sample summaries found in {reports}')
    _write_rows(reports/'all_samples.codon_match_summary.tsv',schema,[by_sample[x] for x in sorted(by_sample)])
    print(f'Merged {len(files)} input files and {len(by_sample)} samples')

def _load_inputs(paths,strict):
    reference=paths.get('reference_codon_table'); mapping=paths.get('sample_reference_map'); refmode=bool(reference and mapping)
    primate=load_codon_index(reference,'reference_key',strict=strict) if refmode else load_codon_index(paths['all_primate_position_codon_table'],'sample',strict=strict)
    map_stats={}; maps=load_sample_reference_map(mapping,strict,map_stats) if refmode else {}
    human=load_codon_index(paths['human_codon_table'],strict=strict)
    return refmode,primate,maps,human,map_stats

def _diagnostics(paths,primate,human,report_overlaps=None):
    overlap=[]
    for table,index in [('source',primate),('human',human)]:
        for (key,pos),ann in index.items():
            genes=_genes(ann)
            if len(genes)>1:
                overlap.append({'table':table,'reference_key':key,'position':pos,'genes':','.join(genes),'annotation_count':len(ann),'unique_gene_count':len(genes),'duplicate_count':index.duplicates[(key,pos)]})
    if report_overlaps:_write_rows(report_overlaps,['table','reference_key','position','genes','annotation_count','unique_gene_count','duplicate_count'],overlap)
    details=[]
    for table,index in [('source',primate),('human',human)]:
        for kind,(key,pos),gene,count,sig in index.duplicate_details: details.append({'table':table,'reference_key':key,'position':pos,'gene':gene,'duplicate_count':count,'signature':sig})
    if details:_write_rows(Path(paths['reports_dir'])/'codon_annotation_duplicate_rows.tsv',['table','reference_key','position','gene','duplicate_count','signature'],details)
    ambiguous=[{**x,'table':'source'} for x in primate.ambiguous_details]+[{**x,'table':'human'} for x in human.ambiguous_details]
    if ambiguous:_write_rows(Path(paths['reports_dir'])/'codon_annotation_ambiguous_codons.tsv',
        ['table','reference_key','position','gene','codon_seq','annotation_source','file_name'],ambiguous)
    ambiguous_refs=[{**x,'table':'source'} for x in primate.ambiguous_ref_details]
    if ambiguous_refs:_write_rows(Path(paths['reports_dir'])/'codon_annotation_ambiguous_reference_bases.tsv',
        ['table','reference_key','position','gene','ref_base_genome','codon_seq','annotation_source','file_name'],ambiguous_refs)
    return overlap

def main():
    p=argparse.ArgumentParser(); p.add_argument('--config',required=True); p.add_argument('--sample'); p.add_argument('--input'); p.add_argument('--output')
    p.add_argument('--merge-summaries',action='store_true'); p.add_argument('--validate-inputs',action='store_true'); p.add_argument('--report-overlaps')
    args=p.parse_args(); cfg=yaml(args.config); section=cfg['codon_match']; paths=section['paths']; settings=section.get('settings',{})
    if args.merge_summaries: merge_summaries(paths['reports_dir']); return
    strict_validation=bool(settings.get('strict_input_validation',True)); strict_status=_strict_setting(settings)
    refmode,primate,maps,human,map_stats=_load_inputs(paths,strict_validation); overlaps=_diagnostics(paths,primate,human,args.report_overlaps)
    if args.validate_inputs:
        print(f'reference annotations loaded: {primate.loaded}\nhuman annotations loaded: {human.loaded}\nunique source index positions: {len(primate)}\nunique human positions: {len(human)}')
        print(f"source overlapping-gene positions: {sum(x['table']=='source' for x in overlaps)}\nhuman overlapping-gene positions: {sum(x['table']=='human' for x in overlaps)}")
        print(f'duplicate rows removed: {sum(primate.duplicates.values())+sum(human.duplicates.values())}\nsamples mapped: {len(maps)}\nconflicting mappings: {map_stats.get("conflicting_mappings",0)}\ninconsistent ref_base_genome positions: 0')
        for label,index in [('source',primate),('human',human)]:
            print(f'{label} ambiguous codon rows: {len(index.ambiguous_details)}')
            print(f'{label} ambiguous codon positions: {len({(x["reference_key"],x["position"]) for x in index.ambiguous_details})}')
        print(f'source_ambiguous_ref_base_rows: {len(primate.ambiguous_ref_details)}')
        print(f'source_ambiguous_ref_base_positions: {len({(x["reference_key"],x["position"]) for x in primate.ambiguous_ref_details})}')
        print(f'source_resolved_ref_base_rows: {primate.resolved_ref_rows}')
        print(f'source_conflicting_resolved_ref_base_positions: {primate.conflicting_resolved_ref_positions}')
        return
    samples=[args.sample] if args.sample else sample_names(cfg)
    if args.input:samples=[args.sample or Path(args.input).name.split('.')[0]]
    if not samples:raise SystemExit('No samples found; supply --sample or --input.')
    for sample in samples:
        if refmode:
            if sample not in maps:raise SystemExit(f'Sample {sample!r} is missing from sample_reference_map: {paths.get("sample_reference_map")}')
            reference_key=maps[sample]
        else:reference_key=sample
        inp=Path(args.input) if args.input else Path(paths['input_vcf_dir'])/str(settings['input_vcf_pattern']).format(sample=sample)
        out=Path(args.output) if args.output else Path(paths['output_dir'])/'vcf_codon'/f"{sample}{settings['output_suffix']}"
        if not inp.exists():raise SystemExit(f'Missing input VCF for {sample}: {inp}')
        header=[]; body=[]; counts=Counter(); metrics=Counter()
        with open_text(inp) as h:
            for line in h:
                if line.startswith('#'):header.append(line);continue
                f=line.rstrip('\n').split('\t'); info=info_parse(f[7]); _,sp,sref,salt=source(info); hp=human_pos(f,info)
                sr=primate.get((reference_key,sp),[]) if sp else []; hr=human.get(('',hp),[]) if hp else []
                sd=bool(sp and primate.duplicates[(reference_key,sp)]); hd=bool(hp and human.duplicates[('',hp)])
                vals,candidates=annotate(sr,hr,salt,strict_status,sref,sd,hd)
                status=determine_status(sp,hp,sr,hr,vals,candidates,strict_status)
                vals['MTCODON_STATUS']=status; info.update(vals); f[7]=info_format(info); body.append('\t'.join(f)+'\n'); counts[status]+=1
                sgo=len(_genes(sr))>1; hgo=len(_genes(hr))>1
                metrics['records_with_overlapping_source_cds']+=sgo; metrics['records_with_overlapping_human_cds']+=hgo; metrics['records_with_overlapping_cds']+=sgo or hgo
                metrics['records_with_ambiguous_best_match']+=vals['MTCODON_AMBIGUOUS_BEST_MATCH']=='yes'; metrics['records_with_multiple_pair_candidates']+=len(candidates)>1
                metrics['records_with_duplicate_source_annotations']+=sd; metrics['records_with_duplicate_human_annotations']+=hd; metrics['records_with_duplicate_annotations']+=sd or hd
                metrics['records_with_ambiguous_source_codon']+=any(not is_resolved_codon(_get(x,'codon_seq')) for x in sr)
                metrics['records_with_ambiguous_human_codon']+=any(not is_resolved_codon(_get(x,'codon_seq')) for x in hr)
                metrics['records_without_any_resolved_pair']+=vals['MTCODON_ANY_RESOLVED_PAIR']=='no'
        provenance=[f'##MTCODONVersion={CODON_MATCH_VERSION}\n',f'##MTCODONReferenceCodonTable={paths.get("reference_codon_table",paths.get("all_primate_position_codon_table"))}\n',f'##MTCODONHumanCodonTable={paths["human_codon_table"]}\n']
        existing=''.join(header); provenance=[x for x in provenance if x.split('=',1)[0] not in existing]
        _atomic_text(out,lambda h:(h.writelines(inject_headers(provenance+header,FIELDS,'MTCODON')),h.writelines(body)))
        statuses=['PASS','SKIPPED_NONCODING','NO_HUMAN_CODON','SOURCE_REF_MISMATCH','UNSUPPORTED_NON_SNV','AMBIGUOUS_SOURCE_REF','AMBIGUOUS_CODON','GENE_MISMATCH','PHASE_MISMATCH','MISMATCH','MISSING_COORD']
        metric_names=['records_with_overlapping_source_cds','records_with_overlapping_human_cds','records_with_overlapping_cds','records_with_ambiguous_best_match','records_with_multiple_pair_candidates','records_with_duplicate_source_annotations','records_with_duplicate_human_annotations','records_with_duplicate_annotations']
        metric_names += ['records_with_ambiguous_source_codon','records_with_ambiguous_human_codon','records_without_any_resolved_pair']
        row={'sample':sample,'input_vcf':str(inp),'output_vcf':str(out),'total_records':len(body),**{f'status_{x}':counts[x] for x in statuses},**{x:metrics[x] for x in metric_names},
        'codon_match_script_version':CODON_MATCH_VERSION,'reference_codon_table':str(paths.get('reference_codon_table',paths.get('all_primate_position_codon_table'))),'human_codon_table':str(paths['human_codon_table']),
        'sample_reference_map':str(paths.get('sample_reference_map','')),'reference_key':reference_key,'strict_gene_phase_status':strict_status,'strict_input_validation':strict_validation,'strict_phase_match':strict_status,
        'source_ambiguous_codon_rows':len(primate.ambiguous_details),'human_ambiguous_codon_rows':len(human.ambiguous_details),
        'source_ambiguous_codon_positions':len({(x['reference_key'],x['position']) for x in primate.ambiguous_details}),
        'human_ambiguous_codon_positions':len({(x['reference_key'],x['position']) for x in human.ambiguous_details}),
        'source_ambiguous_ref_base_rows':len(primate.ambiguous_ref_details),
        'source_ambiguous_ref_base_positions':len({(x['reference_key'],x['position']) for x in primate.ambiguous_ref_details}),
        'source_resolved_ref_base_rows':primate.resolved_ref_rows,
        'source_conflicting_resolved_ref_base_positions':primate.conflicting_resolved_ref_positions,'status':'completed'}
        write_summary(Path(paths['reports_dir'])/f'{sample}.codon_match_summary.tsv',row)
if __name__=='__main__':main()
