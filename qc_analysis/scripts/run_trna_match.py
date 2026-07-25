#!/usr/bin/env python3
"""Annotate lifted VCFs with tRNA structural comparisons.

Index bases are genomic-orientation bases.  For a ``-`` strand record both the
variant allele and ``paired_base`` are reverse-complemented before RNA pairing
is evaluated.
"""
import argparse, csv, sys, warnings
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from qc_analysis.lib.match_utils import *

REQUIRED = {'chrom','pos','trna_id','local_pos','struct_class','struct_element','pair_type',
            'pair_state','paired_local_pos','paired_genomic_pos','paired_base','strand'}
N=['STATUS','S_ID','H_ID','S_LOCAL','H_LOCAL','S_CLASS','H_CLASS','REGION_MATCH','S_ELEMENT','H_ELEMENT','ELEMENT_MATCH','S_PAIR_TYPE','H_PAIR_TYPE','PAIR_TYPE_MATCH','S_PAIR_STATE','H_PAIR_STATE','PAIR_STATE_MATCH','S_PAIR_LOCAL','H_PAIR_LOCAL','PAIR_LOCAL_MATCH','S_PAIR_POS','H_PAIR_POS','S_PAIR_LIFTED_HPOS','PAIR_POS_MATCH','H_ALT_PAIR_TYPE','S_ALT_PAIR_TYPE','H_ALT_EFFECT','S_ALT_EFFECT','ALLELE_EFFECT_MATCH','COMPENSATED','STRICT_MATCH','S_COORD_SPACE','S_LOOKUP_CHROM','S_LOOKUP_POS']
FIELDS=[('MTTRNA_'+n,'tRNA structural match annotation') for n in N]

def normalize_chrom(chrom, mode):
    """Normalize a chromosome using none, strip_chr, add_chr, or mitochondrial_alias."""
    value=str(chrom or '').strip(); mode=mode or 'none'
    if mode == 'none': return value
    if mode == 'strip_chr': return value[3:] if value.lower().startswith('chr') else value
    if mode == 'add_chr': return value if value.lower().startswith('chr') else 'chr'+value
    if mode == 'mitochondrial_alias':
        return 'MT' if value.lower().removeprefix('chr') in {'m','mt','mitochondria','mitochondrion'} else value
    raise ValueError(f'Unsupported chromosome normalization mode: {mode}')

@dataclass
class TrnaIndex:
    exact: dict
    positions: dict
    ambiguous_positions: set
    duplicate_keys: set
    chroms_by_pos: dict

    def lookup(self, chrom, pos, ignore_chrom=False):
        if pos is None: return None, 'missing_coordinate'
        record=self.exact.get((chrom,pos))
        if record is not None: return record, 'exact'
        if not ignore_chrom:
            return None, 'chromosome_mismatch' if pos in self.chroms_by_pos else 'not_found'
        if pos in self.ambiguous_positions: return None, 'ambiguous'
        return self.positions.get(pos), 'position' if pos in self.positions else 'not_found'

def index(path, chrom_norm='none'):
    data=rows(path)
    with open_text(path) as handle: fieldnames=next(csv.reader(handle,delimiter='\t'),[])
    missing=REQUIRED-set(fieldnames)
    if missing: raise ValueError(f'tRNA index {path} missing required columns: {", ".join(sorted(missing))}')
    exact={}; duplicates=set(); conflicts=[]; bypos=defaultdict(list); chroms=defaultdict(set)
    for number,r in enumerate(data,2):
        try: pos=int(r['pos'])
        except (TypeError,ValueError): raise ValueError(f'Invalid pos in {path} line {number}: {r.get("pos")!r}')
        chrom=normalize_chrom(r['chrom'],chrom_norm)
        if not chrom: raise ValueError(f'Blank chrom in {path} line {number}')
        r=dict(r); r['chrom']=chrom; key=(chrom,pos)
        if key in exact:
            if exact[key] == r: duplicates.add(key); continue
            conflicts.append(key); continue
        exact[key]=r; bypos[pos].append(r); chroms[pos].add(chrom)
    if conflicts: raise ValueError(f'Conflicting duplicate chrom/pos keys in {path}: {sorted(set(conflicts))}')
    if duplicates: warnings.warn(f'Collapsed {len(duplicates)} exact duplicate chrom/pos key(s) in {path}', RuntimeWarning)
    ambiguous={pos for pos,rs in bypos.items() if len(rs)>1}
    positions={pos:rs[0] for pos,rs in bypos.items() if len(rs)==1}
    return TrnaIndex(exact,positions,ambiguous,duplicates,dict(chroms))

def map_for(directory,sample):
    for suffix in ('.coordinate_map.tsv','.coordinate_map.tsv.gz'):
        path=Path(directory)/f'{sample}{suffix}'
        if path.exists(): return load_coordinate_map(path)
    return {}

def oriented(base, record):
    b=normalize_rna_base(base)
    if not b: return base
    return {'A':'U','U':'A','C':'G','G':'C'}[b] if record.get('strand') == '-' else b

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--config',required=True); ap.add_argument('--sample'); ap.add_argument('--input'); ap.add_argument('--output'); a=ap.parse_args()
    c=yaml(a.config); sec=c['trna_match']; p,s=sec['paths'],sec['settings']
    if s.get('run_trnascan_if_missing',False): raise SystemExit('Unsupported feature: run_trnascan_if_missing=true; tRNAscan index generation is unavailable')
    if s.get('run_trna_gene_qc',False): raise SystemExit('Unsupported feature: run_trna_gene_qc=true; interval-level tRNA gene QC is not implemented')
    hpath=Path(p['human_trna_index'])
    if not hpath.exists(): raise SystemExit(f'Missing human tRNA index: {hpath}')
    hi=index(hpath,s.get('human_trna_chrom_norm','none'))
    samples=[a.sample] if a.sample else sample_names(c)
    if a.input: samples=[a.sample or Path(a.input).name.split('.')[0]]
    if not samples or any(not x for x in samples): raise SystemExit('No samples selected; provide --sample/--input or a non-empty sample reference file')
    allrows=[]
    for sample in samples:
        primary=Path(p['input_vcf_dir'])/str(s['input_vcf_pattern']).format(sample=sample); fallback=Path(p['fallback_input_vcf_dir'])/str(s['fallback_input_vcf_pattern']).format(sample=sample)
        inp=Path(a.input) if a.input else (primary if primary.exists() else fallback); codon=inp==primary
        if not inp.exists(): raise SystemExit(f'Missing input VCF for sample {sample}; attempted primary {primary} and fallback {fallback}' + (f'; explicit input {inp}' if a.input else ''))
        spi=Path(str(p['species_trna_index_template']).format(species_trna_index_dir=p['species_trna_index_dir'],sample=sample))
        if not spi.exists(): raise SystemExit(f'Missing species tRNA index for {sample}: {spi}')
        si=index(spi,s.get('species_trna_chrom_norm','none')); cmap=map_for(p.get('coordinate_map_dir',''),sample)
        out=Path(a.output) if a.output else Path(p['output_dir'])/'vcf_trna'/f"{sample}{s['output_suffix'] if codon else '.lifted.trna.vcf'}"; out.parent.mkdir(parents=True,exist_ok=True)
        head=[]; body=[]; counts=Counter()
        for lineno,line in enumerate(open_text(inp),1):
            if line.startswith('#'): head.append(line); continue
            x=line.rstrip('\n').split('\t')
            if len(x)<8: raise ValueError(f'VCF {inp} line {lineno} has {len(x)} fields; expected at least 8')
            try: int(x[1])
            except ValueError: raise ValueError(f'Invalid VCF POS in {inp} line {lineno}: {x[1]!r}')
            inf=info_parse(x[7]); sch,pos,_,alt=source(inf); hp=human_pos(x,inf)
            source_pos_raw=inf.get('SRC_POS',inf.get('MTLIFT_ORIG_POS'))
            if source_pos_raw not in (None,'') and pos is None: raise ValueError(f'Invalid INFO source coordinate in {inp} line {lineno}: {source_pos_raw!r}')
            human_pos_raw=inf.get('MTLIFT_HUMAN_POS')
            if human_pos_raw not in (None,'') and hp is None: raise ValueError(f'Invalid INFO human coordinate in {inp} line {lineno}: {human_pos_raw!r}')
            sch=normalize_chrom(sch,s.get('species_vcf_chrom_norm','none')); hch=normalize_chrom(x[0],s.get('human_vcf_chrom_norm','none'))
            sr,skind=si.lookup(sch,pos,bool(s.get('species_trna_lookup_ignore_chrom',False))); hr,hkind=hi.lookup(hch,hp,bool(s.get('human_trna_lookup_ignore_chrom',False)))
            counts['ambiguous_species_index_lookup'] += skind=='ambiguous'; counts['ambiguous_human_index_lookup'] += hkind=='ambiguous'; counts['chromosome_mismatch'] += skind=='chromosome_mismatch' or hkind=='chromosome_mismatch'
            if not pos: status='MISSING_SPECIES_COORD'
            elif not sr and not hr: status='NO_SPECIES_OR_HUMAN_TRNA'
            elif not sr: status='NO_SPECIES_TRNA'
            elif not hr: status='NO_HUMAN_TRNA'
            else: status='OK'
            v={'MTTRNA_STATUS':status,'MTTRNA_S_COORD_SPACE':s.get('species_trna_coord_space','original'),'MTTRNA_S_LOOKUP_CHROM':sch or '.','MTTRNA_S_LOOKUP_POS':pos or '.'}
            for short,col in [('S_ID','trna_id'),('H_ID','trna_id'),('S_LOCAL','local_pos'),('H_LOCAL','local_pos'),('S_CLASS','struct_class'),('H_CLASS','struct_class'),('S_ELEMENT','struct_element'),('H_ELEMENT','struct_element'),('S_PAIR_TYPE','pair_type'),('H_PAIR_TYPE','pair_type'),('S_PAIR_STATE','pair_state'),('H_PAIR_STATE','pair_state'),('S_PAIR_LOCAL','paired_local_pos'),('H_PAIR_LOCAL','paired_local_pos'),('S_PAIR_POS','paired_genomic_pos'),('H_PAIR_POS','paired_genomic_pos')]:
                record=(sr or {}) if short.startswith('S') else (hr or {}); v['MTTRNA_'+short]=record.get(col,'.')
            for key,col in [('REGION_MATCH','struct_class'),('ELEMENT_MATCH','struct_element'),('PAIR_TYPE_MATCH','pair_type'),('PAIR_STATE_MATCH','pair_state'),('PAIR_LOCAL_MATCH','paired_local_pos')]: v['MTTRNA_'+key]=compare_values(sr.get(col) if sr else '.',hr.get(col) if hr else '.')
            stem=status=='OK' and v['MTTRNA_S_CLASS']=='stem' and v['MTTRNA_H_CLASS']=='stem'; lifted=lift_source_pos_to_human(v['MTTRNA_S_PAIR_POS'],cmap) if stem else '.'
            if stem and lifted=='.': counts['missing_coordinate_map']+=1
            posmatch=compare_values(lifted,v['MTTRNA_H_PAIR_POS']) if stem else '.'
            salt=alt or x[4]; halt=x[4]
            spt=pair_type(oriented(salt,sr),oriented(sr.get('paired_base'),sr)) if stem else '.'; hpt=pair_type(oriented(halt,hr),oriented(hr.get('paired_base'),hr)) if stem else '.'
            if (sr and sr.get('strand')=='-') or (hr and hr.get('strand')=='-'): counts['negative_strand_records']+=1
            seffect=pair_effect(v['MTTRNA_S_PAIR_TYPE'],spt) if stem else '.'; heffect=pair_effect(v['MTTRNA_H_PAIR_TYPE'],hpt) if stem else '.'; effectmatch=compare_values(seffect,heffect) if stem else '.'
            compatible={'WC','GU_wobble'}; compensated='yes' if stem and spt in compatible and hpt in compatible else 'no' if stem else '.'; strict='no'
            if status=='OK' and v['MTTRNA_S_CLASS']=='loop' and v['MTTRNA_H_CLASS']=='loop': strict='yes' if v['MTTRNA_REGION_MATCH']=='yes' and v['MTTRNA_ELEMENT_MATCH']=='yes' and compare_values(v['MTTRNA_S_LOCAL'],v['MTTRNA_H_LOCAL'])=='yes' else 'no'
            elif stem:
                checks=[v['MTTRNA_REGION_MATCH']=='yes',v['MTTRNA_ELEMENT_MATCH']=='yes',v['MTTRNA_PAIR_STATE_MATCH']=='yes',posmatch=='yes',effectmatch=='yes']
                if s.get('require_compensated_for_strict_stem',True): checks.append(compensated=='yes')
                strict='yes' if all(checks) else 'no'
            v.update({'MTTRNA_S_PAIR_LIFTED_HPOS':lifted,'MTTRNA_PAIR_POS_MATCH':posmatch,'MTTRNA_H_ALT_PAIR_TYPE':hpt,'MTTRNA_S_ALT_PAIR_TYPE':spt,'MTTRNA_H_ALT_EFFECT':heffect,'MTTRNA_S_ALT_EFFECT':seffect,'MTTRNA_ALLELE_EFFECT_MATCH':effectmatch,'MTTRNA_COMPENSATED':compensated,'MTTRNA_STRICT_MATCH':strict})
            inf.update(v); x[7]=info_format(inf)
            if not s.get('pass_only',False) or x[6]=='PASS': body.append('\t'.join(x)+'\n'); counts[status]+=1
        with out.open('w') as f: f.writelines(inject_headers(head,FIELDS,'MTTRNA')); f.writelines(body)
        row={'sample':sample,'input_vcf':str(inp),'output_vcf':str(out),'total_records':len(body),**{f'status_{q}':counts[q] for q in ['OK','NO_SPECIES_TRNA','NO_HUMAN_TRNA','NO_SPECIES_OR_HUMAN_TRNA','MISSING_SPECIES_COORD']},**{q:counts[q] for q in ['ambiguous_species_index_lookup','ambiguous_human_index_lookup','chromosome_mismatch','missing_coordinate_map','negative_strand_records']},'status':'completed'}
        if s.get('write_summary',True): write_summary(Path(p['reports_dir'])/f'{sample}.trna_match_summary.tsv',row)
        allrows.append(row)
    if s.get('write_summary',True) and allrows and not a.sample and not a.input:
        q=Path(p['reports_dir'])/'all_samples.trna_match_summary.tsv'; q.parent.mkdir(parents=True,exist_ok=True)
        with q.open('w',newline='') as f: w=csv.DictWriter(f,fieldnames=list(allrows[0]),delimiter='\t'); w.writeheader(); w.writerows(allrows)
if __name__=='__main__': main()
