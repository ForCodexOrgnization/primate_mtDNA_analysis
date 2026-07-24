#!/usr/bin/env python3
"""Compare GenBank and MITOS2 CDS tables at coordinate-reference level only."""
from __future__ import annotations
import argparse, csv, os, tempfile
from collections import defaultdict
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from qc_analysis.lib.match_utils import yaml
from qc_analysis.lib.reference_utils import normalized_fasta_sequence_sha256

GENES = ('MT-ND1 MT-ND2 MT-ND3 MT-ND4 MT-ND4L MT-ND5 MT-ND6 MT-CO1 MT-CO2 MT-CO3 MT-CYB MT-ATP6 MT-ATP8').split()
GENE_FIELDS = '''genbank_reference_key mitos2_reference_key gene genbank_present mitos2_present genbank_n_rows mitos2_n_rows genbank_n_unique_positions mitos2_n_unique_positions genbank_start genbank_end mitos2_start mitos2_end genbank_strand mitos2_strand position_overlap position_union position_jaccard genbank_only_positions mitos2_only_positions same_position_set same_strand same_codon_triplets same_codon_position_mapping same_reference_bases start_delta end_delta length_delta wraps_origin_genbank wraps_origin_mitos2 ordered_coordinate_match sequence_compatibility_category gene_comparison_category'''.split()
SUMMARY_FIELDS = '''genbank_reference_key mitos2_reference_key target_species reference_species genbank_accession mitos2_accession genbank_coordinate_fasta mitos2_coordinate_fasta genbank_coordinate_sequence_length mitos2_coordinate_sequence_length genbank_coordinate_sequence_sha256 mitos2_coordinate_sequence_sha256 genbank_record_sequence_sha256 mitos2_input_sequence_sha256 sequence_compatibility_category sequence_match rotation_equivalent rotation_offset coordinate_comparison_performed n_genes_genbank n_genes_mitos2 n_exact_genes n_minor_difference_genes n_moderate_difference_genes n_major_difference_genes n_missing_genes n_strand_mismatches all_13_genes_present_genbank all_13_genes_present_mitos2 all_13_exact_match reference_comparison_category note'''.split()
DIAG_FIELDS = '''genbank_reference_key mitos2_reference_key target_species reference_species genbank_accession mitos2_accession genbank_coordinate_fasta mitos2_coordinate_fasta genbank_coordinate_sequence_length mitos2_coordinate_sequence_length genbank_coordinate_sequence_sha256 mitos2_coordinate_sequence_sha256 genbank_record_sequence_sha256 mitos2_input_sequence_sha256 sequence_compatibility_category candidate_matching_basis rotation_check_result rotation_offset reason_comparison_skipped'''.split()
def v(r,k): return (r.get(k) or '').strip()
def yes(x): return 'yes' if x else 'no'
def read_groups(path):
    p=Path(path)
    if not p.is_file(): raise FileNotFoundError(f'Required reference codon table is missing: {p}')
    groups=defaultdict(list)
    with p.open(newline='') as h:
        for r in csv.DictReader(h,delimiter='\t'):
            key=v(r,'reference_key') or (v(r,'coordinate_reference_accession')+'|'+v(r,'coordinate_reference_fasta'))
            r['reference_key']=key; groups[key].append(r)
    return groups
def meta(rows, mitos=False):
    r=rows[0] if rows else {}; path=v(r,'coordinate_reference_fasta') or (v(r,'mitos2_input_fasta') if mitos else '')
    h=v(r,'coordinate_reference_sequence_sha256') or (v(r,'mitos2_input_sequence_sha256') if mitos else '')
    length=v(r,'mitos2_input_sequence_length') if mitos else v(r,'coordinate_reference_sequence_length')
    # Tables from older pipeline versions are upgraded from their explicitly recorded FASTA.
    if not h and path:
        info=normalized_fasta_sequence_sha256(path); h=info['sequence_sha256']; length=str(info['sequence_length'])
    return {'key':v(r,'reference_key'),'path':path,'hash':h,'length':length,
            'accession':v(r,'coordinate_reference_accession') or v(r,'accession_version') or v(r,'accession'),
            'species':v(r,'reference_species') or v(r,'species') or v(r,'target_species'),
            'target':v(r,'target_species') or v(r,'species'),
            'record_hash':v(r,'genbank_record_sequence_sha256'), 'input_hash':v(r,'mitos2_input_sequence_sha256') or h}
def rotation(a,b):
    if not a or not b or len(a)!=len(b): return None
    return (b+b).find(a)
def ordered(rows):
    return [(v(r,'pos'),v(r,'codon_index'),v(r,'codon_pos_in_triplet'),v(r,'strand')) for r in sorted(rows,key=lambda x:(int(v(x,'codon_index') or 0),int(v(x,'codon_pos_in_triplet') or 0),int(v(x,'pos') or 0)))]
def gene_row(g,m,gene,cat,minor,moderate):
    a=[r for r in g if v(r,'gene')==gene]; b=[r for r in m if v(r,'gene')==gene]; ap={v(r,'pos') for r in a}; bp={v(r,'pos') for r in b}
    strand_a={v(r,'strand') for r in a}; strand_b={v(r,'strand') for r in b}; overlap=len(ap&bp); union=len(ap|bp); jac=overlap/union if union else 1.0
    wrapa=bool(a and ordered(a)[0][0] != min(ap,key=lambda x:int(x))); wrapb=bool(b and ordered(b)[0][0] != min(bp,key=lambda x:int(x)))
    base=lambda xs:{v(x,'pos'):v(x,'ref_base_genome') for x in xs}
    samepos=ap==bp; samestrand=strand_a==strand_b and len(strand_a)==1; ordermatch=ordered(a)==ordered(b)
    trip=lambda xs:[(v(x,'codon_index'),v(x,'codon_pos_in_triplet'),v(x,'codon_seq')) for x in sorted(xs,key=lambda x:(int(v(x,'codon_index') or 0),int(v(x,'codon_pos_in_triplet') or 0))) ]
    mapping=lambda xs:[(v(x,'pos'),v(x,'codon_index'),v(x,'codon_pos_in_triplet')) for x in xs]
    if not a: kind='missing_in_genbank'
    elif not b: kind='missing_in_mitos2'
    elif not samestrand: kind='strand_mismatch'
    elif samepos and ordermatch and trip(a)==trip(b) and base(a)==base(b): kind='exact_match'
    elif samepos: kind='coordinate_match_codon_difference'
    elif jac>=minor: kind='minor_boundary_difference'
    elif jac>=moderate: kind='moderate_boundary_difference'
    else: kind='major_coordinate_difference'
    def bounds(x): return (min((int(v(r,'pos')) for r in x),default=''),max((int(v(r,'pos')) for r in x),default=''))
    sa,ea=bounds(a); sb,eb=bounds(b)
    return dict(zip(GENE_FIELDS,[v(g[0],'reference_key'),v(m[0],'reference_key'),gene,yes(a),yes(b),len(a),len(b),len(ap),len(bp),sa,ea,sb,eb,','.join(sorted(strand_a)),','.join(sorted(strand_b)),overlap,union,f'{jac:.12g}',','.join(sorted(ap-bp,key=int)),','.join(sorted(bp-ap,key=int)),yes(samepos),yes(samestrand),yes(trip(a)==trip(b)),yes(mapping(a)==mapping(b)),yes(base(a)==base(b)),(sb-sa if a and b else ''),(eb-ea if a and b else ''),(len(bp)-len(ap)),yes(wrapa),yes(wrapb),yes(ordermatch),cat,kind]))
def atomic(path, fields, rows):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    fd,tmp=tempfile.mkstemp(dir=path.parent,prefix='.'+path.name+'.',text=True)
    with os.fdopen(fd,'w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
    if not Path(tmp).read_text().splitlines()[0].split('\t')==fields: raise RuntimeError('Atomic output validation failed')
    os.replace(tmp,path)
def compare(genbank_table,mitos_table,output,summary_output,mismatch_output,strict=True,allow_rotation=False,fail_no_shared=True,minor=.99,moderate=.90):
    gb=read_groups(genbank_table); mi=read_groups(mitos_table); gm={k:meta(x) for k,x in gb.items()}; mm={k:meta(x,True) for k,x in mi.items()}
    candidates=[]
    for g, x in gm.items():
        for m,y in mm.items():
            if x['hash'] and y['hash'] and x['hash']==y['hash']: candidates.append((g,m,'sequence_sha256'))
            elif g==m or (x['accession'] and x['accession']==y['accession']): candidates.append((g,m,'reference_key_or_accession'))
    # no species-only matching
    ambiguous={pair for pair in candidates if sum(a==pair[0] for a,_,_ in candidates)>1 or sum(b==pair[1] for _,b,_ in candidates)>1}
    gene_rows=[]; summaries=[]; diagnostics=[]
    for g,m,basis in candidates:
        x,y=gm[g],mm[m]; cat=''
        if (g,m) in ambiguous: cat='ambiguous_reference_pair'
        elif not x['hash'] and not y['hash']: cat='missing_both_sequence_hash'
        elif not x['hash']: cat='missing_genbank_sequence_hash'
        elif not y['hash']: cat='missing_mitos2_sequence_hash'
        elif x['hash']==y['hash']: cat='exact_sequence_match' if x['record_hash'] in ('',y['input_hash']) else 'coordinate_fasta_match'
        else: cat='sequence_mismatch'
        rot='no'; offset=''
        if cat=='sequence_mismatch' and allow_rotation and x['path'] and y['path']:
            # Only read after hash mismatch; raw coordinates remain intentionally un-compared.
            def seq(p):
                import gzip
                op=gzip.open if str(p).endswith('.gz') else open
                with op(p,'rt') as h:return ''.join(''.join(z.split()) for z in h if not z.startswith('>')).upper()
            offset=rotation(seq(x['path']),seq(y['path']))
            if offset is not None and offset>=0: cat='rotation_equivalent';rot='yes'
        perform=cat in ('exact_sequence_match','coordinate_fasta_match') or (not strict and cat=='sequence_mismatch')
        base={'genbank_reference_key':g,'mitos2_reference_key':m,'target_species':x['target'] or y['target'],'reference_species':x['species'] or y['species'],'genbank_accession':x['accession'],'mitos2_accession':y['accession'],'genbank_coordinate_fasta':x['path'],'mitos2_coordinate_fasta':y['path'],'genbank_coordinate_sequence_length':x['length'],'mitos2_coordinate_sequence_length':y['length'],'genbank_coordinate_sequence_sha256':x['hash'],'mitos2_coordinate_sequence_sha256':y['hash'],'genbank_record_sequence_sha256':x['record_hash'],'mitos2_input_sequence_sha256':y['input_hash'],'sequence_compatibility_category':cat,'sequence_match':yes(cat in ('exact_sequence_match','coordinate_fasta_match')),'rotation_equivalent':rot,'rotation_offset':offset,'coordinate_comparison_performed':yes(perform)}
        if not perform:
            diagnostics.append({**base,'candidate_matching_basis':basis,'rotation_check_result':rot,'reason_comparison_skipped':'strict sequence identity is required before raw coordinate comparison'})
            summaries.append({**base,**{k:0 for k in SUMMARY_FIELDS if k.startswith('n_')},'all_13_genes_present_genbank':'no','all_13_genes_present_mitos2':'no','all_13_exact_match':'no','reference_comparison_category':'ambiguous_reference_pair' if cat=='ambiguous_reference_pair' else 'sequence_not_comparable','note':'Coordinate comparison skipped.'}); continue
        rows=[gene_row(gb[g],mi[m],z,cat,minor,moderate) for z in GENES];gene_rows.extend(rows); kinds=[r['gene_comparison_category'] for r in rows]
        category='all_13_exact' if all(k=='exact_match' for k in kinds) else ('strand_mismatch' if 'strand_mismatch' in kinds else 'gene_missing' if any(k.startswith('missing') for k in kinds) else 'only_minor_boundary_differences' if all(k in ('exact_match','minor_boundary_difference') for k in kinds) else 'moderate_or_major_differences')
        summaries.append({**base,'n_genes_genbank':len({v(r,'gene') for r in gb[g]}),'n_genes_mitos2':len({v(r,'gene') for r in mi[m]}),'n_exact_genes':kinds.count('exact_match'),'n_minor_difference_genes':kinds.count('minor_boundary_difference'),'n_moderate_difference_genes':kinds.count('moderate_boundary_difference'),'n_major_difference_genes':kinds.count('major_coordinate_difference')+kinds.count('coordinate_match_codon_difference'),'n_missing_genes':sum(k.startswith('missing') for k in kinds),'n_strand_mismatches':kinds.count('strand_mismatch'),'all_13_genes_present_genbank':yes(all(r['genbank_present']=='yes' for r in rows)),'all_13_genes_present_mitos2':yes(all(r['mitos2_present']=='yes' for r in rows)),'all_13_exact_match':yes(category=='all_13_exact'),'reference_comparison_category':category,'note':''})
    if not candidates and fail_no_shared: raise RuntimeError('No candidate shared references; existing outputs were left untouched.')
    if not gene_rows and fail_no_shared: raise RuntimeError('No sequence-compatible reference pairs; existing outputs were left untouched.')
    atomic(output,GENE_FIELDS,gene_rows);atomic(summary_output,SUMMARY_FIELDS,summaries);atomic(mismatch_output,DIAG_FIELDS,diagnostics)
    print(f'GenBank reference groups: {len(gb)}\nMITOS2 reference groups: {len(mi)}\ncandidate shared references: {len(candidates)}\nexact sequence matches: {sum(x["sequence_match"]=="yes" for x in summaries)}\nsequence mismatches: {sum(x["sequence_compatibility_category"]=="sequence_mismatch" for x in summaries)}\nmissing sequence hashes: {sum(x["sequence_compatibility_category"].startswith("missing") for x in summaries)}\ncoordinate comparisons performed: {len({r["genbank_reference_key"] for r in gene_rows})}\ncomparison output: {output}')
def main():
 p=argparse.ArgumentParser();p.add_argument('--config',required=True)
 for x in ('genbank-table','mitos2-table','sample-reference-map','output','reference-summary-output','mismatch-output'):p.add_argument('--'+x)
 p.add_argument('--strict-sequence-match',action=argparse.BooleanOptionalAction,default=None);p.add_argument('--allow-rotation-equivalent',action=argparse.BooleanOptionalAction,default=None);p.add_argument('--fail-on-no-shared-references',action=argparse.BooleanOptionalAction,default=None);a=p.parse_args(); sec=yaml(a.config).get('genbank_mitos2_comparison',{});paths=sec.get('paths',{});s=sec.get('settings',{})
 try: compare(a.genbank_table or paths.get('genbank_reference_codon_table'),a.mitos2_table or paths.get('mitos2_reference_codon_table'),a.output or paths.get('gene_comparison_table'),a.reference_summary_output or paths.get('reference_summary_table'),a.mismatch_output or paths.get('sequence_mismatch_table'),s.get('strict_sequence_match',True) if a.strict_sequence_match is None else a.strict_sequence_match,s.get('allow_rotation_equivalent',False) if a.allow_rotation_equivalent is None else a.allow_rotation_equivalent,s.get('fail_on_no_shared_references',True) if a.fail_on_no_shared_references is None else a.fail_on_no_shared_references,float(s.get('minor_position_jaccard_threshold',.99)),float(s.get('moderate_position_jaccard_threshold',.9)))
 except Exception as e: raise SystemExit(f'GenBank-versus-MITOS2 comparison failed: {e}')
if __name__=='__main__':main()
