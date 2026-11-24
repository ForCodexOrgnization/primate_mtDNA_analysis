#!/usr/bin/env python3
"""Audit species-union vs sample-specific NUMT overlap from existing per-variant reports.

Does NOT modify filtering. Requires the per-sample upstream besthit and highconf BED
files already used by qc_analysis. All coordinates are native mtDNA coordinates.

The default --mt-length 16569 is a COMPARISON PROXY, not an assertion that every
primate reference has that length. Provide --lengths-tsv for exact fractions.
"""
from __future__ import annotations
import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if not (ROOT / 'qc_analysis/lib').exists():
    ROOT = Path.cwd()
sys.path.insert(0, str(ROOT))
from qc_analysis.lib.numt_annotation import load_sample_numts
from qc_analysis.lib.simple_yaml import read_simple_yaml

COHORTS = {
    'hom10': ('config/qc_cohort_hom10_mean_mtcn.yaml', 'results/qc_cohort_hom10_mean_mtcn'),
    'homhet40': ('config/qc_cohort_homhet40_mean_mtcn.yaml', 'results/qc_cohort_homhet40_mean_mtcn'),
}

def resolve(v):
    p = Path(v).expanduser()
    return p if p.is_absolute() else ROOT / p

def rows(path):
    with path.open(newline='', encoding='utf8') as f:
        yield from csv.DictReader(f, delimiter='\t')

def union_length(intervals, size):
    spans = sorted((max(1, min(a,b)), min(size, max(a,b))) for a,b in intervals)
    total, last_a, last_b = 0, None, None
    for a,b in spans:
        if a>b: continue
        if last_a is None: last_a,last_b = a,b
        elif a <= last_b+1: last_b=max(last_b,b)
        else: total+=last_b-last_a+1; last_a,last_b=a,b
    if last_a is not None: total+=last_b-last_a+1
    return total

def interval_list(catalog, tier):
    result=[]
    for t, vals in catalog.items():
        if tier=='any' or t=='HIGH_CONF_NUMT': result.extend(vals)
    return result

def in_spans(pos, spans):
    return any(a<=pos<=b for a,b in spans)

def read_lengths(path):
    if not path:return {}
    result={}
    for r in rows(path):
        k=(r['species'],r['reference_key'])
        n=int(r['mt_length'])
        if n<1000:raise ValueError(f'Invalid mt length for {k}: {n}')
        result[k]=n
    return result

def load_species_catalog(path):
    d=defaultdict(lambda:defaultdict(list))
    for r in rows(path):
        k=(r['species'],r['reference_key'])
        try: a,b=int(r['chrm_start']),int(r['chrm_end'])
        except (KeyError,ValueError):continue
        d[k][r['tier']].append((min(a,b),max(a,b)))
    return d

def names_from_config(sec):
    best=resolve(sec.get('numt_besthit_dir','data/numt_besthit'))
    high=resolve(sec.get('numt_highconf_dir', str(best)))
    # The committed config may contain an old HPC absolute path; use this checkout's data.
    fallback=ROOT/'data/numt_besthit'
    if not best.is_dir() and fallback.is_dir():best=fallback
    if not high.is_dir() and fallback.is_dir():high=fallback
    if not best.is_dir():
        raise FileNotFoundError(f'NUMT besthit directory not found: {best}')
    # Same priority as run_local_heteroplasmy_qc.py defaults
    best_pats=sec.get('numt_besthit_pattern') or [
        '{sample}.numt_vs_chrM.besthit.tsv','{sample}.numt_besthit.tsv',
        '{sample}.besthit.tsv','{sample}.tsv']
    high_pats=sec.get('numt_highconf_pattern') or ['{sample}.highconf_numt.bed']
    def clean(x):
        if isinstance(x,str):return [z.strip() for z in x.split(',') if z.strip()]
        return list(x)
    return best,high,clean(best_pats),clean(high_pats)

def find_file(directory,patterns,sample):
    for pat in patterns:
        p=directory / pat.replace('{sample}',sample)
        if p.is_file():return p
    return None

def emit(path, header, records):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,delimiter='\t',fieldnames=header)
        w.writeheader()
        for record in records:w.writerow(record)


def analyze(cohort, args, lengths):
    cfg_rel, root_rel=COHORTS[cohort]
    cfg=read_simple_yaml(resolve(cfg_rel))
    sec=cfg.get('local_heteroplasmy_qc') or {}
    root=resolve(root_rel)
    variants=resolve(args.variant_dir)/f'{cohort}_numt_overlap_per_variant.tsv'
    catalogue=root/'local_heteroplasmy_qc/reports/species_numt_intervals.tsv'
    if not variants.exists():raise FileNotFoundError(variants)
    if not catalogue.exists():raise FileNotFoundError(catalogue)
    species_catalog=load_species_catalog(catalogue)
    best_dir,high_dir,best_pats,high_pats=names_from_config(sec)
    # Prefer a directory with matching sample files, as in production NUMT QC.
    from qc_analysis.scripts.run_local_heteroplasmy_qc import choose_numt_dir
    sample_ids = sorted({r['sample'] for r in rows(variants)})
    best_dir = choose_numt_dir(sec.get('numt_besthit_dir'),best_pats,sample_ids,
                               [ROOT/'data/numt_besthit'])
    high_dir = choose_numt_dir(sec.get('numt_highconf_dir'),high_pats,sample_ids,
                               [best_dir,ROOT/'data/numt_besthit'])
    print(cohort,'NUMT besthit directory:',best_dir,'highconf directory:',high_dir,
          'variant sample count:',len(sample_ids))
    # IMPORTANT: no sample without a besthit file is counted as 0 overlap.
    sample_cache={}
    input_rows=0
    eligible_rows=0
    counts=defaultdict(Counter)
    per_species=defaultdict(Counter)
    missing=Counter()
    sample_meta={}
    for row in rows(variants):
        input_rows+=1
        sample=row['sample']; species=row['species']; reference=row['reference_key']
        spec_key=(species,reference)
        if sample not in sample_cache:
            bh=find_file(best_dir,best_pats,sample)
            hc=find_file(high_dir,high_pats,sample)
            if not bh:
                sample_cache[sample]=None
            else:
                # A missing highconf BED does not invalidate best-hit overlaps.
                # Preserve missingness instead of converting unknown to 0%.
                ivals=load_sample_numts(sample,species,reference,bh,
                                        hc if hc else high_dir / f'{sample}.highconf_numt.bed')
                sample_cache[sample]={
                    'any':[(v.chrm_start,v.chrm_end) for v in ivals],
                    'high':[(v.chrm_start,v.chrm_end) for v in ivals if v.tier=='HIGH_CONF_NUMT']
                             if hc else None}
                sample_meta[sample]=spec_key
        sample_spans=sample_cache[sample]
        if sample_spans is None:
            missing['missing_besthit_file']+=1
            continue
        if spec_key not in species_catalog:
            missing['missing_species_catalog']+=1
            continue
        try:pos=int(row['source_pos'])
        except ValueError:
            missing['invalid_source_pos']+=1
            continue
        eligible_rows+=1
        stage='AFTER' if row['final_status']=='PASS' else 'REMOVED'
        cls=row['call_class']
        for kind in ('any','high'):
            if sample_spans[kind] is None:
                missing['missing_highconf_file']+=1
                continue
            sample_hit=in_spans(pos,sample_spans[kind])
            species_hit=in_spans(pos,interval_list(species_catalog[spec_key],kind))
            old_col='numt_besthit_overlap' if kind=='any' else 'numt_highconf_overlap'
            if row.get(old_col) in ('0','1') and species_hit != (row[old_col]=='1'):
                missing[f'species_catalog_mismatch_{kind}']+=1
            if sample_hit and not species_hit:
                missing['sample_hit_not_in_species_union']+=1
            for label,status in ((cls,stage),('ALL',stage),(cls,'BEFORE'),('ALL','BEFORE')):
                c=counts[(kind,label,status)]
                c['n']+=1;c['sample_hits']+=int(sample_hit);c['species_hits']+=int(species_hit)
            per_species[(kind,species,reference,cls,stage)]['n']+=1
            per_species[(kind,species,reference,cls,stage)]['sample_hits']+=int(sample_hit)
            per_species[(kind,species,reference,cls,stage)]['species_hits']+=int(species_hit)
    result=[]
    for (kind,cls,stage),c in sorted(counts.items()):
        n=c['n'];result.append(dict(cohort=cohort,interval_tier=kind,call_class=cls,stage=stage,
                                  variants=n,sample_hit=c['sample_hits'],sample_pct=round(100*c['sample_hits']/n,4),
                                  species_union_hit=c['species_hits'],species_union_pct=round(100*c['species_hits']/n,4),
                                  reduction_pp=round(100*(c['species_hits']-c['sample_hits'])/n,4)))
    if not result:
        print(f'WARNING {cohort}: no comparable sample-specific overlap rows; check '
              f'{args.output / (cohort + "_missing_numt_inputs.tsv")} and NUMT input directory',
              file=sys.stderr)
    emit(args.output/f'{cohort}_sample_vs_species_numt_overlap.tsv',
         ['cohort','interval_tier','call_class','stage','variants','sample_hit','sample_pct','species_union_hit','species_union_pct','reduction_pp'],result)
    coverage=[]
    for (species,ref), catalog in sorted(species_catalog.items()):
        n=lengths.get((species,ref),args.mt_length)
        for kind in ('any','high'):
            bp=union_length(interval_list(catalog,kind),n)
            coverage.append(dict(cohort=cohort,species=species,reference_key=ref,interval_tier=kind,
                                 mt_length=n,length_source=('species_specific' if (species,ref) in lengths else 'fixed_proxy'),
                                 numt_union_bp=bp,numt_union_pct=round(bp/n*100,4)))
    emit(args.output/f'{cohort}_species_numt_coverage.tsv',
         ['cohort','species','reference_key','interval_tier','mt_length','length_source','numt_union_bp','numt_union_pct'],coverage)
    per_sample=[]
    for sample,spans in sorted(sample_cache.items()):
        if spans is None:continue
        spec_key=sample_meta[sample]; n=lengths.get(spec_key,args.mt_length)
        for kind in ('any','high'):
            if spans[kind] is None:continue
            bp=union_length(spans[kind],n)
            per_sample.append(dict(cohort=cohort,sample=sample,species=spec_key[0],reference_key=spec_key[1],
                                   interval_tier=kind,mt_length=n,
                                   length_source=('species_specific' if spec_key in lengths else 'fixed_proxy'),
                                   numt_bp=bp,numt_pct=round(bp/n*100,4)))
    emit(args.output/f'{cohort}_sample_numt_coverage.tsv',
         ['cohort','sample','species','reference_key','interval_tier','mt_length','length_source','numt_bp','numt_pct'],per_sample)
    emit(args.output/f'{cohort}_missing_numt_inputs.tsv',['reason','records'],
         [{'reason':k,'records':v} for k,v in sorted(missing.items())])
    print(cohort,'source rows=',input_rows,'assessable besthit rows=',eligible_rows,
          'samples checked=',len(sample_cache),'missing data=',dict(missing))
    if input_rows == 0:
        print('ERROR: upstream per-variant file contains no data rows:',variants,file=sys.stderr)
    elif eligible_rows == 0:
        print('ERROR: all rows dropped before overlap assessment. Check missing_numt_inputs.tsv '
              'and besthit directory/pattern in config.',file=sys.stderr)
    for r in result:
        if r['call_class'] in ('HOM','HET') and r['stage'] in ('BEFORE','AFTER'):
            print(' ',r['interval_tier'],r['call_class'],r['stage'],
                  'sample:',r['sample_pct'],'species:',r['species_union_pct'])
    return result,coverage

def plot(results,coverage,out):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib not available; TSV outputs were written',file=sys.stderr)
        return
    for kind in ('any','high'):
        entries=[r for r in results if r['interval_tier']==kind and r['call_class'] in ('HOM','HET')
                 and r['stage'] in ('BEFORE','AFTER')]
        if not entries:continue
        labels=[f"{x['cohort']} {x['call_class']}\n{x['stage']}" for x in entries]
        xs=list(range(len(entries)))
        width=.38
        fig,ax=plt.subplots(figsize=(max(9,len(xs)*1.1),5.5))
        ax.bar([x-width/2 for x in xs],[r['species_union_pct'] for r in entries],width,label='Species-level union')
        ax.bar([x+width/2 for x in xs],[r['sample_pct'] for r in entries],width,label='Sample-specific')
        ax.set_ylabel('Variant overlap (%)');ax.set_ylim(0,105);ax.set_xticks(xs,labels)
        ax.legend(frameon=False);ax.set_title(f'NUMT overlap before vs after — {kind}')
        fig.tight_layout();fig.savefig(out/f'numt_{kind}_sample_vs_species.png',dpi=240);plt.close(fig)
    for kind in ('any','high'):
        vals=[r['numt_union_pct'] for r in coverage if r['interval_tier']==kind]
        if not vals:continue
        fig,ax=plt.subplots(figsize=(7,4.5));ax.hist(vals,bins=20)
        ax.set_xlabel('Fraction of mtDNA bases covered by species NUMT intervals (%)')
        ax.set_ylabel('Species/reference catalogues');ax.set_title(f'Species union NUMT coverage — {kind}')
        fig.tight_layout();fig.savefig(out/f'species_numt_coverage_{kind}.png',dpi=240);plt.close(fig)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--cohort',choices=('hom10','homhet40','both'),default='both')
    p.add_argument('--variant-dir',type=Path,default=ROOT/'results/numt_overlap_before_after')
    p.add_argument('--output',type=Path,default=ROOT/'results/numt_union_audit')
    p.add_argument('--mt-length',type=int,default=16569,help='Proxy only; use --lengths-tsv for species-reference-specific lengths')
    p.add_argument('--lengths-tsv',type=Path,help='TSV: species, reference_key, mt_length')
    args=p.parse_args();args.variant_dir=resolve(args.variant_dir);args.output=resolve(args.output)
    args.output.mkdir(parents=True,exist_ok=True)
    lens=read_lengths(args.lengths_tsv)
    results=[];coverage=[]
    for c in (COHORTS if args.cohort=='both' else (args.cohort,)):
        r,x=analyze(c,args,lens);results.extend(r);coverage.extend(x)
    plot(results,coverage,args.output)
    print('OUTPUT',args.output)
    if not lens:print('WARNING: coverage percentages use fixed 16569-bp denominator; provide --lengths-tsv for exact species-specific results')

if __name__=='__main__':main()