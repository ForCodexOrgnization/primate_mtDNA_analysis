#!/usr/bin/env python3
"""Report native-coordinate clustering of source-classified heteroplasmies."""
import argparse,csv,random,sys
from collections import Counter,defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from qc_analysis.lib.simple_yaml import read_simple_yaml

def max_cluster(pos,length,window):
    if not pos:return 0
    p=sorted(pos); q=p+[x+length for x in p]; j=0; best=0
    for i,x in enumerate(q[:len(p)]):
        while j<len(q) and q[j]-x<=window:j+=1
        best=max(best,j-i)
    return best
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--config',type=Path,required=True);a=ap.parse_args();cfg=read_simple_yaml(a.config);sec=cfg.get('local_heteroplasmy_qc',{})
    if sec.get('enabled',True) is False:return 0
    source=ROOT/str(sec.get('source_report','results/qc/pre_liftover_variant_qc/reports/source_variant_qc.tsv'));out=ROOT/str(sec.get('output_dir','results/qc/local_heteroplasmy_qc'))/'reports';out.mkdir(parents=True,exist_ok=True)
    with source.open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
    eligible=[r for r in rows if r['source_call_class']=='HET' and r['source_filter']=='PASS' and float(r['source_dp'])>=float(sec.get('dp_min',100)) and r['source_variant_qc']=='PASS']
    by=defaultdict(list)
    for r in eligible:by[r['sample']].append(int(r['source_pos']))
    recurrence=Counter(int(r['source_pos']) for r in eligible); population=list(recurrence);weights=list(recurrence.values());length=int(sec.get('mt_length',16569));window=int(sec.get('window_bp',250));sim=int(sec.get('simulations',2000));minimum=int(sec.get('min_het_variants',11));rng=random.Random(int(sec.get('random_seed',1729))); summaries=[]; details=[]
    for sample,pos in sorted(by.items()):
        observed=max_cluster(pos,length,window);assessed=len(pos)>=minimum
        null=[max_cluster(rng.choices(population,weights=weights,k=len(pos)),length,window) for _ in range(sim)] if assessed and population else []
        p=(1+sum(x>=observed for x in null))/(len(null)+1) if null else None; clustered=bool(p is not None and p<=float(sec.get('empirical_p_max',.01)))
        summaries.append({'sample':sample,'n_het':len(pos),'assessed':str(assessed).lower(),'max_variants_in_250bp_window':observed,'empirical_p':'NA' if p is None else f'{p:.6g}','locally_clustered_heteroplasmy':str(clustered).lower()})
        for x in pos:details.append({'sample':sample,'source_pos':x,'source_call_class':'HET','eligible':True,'locally_clustered_heteroplasmy':str(clustered).lower()})
    for name,fields,data in [('local_heteroplasmy_sample_summary.tsv',['sample','n_het','assessed','max_variants_in_250bp_window','empirical_p','locally_clustered_heteroplasmy'],summaries),('local_heteroplasmy_variant_detail.tsv',['sample','source_pos','source_call_class','eligible','locally_clustered_heteroplasmy'],details)]:
        with (out/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(data)
    return 0
if __name__=='__main__':raise SystemExit(main())
