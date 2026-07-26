#!/usr/bin/env python3
"""Compare source tRNA genes with best-overlapping human genes after liftover."""
import argparse,csv,gzip,sys
from collections import defaultdict
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from qc_analysis.lib.match_utils import load_coordinate_map

def read(path):
    op=gzip.open if str(path).endswith('.gz') else open
    with op(path,'rt',newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
def genes(path):
    grouped=defaultdict(list)
    for r in read(path):grouped[r['trna_id']].append(r)
    return grouped
def main():
    p=argparse.ArgumentParser();p.add_argument('--source-index',required=True);p.add_argument('--human-index',required=True);p.add_argument('--coordinate-map',required=True);p.add_argument('--output',required=True);p.add_argument('--min-overlap-ratio',type=float,default=.8);a=p.parse_args()
    source,human=genes(a.source_index),genes(a.human_index); cmap=load_coordinate_map(a.coordinate_map)
    hpos={k:{int(r['pos']) for r in rs} for k,rs in human.items()};rows=[]
    for sid,srs in source.items():
        mapped={int(cmap[p]) for p in (int(r['pos']) for r in srs) if cmap.get(p) not in ('',None,'.')}
        best=max(human,key=lambda h:len(mapped&hpos[h])) if mapped and human else None; overlap=len(mapped&hpos[best]) if best else 0
        sr=srs[0];hr=human[best][0] if best else {};ratio=overlap/len(srs) if srs else 0
        orient=str(sr.get('strand'))==str(hr.get('strand'));iso=sr.get('aa')==hr.get('aa');anti=sr.get('anticodon')==hr.get('anticodon')
        rows.append({'source_trna_id':sid,'human_trna_id':best or '.','mapped_human_begin':min(mapped) if mapped else '.','mapped_human_end':max(mapped) if mapped else '.','overlap_bp':overlap,'overlap_ratio':f'{ratio:.6f}','orientation_match':'yes' if orient else 'no','isotype_match':'yes' if iso else 'no','anticodon_match':'yes' if anti else 'no','status':'pass' if ratio>=a.min_overlap_ratio and orient and iso and anti else 'fail'})
    Path(a.output).parent.mkdir(parents=True,exist_ok=True);cols=list(rows[0]) if rows else ['source_trna_id','human_trna_id','mapped_human_begin','mapped_human_end','overlap_bp','overlap_ratio','orientation_match','isotype_match','anticodon_match','status']
    with open(a.output,'w',newline='') as h:w=csv.DictWriter(h,cols,delimiter='\t');w.writeheader();w.writerows(rows)
if __name__=='__main__':main()
