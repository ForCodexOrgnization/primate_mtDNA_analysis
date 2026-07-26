#!/usr/bin/env python3
"""Build one tRNAscan position index per unique mitochondrial reference."""
import argparse, csv, sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from qc_analysis.lib.simple_yaml import read_simple_yaml
from qc_analysis.scripts.build_trna_position_index import build

def table(path):
    with open(path,newline="") as h: return list(csv.DictReader(h,delimiter="\t"))

def resolve_fastas(manifest):
    result={}
    for r in table(manifest):
        key=r.get("reference_key") or r.get("reference_id") or r.get("accession")
        path=next((r.get(k) for k in ("fasta","fasta_path","chrM_fasta_path","chrM_expected_output_fasta","wg_expected_output_fasta") if r.get(k)),None)
        if key and path: result[key]=path
    return result

def valid(path):
    import gzip
    p=Path(path)
    if not p.is_file() or p.stat().st_size==0:return False
    op=gzip.open if str(p).endswith(".gz") else open
    try:
        with op(p,"rt") as h: return "index_format_version" in h.readline().split("\t")
    except OSError:return False

def write(path,rows,columns):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",newline="") as h:
        w=csv.DictWriter(h,columns,delimiter="\t",extrasaction="ignore");w.writeheader();w.writerows(rows)

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--config",default="config/qc_preprocessing.yaml");ap.add_argument("--workers",type=int,default=1)
    ap.add_argument("--overwrite",action="store_true");ap.add_argument("--reference-key");ap.add_argument("--task-manifest");a=ap.parse_args()
    cfg=read_simple_yaml(Path(a.config));sec=cfg["trna_match"];p,s=sec["paths"],sec["settings"]
    samples=table(p["sample_reference_map"]); fastas=resolve_fastas(p["reference_fasta_manifest"])
    keys=sorted({r["reference_key"] for r in samples if r.get("reference_key")})
    if a.reference_key: keys=[a.reference_key]
    outdir=Path(p["reference_trna_index_dir"]); reportdir=Path(p["index_build_reports_dir"]); scan=Path(p["trnascan_output_dir"])
    template=p["reference_trna_index_template"]
    tasks=[]
    human_key="human"; human_fasta=p["human_fasta"]; human_output=p["human_trna_index"]
    if not a.reference_key: tasks.append((human_key,human_fasta,human_output))
    for key in keys:
        if key not in fastas: tasks.append((key,None,str(template).format(reference_trna_index_dir=outdir,reference_key=key)))
        else: tasks.append((key,fastas[key],str(template).format(reference_trna_index_dir=outdir,reference_key=key)))
    if a.task_manifest:
        write(a.task_manifest,[{"task_id":i,"reference_key":k,"fasta":f or "","output_index":o} for i,(k,f,o) in enumerate(tasks,1)],
              ["task_id","reference_key","fasta","output_index"]);return
    def one(task):
        key,fasta,output=task
        if not fasta:return {"reference_key":key,"fasta":"","output_index":output,"status":"failed","notes":"FASTA not resolved"}
        if valid(output) and not a.overwrite:return {"reference_key":key,"fasta":fasta,"output_index":output,"status":"skipped","notes":"valid existing index"}
        ns=SimpleNamespace(reference_key=key,fasta=fasta,trnascan_out=None,trnascan_ss=None,run_trnascan=True,
          trnascan_bin=s.get("trnascan_bin","tRNAscan-SE"),trnascan_mode=s.get("trnascan_mode","mito_mammal"),threads=s.get("trnascan_threads",1),
          trnascan_extra_args=s.get("trnascan_extra_args","") or "",trnascan_prefix=str(scan/key),output=output,overwrite=True,summary=None,
          chrom_normalization=s.get("species_trna_chrom_norm","none"),max_sequence_mismatch_rate=s.get("max_sequence_mismatch_rate",0.0))
        try:return build(ns)
        except Exception as e:return {"reference_key":key,"fasta":fasta,"output_index":output,"status":"failed","notes":str(e)}
    results=[]
    with ThreadPoolExecutor(max_workers=max(1,a.workers)) as pool:
        futures={pool.submit(one,t):t[0] for t in tasks}
        for future in as_completed(futures):results.append(future.result())
    results.sort(key=lambda r:r["reference_key"])
    columns=["reference_key","fasta","fasta_length","trnascan_out","trnascan_ss","output_index","n_trna_records","n_index_rows","n_positive_strand_trna","n_negative_strand_trna","n_stem_positions","n_loop_positions","n_missing_structure","n_fasta_sequence_mismatch","status","notes"]
    write(reportdir/"trna_index_build_summary.tsv",results,columns);write(reportdir/"trna_index_build_failures.tsv",[r for r in results if r["status"]=="failed"],columns)
    bykey={r["reference_key"]:r for r in results}; maps=[]
    for row in samples:
        r=bykey.get(row.get("reference_key"),{});maps.append({"sample":row.get("sample",row.get("sample_name","")),"reference_key":row.get("reference_key",""),"trna_index":r.get("output_index",""),"status":r.get("status","failed"),"notes":r.get("notes","")})
    write(reportdir/"sample_trna_index_map.tsv",maps,["sample","reference_key","trna_index","status","notes"])
    if any(r["status"]=="failed" for r in results):raise SystemExit("One or more tRNA indexes failed; see failure report")
if __name__=="__main__":main()
