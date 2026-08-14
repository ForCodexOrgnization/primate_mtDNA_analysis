#!/usr/bin/env python3
"""Build one tRNAscan position index per unique mitochondrial reference."""
import argparse, csv, sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from qc_analysis.lib.simple_yaml import read_simple_yaml
from qc_analysis.lib.reference_utils import normalized_fasta_sequence_sha256
from qc_analysis.scripts.build_trna_position_index import build, SUPPORTED_CHROM_NORMALIZATIONS
from qc_analysis.lib.trnascan_utils import read_fasta, validate_trna_index

def table(path):
    with open(path,newline="") as h: return list(csv.DictReader(h,delimiter="\t"))

def resolve_fastas(manifest, min_length=14000, max_length=19000):
    """Resolve only mitochondrial FASTAs, never a whole-genome fallback."""
    result={}
    for r in table(manifest):
        key=r.get("reference_key") or r.get("reference_id") or r.get("accession")
        path=next((r.get(k) for k in ("chrM_fasta_path","chrM_expected_output_fasta","fasta","fasta_path") if r.get(k)),None)
        if not key or not path: continue
        seqs=read_fasta(path); target=next((r.get(k) for k in ("target_sequence_id","chrM_record_id","sequence_id") if r.get(k)),None)
        if len(seqs)>1 and not target:
            raise ValueError(f"Reference {key}: multi-contig FASTA requires target_sequence_id: {path}")
        if target and target not in seqs: raise ValueError(f"Reference {key}: target sequence {target!r} absent from {path}")
        selected=target or next(iter(seqs)); length=len(seqs[selected])
        if not min_length <= length <= max_length:
            raise ValueError(f"Reference {key}: mitochondrial length {length} outside [{min_length}, {max_length}]")
        result[key]={"path":path,"target_sequence_id":target}
    return result

def add_coordinate_fastas(result, samples, min_length=14000, max_length=19000):
    """Resolve and validate aliases of each variant-calling reference FASTA.

    Values from the coordinate map deliberately replace legacy manifest values.
    A reference key identifies normalized sequence content, not a physical path.
    Consequently, multiple paths are accepted only after their content hashes
    have been shown to be identical.
    """
    grouped={}
    for row in samples:
        key=(row.get("reference_key") or "").strip()
        if not key:
            raise ValueError("Coordinate reference map row is missing reference_key")
        grouped.setdefault(key,[]).append(row)
    for key,rows in grouped.items():
        paths={(row.get("coordinate_reference_fasta") or "").strip() for row in rows}
        if "" in paths:
            raise ValueError(f"Reference {key}: FAIL_MISSING_FASTA: missing coordinate_reference_fasta")
        paths=sorted(paths)
        hashes={}
        lengths={}
        for path in paths:
            if not Path(path).is_file():
                raise ValueError(f"Reference {key}: FAIL_MISSING_FASTA: coordinate FASTA does not exist: {path}")
            # Enforce the coordinate-FASTA contract independently of hashing:
            # the shared hash helper deliberately concatenates all FASTA records.
            seqs=read_fasta(path)
            if len(seqs) != 1:
                raise ValueError(f"Reference {key}: coordinate FASTA must contain exactly one record: {path}")
            info=normalized_fasta_sequence_sha256(path)
            length=info["sequence_length"]
            if not min_length <= length <= max_length:
                raise ValueError(f"Reference {key}: mitochondrial length {length} outside [{min_length}, {max_length}]: {path}")
            hashes[path]=info["sequence_sha256"]
            lengths[path]=length
        observed=set(hashes.values())
        if len(observed) != 1:
            details={path:hashes[path] for path in paths}
            raise ValueError(f"Reference {key}: FAIL_HASH_CONFLICT: coordinate FASTAs have different normalized sequence hashes: {details}")
        actual=next(iter(observed))
        encoded=key[6:].lower() if key.startswith("mtref_") and len(key)==70 else None
        if encoded and all(c in "0123456789abcdef" for c in encoded) and actual != encoded:
            raise ValueError(f"Reference {key}: FAIL_REFERENCE_KEY_HASH_MISMATCH: FASTA normalized SHA256 is {actual}")
        declared={(row.get("coordinate_reference_sequence_sha256") or "").strip().lower() for row in rows}
        declared.discard("")
        if len(declared) > 1:
            raise ValueError(f"Reference {key}: conflicting coordinate_reference_sequence_sha256 values")
        if declared and actual != next(iter(declared)):
            raise ValueError(f"Reference {key}: coordinate FASTA SHA256 mismatch: expected {next(iter(declared))}, observed {actual}")
        canonical=paths[0]
        result[key]={"path":canonical,"target_sequence_id":None,
                     "canonical_fasta":canonical,"n_fasta_aliases":len(paths),
                     "fasta_aliases":";".join(paths),
                     "normalized_sequence_sha256":actual,
                     "sequence_length":lengths[canonical],
                     "validation_status":("VALID_SINGLE_FASTA" if len(paths)==1
                                          else "VALID_IDENTICAL_FASTA_ALIASES")}
    return result

def valid(path, reference_key):
    p=Path(path)
    if not p.is_file() or p.stat().st_size==0:return False
    try:
        validate_trna_index(p, reference_key); return True
    except (OSError, EOFError, ValueError):return False

def trna_chrom_normalization(reference_key, settings):
    key = "human_trna_chrom_norm" if reference_key == "human" else "species_trna_chrom_norm"
    configured = settings.get(key)
    value = str(configured or "none").strip().lower()
    if value not in SUPPORTED_CHROM_NORMALIZATIONS:
        raise ValueError(f"Unsupported chromosome normalization for {key}: {configured!r}")
    return value

def write(path,rows,columns):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",newline="") as h:
        w=csv.DictWriter(h,columns,delimiter="\t",extrasaction="ignore");w.writeheader();w.writerows(rows)

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--config",default="config/qc_preprocessing.yaml");ap.add_argument("--workers",type=int,default=1)
    ap.add_argument("--overwrite",action="store_true");ap.add_argument("--reference-key");ap.add_argument("--task-manifest");a=ap.parse_args()
    cfg=read_simple_yaml(Path(a.config));sec=cfg["trna_match"];p,s=sec["paths"],sec["settings"]
    samples=table(p["sample_reference_map"]); minimum=int(s.get("min_mitochondrial_reference_length",14000)); maximum=int(s.get("max_mitochondrial_reference_length",19000))
    # The sample coordinate map is authoritative.  The old manifest remains an
    # optional source only for legacy configurations and is overridden below.
    fastas={}
    if p.get("reference_fasta_manifest"):
        fastas.update(resolve_fastas(p["reference_fasta_manifest"],minimum,maximum))
    fastas=add_coordinate_fastas(fastas,samples,minimum,maximum)
    keys=sorted({r["reference_key"] for r in samples if r.get("reference_key")})
    if a.reference_key: keys=[a.reference_key]
    outdir=Path(p["reference_trna_index_dir"]); reportdir=Path(p["index_build_reports_dir"]); scan=Path(p["trnascan_output_dir"])
    template=p["reference_trna_index_template"]
    tasks=[]
    human_key="human"; human_fasta=p["human_fasta"]; human_output=p["human_trna_index"]
    if not a.reference_key: tasks.append((human_key,{"path":human_fasta,"target_sequence_id":s.get("human_target_sequence_id")},human_output))
    for key in keys:
        if key not in fastas: tasks.append((key,None,str(template).format(reference_trna_index_dir=outdir,reference_key=key)))
        else: tasks.append((key,fastas[key],str(template).format(reference_trna_index_dir=outdir,reference_key=key)))
    if a.task_manifest:
        write(a.task_manifest,[{"task_id":i,"reference_key":k,"fasta":(f or {}).get("path", ""),
              "canonical_fasta":(f or {}).get("canonical_fasta",(f or {}).get("path", "")),
              "n_fasta_aliases":(f or {}).get("n_fasta_aliases",1 if f else 0),
              "fasta_aliases":(f or {}).get("fasta_aliases",(f or {}).get("path", "")),
              "normalized_sequence_sha256":(f or {}).get("normalized_sequence_sha256",""),
              "sequence_length":(f or {}).get("sequence_length",""),
              "validation_status":(f or {}).get("validation_status",""),"output_index":o}
              for i,(k,f,o) in enumerate(tasks,1)],
              ["task_id","reference_key","fasta","canonical_fasta","n_fasta_aliases","fasta_aliases",
               "normalized_sequence_sha256","sequence_length","validation_status","output_index"]);return
    def one(task):
        key,fasta_info,output=task; fasta=fasta_info["path"] if fasta_info else None
        if not fasta:return {"reference_key":key,"fasta":"","output_index":output,"status":"failed","notes":"FASTA not resolved"}
        provenance={name:fasta_info.get(name,"") for name in
                    ("canonical_fasta","n_fasta_aliases","fasta_aliases","normalized_sequence_sha256",
                     "sequence_length","validation_status")}
        if valid(output,key) and not a.overwrite:return {"reference_key":key,"fasta":fasta,"output_index":output,"status":"skipped","notes":"valid existing index",**provenance}
        ns=SimpleNamespace(reference_key=key,fasta=fasta,trnascan_out=None,trnascan_ss=None,run_trnascan=True,
          trnascan_bin=s.get("trnascan_bin","tRNAscan-SE"),trnascan_mode=s.get("trnascan_mode","mito_mammal"),threads=s.get("trnascan_threads",1),
          trnascan_extra_args=s.get("trnascan_extra_args","") or "",trnascan_prefix=str(scan/key),output=output,overwrite=True,summary=None,
          chrom_normalization=trna_chrom_normalization(key,s),
          target_sequence_id=fasta_info.get("target_sequence_id"),allow_ss_order_fallback=bool(s.get("allow_ss_order_fallback",False)),
          max_sequence_mismatch_rate=s.get("max_sequence_mismatch_rate",0.0))
        try:return {**build(ns),**provenance}
        except Exception as e:return {"reference_key":key,"fasta":fasta,"output_index":output,"status":"failed","notes":str(e),**provenance}
    results=[]
    with ThreadPoolExecutor(max_workers=max(1,a.workers)) as pool:
        futures={pool.submit(one,t):t[0] for t in tasks}
        for future in as_completed(futures):results.append(future.result())
    results.sort(key=lambda r:r["reference_key"])
    columns=["reference_key","fasta","canonical_fasta","n_fasta_aliases","fasta_aliases","normalized_sequence_sha256","sequence_length","validation_status","selected_record_id","fasta_length","fasta_sha256","trnascan_out","trnascan_ss","output_index","n_trna_records","n_index_rows","n_positive_strand_trna","n_negative_strand_trna","n_stem_positions","n_loop_positions","n_missing_structure","n_fasta_sequence_mismatch","status","notes"]
    write(reportdir/"trna_index_build_summary.tsv",results,columns);write(reportdir/"trna_index_build_failures.tsv",[r for r in results if r["status"]=="failed"],columns)
    bykey={r["reference_key"]:r for r in results}; maps=[]
    for row in samples:
        r=bykey.get(row.get("reference_key"),{});maps.append({"sample":row.get("sample",row.get("sample_name","")),"reference_key":row.get("reference_key",""),"trna_index":r.get("output_index",""),"status":r.get("status","failed"),"notes":r.get("notes","")})
    write(reportdir/"sample_trna_index_map.tsv",maps,["sample","reference_key","trna_index","status","notes"])
    if any(r["status"]=="failed" for r in results):raise SystemExit("One or more tRNA indexes failed; see failure report")
if __name__=="__main__":main()
