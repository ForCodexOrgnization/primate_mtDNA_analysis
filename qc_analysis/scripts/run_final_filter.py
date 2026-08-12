#!/usr/bin/env python3
"""Strict terminal sample/variant filtering of the most downstream VCF."""
from __future__ import annotations
import argparse,csv,gzip,shutil,subprocess,sys,tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from qc_analysis.lib.simple_yaml import read_simple_yaml
SAMPLE_COLUMNS="sample species intraspecies_status human_contamination_status interspecies_status sample_level_qc_status final_sample_status final_sample_fail_reasons final_sample_warnings vcf_source".split()
VARIANT_COLUMNS="sample original_chrom original_pos original_ref original_alt human_chrom human_pos human_ref human_alt liftover_status human_contamination_status interspecies_status sample_variant_qc_status match_status final_variant_status final_variant_fail_reasons".split()
def resolve(v):
 p=Path(str(v)).expanduser();return p if p.is_absolute() else ROOT/p
def read_tsv(p):
 with p.open(newline="",encoding="utf-8") as h:return list(csv.DictReader(h,delimiter="\t"))
def pick(row,names,default="NOT_AVAILABLE"):
 if isinstance(names,str):names=[x.strip() for x in names.split(",") if x.strip()]
 low={k.lower():v for k,v in row.items()};return next((low[x.lower()] for x in names if low.get(x.lower(),"")!=""),default)
def index(p):return {pick(r,["sample","Sample"]):r for r in read_tsv(p)}
def names(v):return [x.strip() for x in v.split(",")] if isinstance(v,str) else list(v or [])
def is_fail(v,configured):return v.strip().lower() in {str(x).strip().lower() for x in names(configured)}
def sample_name(p):return p.name.split(".lifted",1)[0].split(".vcf",1)[0]
def open_vcf(p):return gzip.open(p,"rt") if p.suffix==".gz" else p.open(encoding="utf-8")
def find_vcf(directory,sample):
 candidates=sorted(set(directory.glob(f"{sample}*.vcf"))|set(directory.glob(f"{sample}*.vcf.gz")))
 return candidates[0] if candidates else None
def bgzip_and_index(plain,dest):
 """Use htslib-compatible tooling only; ordinary gzip is never accepted."""
 try:
  import pysam
  pysam.tabix_compress(str(plain),str(dest),force=True);pysam.tabix_index(str(dest),preset="vcf",force=True);return
 except ImportError:pass
 bgzip,tabix=shutil.which("bgzip"),shutil.which("tabix")
 if not bgzip or not tabix:raise RuntimeError("BGZF/index output requires pysam or both bgzip and tabix")
 with dest.open("wb") as h:subprocess.run([bgzip,"-c",str(plain)],stdout=h,check=True)
 subprocess.run([tabix,"-f","-p","vcf",str(dest)],check=True)
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("--config",type=Path,required=True);ap.add_argument("--overwrite",action="store_true");a=ap.parse_args();sec=read_simple_yaml(a.config).get("final_filter") or {}
 if sec.get("enabled",True) is False:print("[final_filter] disabled; skipping.");return 0
 out=resolve(sec.get("output_dir","results/qc/final_filter"));collected=resolve(sec.get("collected_dir","results/qc/collected_variant_calling_results"));collection_path=collected/"reports/variant_calling_collection_summary.tsv"
 if not collection_path.is_file():raise ValueError(f"missing collection summary: {collection_path}")
 sources=sec.get("sample_reports") or {};defaults={"intraspecies":("results/qc/intraspecies_contamination/reports/intraspecies_contamination_report.tsv",["contamination_status","qc_status"]),"human":("results/qc/human_contamination/reports/human_contamination_report.tsv",["human_contamination_status","qc_status","status"]),"interspecies":("results/qc/interspecies_contamination/reports/interspecies_contamination_report.tsv",["interspecies_status","qc_status","status"]),"sample_qc":("results/qc/sample_variant_filtering/reports/sample_qc.tsv",["qc_status"])}
 required=set(names(sec.get("required_sample_reports",["intraspecies","sample_qc"])));optional=set(names(sec.get("optional_sample_reports",["human","interspecies"])))
 unknown=(required|optional)-set(defaults)
 if unknown:raise ValueError(f"unknown sample reports: {sorted(unknown)}")
 indexed={};fields={};missing=[]
 for name,(fallback,candidates) in defaults.items():
  spec=sources.get(name,{}) if isinstance(sources,dict) else {};p=resolve(spec.get("path",fallback));fields[name]=spec.get("status_columns",candidates)
  if name in required and not p.is_file():missing.append(f"{name}={p}")
  indexed[name]=index(p) if p.is_file() else {}
 if missing:raise ValueError("missing required sample report(s): "+", ".join(missing))
 collection=index(collection_path);strict=sec.get("strict_missing_samples",True) is not False
 for report in required:
  absent=sorted(set(collection)-set(indexed[report]))
  if absent and strict:raise ValueError(f"required report {report} is missing samples: {', '.join(absent)}")
 if out.exists() and a.overwrite:shutil.rmtree(out)
 for d in (out/"reports",out/"logs",out/"final_vcf",out/"final_cov",out/"final_mtcn"):d.mkdir(parents=True,exist_ok=True)
 fail_cfg=sec.get("sample_fail_status") or {};fail_defaults={"intraspecies":["high_confidence_contaminated"],"human":["FAIL"],"interspecies":["FAIL"],"sample_qc":["FAIL"]};vcf_dirs=[resolve(x) for x in names(sec.get("vcf_sources",["results/qc/rrna_match/vcf_rrna","results/qc/trna_match/vcf_trna","results/qc/codon_match/vcf_codon","results/qc/coordinate_liftover/vcf_lifted_raw"]))]
 sample_rows=[];passing={};sample_qc_rows=indexed["sample_qc"]
 for sample,row in sorted(collection.items()):
  statuses={n:pick(indexed[n].get(sample,{}),fields[n]) for n in defaults};reasons=[];warnings=[]
  for n,v in statuses.items():
   if is_fail(v,fail_cfg.get(n,fail_defaults[n])):
    if n=="sample_qc":
     failed=pick(sample_qc_rows.get(sample,{}),["failed_criteria"],"failed")
     reasons.extend("sample_qc:"+x for x in failed.split(";") if x)
    else:reasons.append(n+":"+v)
   elif n=="intraspecies" and (v.startswith("insufficient_") or v=="candidate_contaminated"):warnings.append(n+":"+v)
   elif v=="NOT_AVAILABLE" and n in required:reasons.append(n+":missing")
  src=next((x for d in vcf_dirs if (x:=find_vcf(d,sample))),None)
  if not reasons and src is None:reasons.append("vcf:missing_downstream_source")
  status="FAIL" if reasons else "PASS"
  if status=="PASS":passing[sample]=src
  sample_rows.append(dict(sample=sample,species=pick(row,["species","Species"],""),intraspecies_status=statuses["intraspecies"],human_contamination_status=statuses["human"],interspecies_status=statuses["interspecies"],sample_level_qc_status=statuses["sample_qc"],final_sample_status=status,final_sample_fail_reasons=";".join(reasons),final_sample_warnings=";".join(warnings),vcf_source=str(src or "")))
 variant_flags=defaultdict(list)
 for source,spec in (sec.get("variant_reports") or {}).items():
  p=resolve(spec["path"] if isinstance(spec,dict) else spec)
  if not p.is_file():continue
  for r in read_tsv(p):
   k=(pick(r,["sample","Sample"],""),pick(r,["original_chrom","human_chrom","CHROM","chrom"],""),pick(r,["original_pos","human_pos","POS","pos"],""),pick(r,["original_ref","human_ref","REF","ref"],""),pick(r,["original_alt","human_alt","ALT","alt"],""));st=pick(r,(spec.get("status_columns",[]) if isinstance(spec,dict) else [])+["qc_status","status","match_status"],"PASS")
   if is_fail(st,spec.get("fail_status",["FAIL"]) if isinstance(spec,dict) else ["FAIL"]):variant_flags[k].append(source+":"+st)
 variant_rows=[];kept={}
 for sample,src in sorted(passing.items()):
  with tempfile.NamedTemporaryFile("w",suffix=".vcf",delete=False,dir=out) as target:
   plain=Path(target.name);n=0
   with open_vcf(src) as inp:
    for line in inp:
     if line.startswith("#"):target.write(line);continue
     f=line.rstrip("\n").split("\t");k=(sample,f[0],f[1],f[3],f[4]);why=variant_flags.get(k,[]);st="FAIL" if why else "PASS";variant_rows.append(dict(sample=sample,original_chrom=f[0],original_pos=f[1],original_ref=f[3],original_alt=f[4],human_chrom=f[0],human_pos=f[1],human_ref=f[3],human_alt=f[4],liftover_status="PASS",human_contamination_status="NOT_AVAILABLE",interspecies_status="NOT_AVAILABLE",sample_variant_qc_status="PASS" if not why else "FAIL",match_status="NOT_AVAILABLE",final_variant_status=st,final_variant_fail_reasons=";".join(why)))
     if st=="PASS":target.write(line);n+=1
  dest=out/"final_vcf"/f"{sample}.final.vcf.gz"
  try:bgzip_and_index(plain,dest)
  finally:plain.unlink(missing_ok=True)
  kept[sample]=n
  for kind in ("cov","mtcn"):
   for source in sorted((collected/f"collected_{kind}").glob(sample+".*")):shutil.copy2(source,out/f"final_{kind}"/source.name)
 with (out/"reports/final_sample_qc.tsv").open("w",newline="",encoding="utf-8") as h:w=csv.DictWriter(h,fieldnames=SAMPLE_COLUMNS,delimiter="\t");w.writeheader();w.writerows(sample_rows)
 with (out/"reports/final_variant_qc.tsv").open("w",newline="",encoding="utf-8") as h:w=csv.DictWriter(h,fieldnames=VARIANT_COLUMNS,delimiter="\t");w.writeheader();w.writerows(variant_rows)
 with (out/"reports/final_filter_summary.tsv").open("w",newline="",encoding="utf-8") as h:w=csv.writer(h,delimiter="\t");w.writerow(("metric","value"));w.writerows((("n_samples",len(collection)),("n_samples_pass",len(passing)),("n_samples_fail",len(collection)-len(passing)),("n_variants_pass",sum(kept.values())),("n_variants_fail",sum(r["final_variant_status"]=="FAIL" for r in variant_rows))))
 (out/"logs/final_filter.log").write_text(f"samples={len(collection)} pass={len(passing)} variants={len(variant_rows)}\n");print(f"[final_filter] output={out} samples_pass={len(passing)}/{len(collection)}");return 0
if __name__=="__main__":
 try:raise SystemExit(main())
 except (ValueError,KeyError,OSError,RuntimeError,subprocess.CalledProcessError) as e:print(f"ERROR: {e}",file=sys.stderr);raise SystemExit(2)
