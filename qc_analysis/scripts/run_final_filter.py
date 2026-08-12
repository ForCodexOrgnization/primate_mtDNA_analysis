#!/usr/bin/env python3
"""Centralize final sample/variant decisions and materialize final QC files."""
from __future__ import annotations
import argparse, csv, gzip, shutil, sys
from collections import defaultdict
from pathlib import Path
from typing import Any

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from qc_analysis.lib.simple_yaml import read_simple_yaml

SAMPLE_COLUMNS="sample species intraspecies_status human_contamination_status interspecies_status sample_level_qc_status final_sample_status final_sample_fail_reasons".split()
VARIANT_COLUMNS="sample original_chrom original_pos original_ref original_alt human_chrom human_pos human_ref human_alt liftover_status human_contamination_status interspecies_status sample_variant_qc_status match_status final_variant_status final_variant_fail_reasons".split()

def resolve(value: Any)->Path:
 p=Path(str(value)).expanduser();return p if p.is_absolute() else ROOT/p
def read_tsv(p:Path):
 with p.open(newline="",encoding="utf-8") as h:return list(csv.DictReader(h,delimiter="\t"))
def pick(row,names,default="NOT_AVAILABLE"):
 if isinstance(names,str):names=[x.strip() for x in names.split(",") if x.strip()]
 low={k.lower():v for k,v in row.items()}
 return next((low[x.lower()] for x in names if low.get(x.lower(),"")!=""),default)
def sample_name(p:Path)->str:return p.name.split(".vcf",1)[0].split(".",1)[0]
def open_vcf(p:Path,mode="rt"):
 return gzip.open(p,mode) if p.suffix==".gz" else p.open(mode,encoding="utf-8")
def report_index(path:Path)->dict[str,dict[str,str]]:
 if not path.is_file():return {}
 return {pick(r,["sample","Sample"]):r for r in read_tsv(path)}
def is_fail(value:str, configured:list[Any])->bool:
 if isinstance(configured,str):configured=[x.strip() for x in configured.split(",") if x.strip()]
 return value.strip().lower() in {str(x).strip().lower() for x in configured}

def main()->int:
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("--config",type=Path,required=True);ap.add_argument("--overwrite",action="store_true");a=ap.parse_args()
 cfg=read_simple_yaml(a.config);sec=cfg.get("final_filter") or {}
 if sec.get("enabled",True) is False:print("[final_filter] disabled; skipping.");return 0
 out=resolve(sec.get("output_dir","results/qc/final_filter")); reports=out/"reports"
 if out.exists() and a.overwrite: shutil.rmtree(out)
 for d in (reports,out/"logs",out/"final_vcf",out/"final_cov",out/"final_mtcn"):d.mkdir(parents=True,exist_ok=True)
 collected=resolve(sec.get("collected_dir","results/qc/collected_variant_calling_results")); vcf_dir=collected/"collected_vcf"
 collection=report_index(collected/"reports/variant_calling_collection_summary.tsv")
 sources=sec.get("sample_reports") or {}
 defaults={"intraspecies":("results/qc/intraspecies_contamination/reports/intraspecies_contamination_report.tsv",["contamination_status","qc_status"]),
  "human":("results/qc/human_contamination/reports/human_contamination_report.tsv",["human_contamination_status","qc_status","status"]),
  "interspecies":("results/qc/interspecies_contamination/reports/interspecies_contamination_report.tsv",["interspecies_status","qc_status","status"]),
  "sample_qc":("results/qc/sample_variant_filtering/reports/sample_qc.tsv",["sample_level_qc_status","qc_status","status"])}
 indexed={}; fields={}
 for name,(fallback,candidates) in defaults.items():
  spec=sources.get(name,{}) if isinstance(sources,dict) else {}; indexed[name]=report_index(resolve(spec.get("path",fallback)));fields[name]=spec.get("status_columns",candidates)
 fail_cfg=sec.get("sample_fail_status") or {}
 fail_defaults={"intraspecies":["high_confidence_contaminated","FAIL"],"human":["fail","FAIL"],"interspecies":["fail","FAIL"],"sample_qc":["fail","FAIL","MISSING_INPUT","FAIL_PROCESSING"]}
 files=sorted(set(vcf_dir.glob("*.vcf"))|set(vcf_dir.glob("*.vcf.gz")))
 samples=sorted(set(collection)|{sample_name(p) for p in files}); vcf_by_sample={sample_name(p):p for p in files}
 sample_rows=[]; passing=set()
 for sample in samples:
  species=pick(collection.get(sample,{}),["species","Species"],"")
  statuses={name:pick(indexed[name].get(sample,{}),fields[name]) for name in defaults}
  if statuses["sample_qc"]=="NOT_AVAILABLE": statuses["sample_qc"]=pick(collection.get(sample,{}),["status"],"NOT_AVAILABLE")
  reasons=[name+":"+value for name,value in statuses.items() if is_fail(value,fail_cfg.get(name,fail_defaults[name]))]
  status="FAIL" if reasons else "PASS"
  if status=="PASS":passing.add(sample)
  sample_rows.append(dict(sample=sample,species=species,intraspecies_status=statuses["intraspecies"],human_contamination_status=statuses["human"],interspecies_status=statuses["interspecies"],sample_level_qc_status=statuses["sample_qc"],final_sample_status=status,final_sample_fail_reasons=";".join(reasons)))
 with (reports/"final_sample_qc.tsv").open("w",newline="",encoding="utf-8") as h:w=csv.DictWriter(h,fieldnames=SAMPLE_COLUMNS,delimiter="\t");w.writeheader();w.writerows(sample_rows)
 # Optional variant reports are joined by sample + original allele. Status FAIL
 # is recognized from a stable family of common status column names.
 variant_specs=sec.get("variant_reports") or {}; variant_flags=defaultdict(list)
 for source,spec in variant_specs.items():
  p=resolve(spec["path"] if isinstance(spec,dict) else spec)
  if not p.is_file():continue
  for row in read_tsv(p):
   key=(pick(row,["sample","Sample"],""),pick(row,["original_chrom","CHROM","chrom"],""),pick(row,["original_pos","POS","pos"],""),pick(row,["original_ref","REF","ref"],""),pick(row,["original_alt","ALT","alt"],""))
   status=pick(row,(spec.get("status_columns",[]) if isinstance(spec,dict) else [])+["qc_status","status","match_status"],"PASS")
   if is_fail(status,(spec.get("fail_status",["FAIL","fail"]) if isinstance(spec,dict) else ["FAIL","fail"])):variant_flags[key].append(source+":"+status)
 variant_rows=[]; kept_counts={}
 for sample in sorted(passing):
  src=vcf_by_sample.get(sample)
  if not src:continue
  dest=out/"final_vcf"/src.name; kept=0
  with open_vcf(src) as inp, open_vcf(dest,"wt") as target:
   for line in inp:
    if line.startswith("#"):target.write(line);continue
    f=line.rstrip("\n").split("\t");key=(sample,f[0],f[1],f[3],f[4]);reasons=variant_flags.get(key,[]);status="FAIL" if reasons else "PASS"
    variant_rows.append(dict(sample=sample,original_chrom=f[0],original_pos=f[1],original_ref=f[3],original_alt=f[4],human_chrom="",human_pos="",human_ref="",human_alt="",liftover_status="NOT_AVAILABLE",human_contamination_status="NOT_AVAILABLE",interspecies_status="NOT_AVAILABLE",sample_variant_qc_status="PASS" if not reasons else "FAIL",match_status="NOT_AVAILABLE",final_variant_status=status,final_variant_fail_reasons=";".join(reasons)))
    if status=="PASS":target.write(line);kept+=1
  kept_counts[sample]=kept
  for kind in ("cov","mtcn"):
   candidates=sorted((collected/f"collected_{kind}").glob(sample+".*"))
   for source in candidates:shutil.copy2(source,out/f"final_{kind}"/source.name)
 with (reports/"final_variant_qc.tsv").open("w",newline="",encoding="utf-8") as h:w=csv.DictWriter(h,fieldnames=VARIANT_COLUMNS,delimiter="\t");w.writeheader();w.writerows(variant_rows)
 with (reports/"final_filter_summary.tsv").open("w",newline="",encoding="utf-8") as h:
  w=csv.writer(h,delimiter="\t");w.writerow(("metric","value"));w.writerows((("n_samples",len(samples)),("n_samples_pass",len(passing)),("n_samples_fail",len(samples)-len(passing)),("n_variants_pass",sum(kept_counts.values())),("n_variants_fail",sum(r["final_variant_status"]=="FAIL" for r in variant_rows))))
 (out/"logs/final_filter.log").write_text(f"samples={len(samples)} pass={len(passing)} variants={len(variant_rows)}\n",encoding="utf-8")
 print(f"[final_filter] output={out} samples_pass={len(passing)}/{len(samples)}");return 0
if __name__=="__main__":
 try:raise SystemExit(main())
 except (ValueError,KeyError,OSError) as e:print(f"ERROR: {e}",file=sys.stderr);raise SystemExit(2)
