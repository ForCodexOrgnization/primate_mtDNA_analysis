#!/usr/bin/env python3
"""Strict terminal sample/variant filtering of the most downstream VCF."""
from __future__ import annotations
import argparse,csv,gzip,math,re,shutil,subprocess,sys,tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from qc_analysis.lib.simple_yaml import read_simple_yaml
SAMPLE_COLUMNS="sample species intraspecies_status human_contamination_status interspecies_status sample_level_qc_status final_sample_status final_sample_fail_reasons final_sample_warnings vcf_source".split()
VARIANT_COLUMNS=("sample species human_chrom human_pos human_ref human_alt source_chrom source_pos source_ref source_alt "
 "AF DP vcf_filter variant_class call_class snv_type mt_median_coverage Percent_100 nuclear_median_coverage mtcn_median MAD "
 "sample_level_qc_status sample_failed_criteria intraspecies_status human_contamination_status interspecies_status "
 "region_type orthology_match_status orthology_fail_reason codon_match_status trna_match_status rrna_match_status "
 "final_variant_status final_variant_fail_reasons original_chrom original_pos original_ref original_alt liftover_status sample_variant_qc_status match_status").split()
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
def parse_info(value):return {item.split("=",1)[0]:item.split("=",1)[1] if "=" in item else True for item in value.split(";") if item and item!="."}
def number(value):
 try:
  result=float(value);return result if math.isfinite(result) else None
 except (TypeError,ValueError):return None
def variant_evidence(fields):
 info=parse_info(fields[7]);format_names=fields[8].split(":") if len(fields)>8 else [];format_values=fields[9].split(":") if len(fields)>9 else [];sample=dict(zip(format_names,format_values))
 af=number(sample.get("AF","").split(",")[0]);dp=number(sample.get("DP"))
 if dp is None:dp=number(info.get("DP"))
 if af is None and "," not in fields[4]:
  ad=[number(value) for value in sample.get("AD","").split(",")]
  if len(ad)==2 and None not in ad and sum(ad)>0:af=ad[1]/sum(ad)
 return info,af,dp
def variant_classes(ref,alt):
 ref,alt=ref.upper(),alt.upper();simple="," not in alt
 if simple and len(ref)==len(alt)==1:variant_class="SNV"
 elif simple and len(ref)!=len(alt) and all(re.fullmatch(r"[ACGTN]+",allele) for allele in (ref,alt)):variant_class="INDEL"
 else:variant_class="OTHER"
 if variant_class=="SNV" and {ref,alt} in ({"A","G"},{"C","T"}):snv_type="SNV_transition"
 elif variant_class=="SNV" and ref in "ACGT" and alt in "ACGT":snv_type="SNV_transversion"
 elif variant_class=="INDEL":snv_type="indel"
 else:snv_type="other"
 return variant_class,snv_type
def call_class(af):
 if af is None:return "UNKNOWN"
 if af>=.95:return "homoplasmic"
 if af>=.10:return "heteroplasmic"
 return "low_af"
def info_value(info,*aliases):return next((str(info[name]) for name in aliases if info.get(name) not in (None,"",".")),"NOT_AVAILABLE")
LEGACY_SUFFIXES=(".lifted.codon.trna.rrna.vcf",".lifted.codon.trna.vcf",".lifted.trna.vcf",".lifted.codon.vcf",".lifted.raw.vcf",".vcf")
def source_specs(value):
 """Normalize named exact-pattern sources and the legacy directory list."""
 if isinstance(value,dict):
  return [(name,resolve(spec.get("dir",spec.get("directory"))),spec.get("pattern")) for name,spec in value.items() if isinstance(spec,dict)]
 return [(str(d),resolve(d),None) for d in names(value)]
def find_vcf(directory,sample,pattern=None):
 """Resolve only exact sample filenames; never use a sample-prefix glob."""
 base=[directory/(pattern.format(sample=sample))] if pattern else [directory/f"{sample}{suffix}" for suffix in LEGACY_SUFFIXES]
 candidates=[]
 for path in base:
  candidates.extend(p for p in (path,Path(str(path)+".gz")) if p.is_file())
 candidates=sorted(set(candidates))
 if len(candidates)>1:raise ValueError(f"ambiguous VCF source for sample {sample} in {directory}: {', '.join(map(str,candidates))}")
 return candidates[0] if candidates else None
def report_variant_key(row,spec,source):
 """Return the canonical sample + post-liftover human allele identity."""
 sample=pick(row,["sample","Sample"],"")
 human=[pick(row,[f"human_{x}"],"") for x in ("chrom","pos","ref","alt")]
 if all(human):return (sample,*human)
 system=str(spec.get("coordinate_system","")).strip().lower()
 if system not in {"human","post-liftover","post_liftover"}:
  raise ValueError(f"variant report {source} uses generic coordinates but coordinate_system is unknown or incompatible: {system or 'not declared'}")
 generic=[pick(row,[x.upper(),x],"") for x in ("chrom","pos","ref","alt")]
 if not all(generic):raise ValueError(f"variant report {source} lacks a complete human-coordinate variant key")
 return (sample,*generic)
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
def sort_plain_vcf(input_vcf:Path,output_vcf:Path)->None:
 """Write a coordinate-sorted VCF while preserving header lines verbatim."""
 contig_order={};headers=[];records=[]
 with input_vcf.open("r",encoding="utf-8",newline="") as source:
  for line_number,line in enumerate(source,1):
   if line.startswith("#"):
    headers.append(line)
    match=re.match(r"^##contig=<ID=([^,>]+)",line)
    if match:
     contig=match.group(1).strip().strip('"')
     if contig not in contig_order:contig_order[contig]=len(contig_order)
    continue
   fields=line.rstrip("\r\n").split("\t")
   if len(fields)<5:raise ValueError(f"invalid VCF data line {line_number} in {input_vcf}")
   try:pos=int(fields[1])
   except ValueError as exc:raise ValueError(f"non-integer VCF POS on line {line_number} in {input_vcf}: {fields[1]!r}") from exc
   chrom,ref,alt=fields[0],fields[3],fields[4]
   contig_key=(0,contig_order[chrom]) if chrom in contig_order else (1,chrom)
   records.append(((contig_key,pos,ref,alt,line),line))
 records.sort(key=lambda item:item[0])
 with output_vcf.open("w",encoding="utf-8",newline="") as target:
  target.writelines(headers)
  target.writelines(line for _key,line in records)
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
 fail_cfg=sec.get("sample_fail_status") or {};fail_defaults={"intraspecies":["high_confidence_contaminated"],"human":["FAIL"],"interspecies":["FAIL"],"sample_qc":["FAIL"]};vcf_sources=source_specs(sec.get("vcf_sources",["results/qc/rrna_match/vcf_rrna","results/qc/trna_match/vcf_trna","results/qc/codon_match/vcf_codon","results/qc/coordinate_liftover/vcf_lifted_raw"]))
 sample_rows=[];sample_context={};passing={};sample_qc_rows=indexed["sample_qc"]
 for sample,row in sorted(collection.items()):
  statuses={n:pick(indexed[n].get(sample,{}),fields[n]) for n in defaults};reasons=[];warnings=[]
  for n,v in statuses.items():
   if is_fail(v,fail_cfg.get(n,fail_defaults[n])):
    if n=="sample_qc":
     failed=pick(sample_qc_rows.get(sample,{}),["failed_criteria"],"failed")
     reasons.extend("sample_qc:"+x for x in failed.split(";") if x)
    else:reasons.append(n+":"+v)
   elif n=="intraspecies" and (v.startswith("insufficient_") or v=="candidate_contaminated"):warnings.append(n+":"+v)
   elif n=="human" and v.upper() in {"CANDIDATE","INSUFFICIENT_DATA"}:warnings.append(n+":"+v)
   elif v=="NOT_AVAILABLE" and n in required:reasons.append(n+":missing")
  src=next((x for _,d,p in vcf_sources if (x:=find_vcf(d,sample,p))),None)
  if not reasons and src is None:reasons.append("vcf:missing_downstream_source")
  status="FAIL" if reasons else "PASS"
  if status=="PASS":passing[sample]=src
  sample_row=dict(sample=sample,species=pick(row,["species","Species"],""),intraspecies_status=statuses["intraspecies"],human_contamination_status=statuses["human"],interspecies_status=statuses["interspecies"],sample_level_qc_status=statuses["sample_qc"],final_sample_status=status,final_sample_fail_reasons=";".join(reasons),final_sample_warnings=";".join(warnings),vcf_source=str(src or ""));sample_rows.append(sample_row);sample_context[sample]=sample_row
 variant_flags=defaultdict(list);variant_annotations=defaultdict(dict)
 for source,spec in (sec.get("variant_reports") or {}).items():
  p=resolve(spec["path"] if isinstance(spec,dict) else spec)
  if not p.is_file():continue
  for r in read_tsv(p):
   if not isinstance(spec,dict):raise ValueError(f"variant report {source} must declare path and coordinate_system")
   k=report_variant_key(r,spec,source);st=pick(r,names(spec.get("status_columns",[]))+["qc_status","status","match_status"],"PASS")
   for field in ("region_type","orthology_match_status","orthology_fail_reason"):
    value=pick(r,[field],"")
    if value:variant_annotations[k][field]=value
   if is_fail(st,spec.get("fail_status",["FAIL"]) if isinstance(spec,dict) else ["FAIL"]):variant_flags[k].append(source+":"+st)
 variant_rows=[];kept={}
 for sample,src in sorted(passing.items()):
  with tempfile.NamedTemporaryFile("w",suffix=".vcf",delete=False,dir=out) as target:
   plain=Path(target.name);n=0
   with open_vcf(src) as inp:
    for line in inp:
     if line.startswith("#"):target.write(line);continue
     f=line.rstrip("\n").split("\t");k=(sample,f[0],f[1],f[3],f[4]);why=list(variant_flags.get(k,[]));vcf_filter=f[6]
     if vcf_filter!="PASS":why.append("vcf_filter:"+vcf_filter)
     info,af,dp=variant_evidence(f);variant_class,snv_type=variant_classes(f[3],f[4]);ref,alt=f[3].upper(),f[4].upper();is_canonical_snv=len(ref)==len(alt)==1 and ref in "ACGT" and alt in "ACGT" and "," not in alt
     if not is_canonical_snv:why.append("variant_class:"+variant_class)
     st="FAIL" if why else "PASS";context=sample_context[sample];qc=sample_qc_rows.get(sample,{});annotation=variant_annotations.get(k,{});source_chrom=info_value(info,"SRC_CHROM","MTLIFT_ORIG_CHROM");source_pos=info_value(info,"SRC_POS","MTLIFT_ORIG_POS");source_ref=info_value(info,"SRC_REF","MTLIFT_ORIG_REF");source_alt=info_value(info,"SRC_ALT","MTLIFT_ORIG_ALT");orthology_status=annotation.get("orthology_match_status","NOT_AVAILABLE")
     variant_rows.append(dict(sample=sample,species=context["species"],human_chrom=f[0],human_pos=f[1],human_ref=f[3],human_alt=f[4],source_chrom=source_chrom,source_pos=source_pos,source_ref=source_ref,source_alt=source_alt,AF=af if af is not None else "NA",DP=dp if dp is not None else "NA",vcf_filter=vcf_filter,variant_class=variant_class,call_class=call_class(af),snv_type=snv_type,mt_median_coverage=pick(qc,["mt_median_coverage"]),Percent_100=pick(qc,["Percent_100"]),nuclear_median_coverage=pick(qc,["nuclear_median_coverage"]),mtcn_median=pick(qc,["mtcn_median"]),MAD=pick(qc,["MAD"]),sample_level_qc_status=context["sample_level_qc_status"],sample_failed_criteria=pick(qc,["failed_criteria"],""),intraspecies_status=context["intraspecies_status"],human_contamination_status=context["human_contamination_status"],interspecies_status=context["interspecies_status"],region_type=annotation.get("region_type","NOT_AVAILABLE"),orthology_match_status=orthology_status,orthology_fail_reason=annotation.get("orthology_fail_reason","NOT_AVAILABLE"),codon_match_status=info_value(info,"MTCODON_STATUS"),trna_match_status=info_value(info,"MTTRNA_STATUS"),rrna_match_status=info_value(info,"MTRRNA_STATUS"),final_variant_status=st,final_variant_fail_reasons=";".join(why),original_chrom=source_chrom,original_pos=source_pos,original_ref=source_ref,original_alt=source_alt,liftover_status="PASS",sample_variant_qc_status="PASS" if not why else "FAIL",match_status=orthology_status))
     if st=="PASS":target.write(line);n+=1
  dest=out/"final_vcf"/f"{sample}.final.vcf.gz";sorted_plain=None
  try:
   with tempfile.NamedTemporaryFile("w",suffix=".sorted.vcf",delete=False,dir=out) as sorted_target:sorted_plain=Path(sorted_target.name)
   sort_plain_vcf(plain,sorted_plain)
   bgzip_and_index(sorted_plain,dest)
  finally:
   plain.unlink(missing_ok=True)
   if sorted_plain is not None:sorted_plain.unlink(missing_ok=True)
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
