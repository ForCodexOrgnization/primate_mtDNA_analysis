#!/usr/bin/env python3
"""Validated original-coordinate intra-species contamination analysis."""
from __future__ import annotations
import argparse,csv,datetime as dt,math,subprocess,sys
from collections import Counter,defaultdict
from pathlib import Path
from typing import Any
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from qc_analysis.lib.simple_yaml import read_simple_yaml
REPORT_COLUMNS="""sample species n_species_samples n_usable_variants n_lowA best_source_sample best_overlap best_frac_lowA_in_highB n_anchor_pool_excluding_A n_anchor_tested_in_A n_depressed_anchor mt_high_hets_contamination mt_high_hets_mode anchor_evidence_level anchor_source_count n_mirror_pairs n_low_variants_with_mirror mirror_low_fraction normalized_mirror_support mirror_p95_threshold mirror_p99_threshold mirror_calibration_status mirror_support_candidate mirror_support_highconf contamination_status contamination_flag_candidate contamination_flag_highconf qc_status qc_reason""".split()
DEFAULTS=dict(dp_min=100,low_vaf_min=.01,low_vaf_max=.20,high_vaf_min=.99,mt_lower=.80,mt_depressed_upper=.998,mt_anchor_upper=1.,min_n_lowA=5,min_overlap=3,min_frac_lowA_in_highB_candidate=.50,min_frac_lowA_in_highB_highconf=.6213636363636358,contam_threshold_candidate=.036420574377757434,contam_threshold_highconf=.07103935483870959,mirror_low_vaf_min=.01,mirror_low_vaf_max=.20,mirror_high_vaf_min=.80,mirror_high_vaf_max=.998,mirror_tolerance=0.,min_negative_control_values=3,target_negative_control_tier="tier2_location_and_batch_different")
def path(v):
 p=Path(str(v)).expanduser();return p if p.is_absolute() else ROOT/p
def truth(v,d=False):
 if v is None:return d
 if isinstance(v,bool):return v
 raise ValueError(f"expected boolean, got {v!r}")
def load_rows(p):
 with p.open(newline="",encoding="utf-8") as h:return list(csv.DictReader(h,delimiter="\t"))
def num(r,k):return float(r[k])
def key(r):return r["CHROM"],r["POS"],r["REF"],r["ALT"]
def quantile7(values,p):
 x=sorted(values);h=(len(x)-1)*p;i=math.floor(h);return x[i]+(x[min(i+1,len(x)-1)]-x[i])*(h-i)
def mirror_stats(rows,p):
 lows=[r for r in rows if p["mirror_low_vaf_min"]<=num(r,"VAF")<=p["mirror_low_vaf_max"]]
 highs=[r for r in rows if p["mirror_high_vaf_min"]<=num(r,"VAF")<=p["mirror_high_vaf_max"]]
 pairs=[(i,j) for i,a in enumerate(lows) for j,b in enumerate(highs) if abs(num(a,"VAF")+num(b,"VAF")-1)<=p["mirror_tolerance"]+1e-12]
 mirrored=len({i for i,_ in pairs}); frac=mirrored/len(lows) if lows else 0.
 normalized=len(pairs)/(len(lows)*len(highs)) if lows and highs else 0.
 return len(pairs),mirrored,frac,normalized

def calibration(rows,p,nc_path):
 if not nc_path:return None,None,"not_calibrated_no_file",0
 q=path(nc_path)
 if not q.is_file():return None,None,"not_calibrated_missing_file",0
 table=load_rows(q);tier=[r for r in table if r.get("negative_control_tier")==str(p["target_negative_control_tier"])]
 if not tier:return None,None,"not_calibrated_no_tier2_pairs",0
 by=defaultdict(list)
 for r in rows:by[r["Sample"]].append(r)
 vals=[]
 for pair in tier:
  a,b=by.get(pair.get("Sample_A",""),[]),by.get(pair.get("Sample_B",""),[])
  if a and b:vals.append(mirror_stats(a+b,p)[3])
 if len(vals)<int(p["min_negative_control_values"]):return None,None,"not_calibrated_insufficient_values",len(vals)
 return quantile7(vals,.95),quantile7(vals,.99),"calibrated",len(vals)

def analyse(rows,p,sample_pairs=None,negative_control_pairs=None):
 samples=sorted({(r["Species"],r["Sample"]) for r in rows}|(sample_pairs or set()));usable=[r for r in rows if num(r,"DP")>=float(p["dp_min"]) and (not p["use_snv_only"] or r["Type"]=="SNV") and (not p["pass_only"] or r["FILTER"]=="PASS")]
 by=defaultdict(list)
 for r in usable:by[r["Species"],r["Sample"]].append(r)
 p95,p99,cal_status,_=calibration(usable,p,negative_control_pairs);counts=Counter(x for x,_ in samples);out=[]
 for species,sample in samples:
  own=by[species,sample];others=[s for sp,s in samples if sp==species and s!=sample];low=[r for r in own if p["low_vaf_min"]<=num(r,"VAF")<=p["low_vaf_max"]];lowkeys={key(r) for r in low};otherhigh={s:{key(r) for r in by[species,s] if num(r,"VAF")>=p["high_vaf_min"]} for s in others}
  ranked=[(len(lowkeys&ks),s) for s,ks in otherhigh.items()];best_overlap,best_source=max(ranked,key=lambda x:(x[0],x[1])) if ranked else (0,"");anchor_pool=set().union(*otherhigh.values()) if otherhigh else set();ownmap=defaultdict(list)
  for r in own:ownmap[key(r)].append(num(r,"VAF"))
  tested=[v for k in anchor_pool for v in ownmap.get(k,[])];dep=[v for v in tested if p["mt_lower"]<=v<=p["mt_depressed_upper"]];fall=[v for v in tested if p["mt_lower"]<=v<=p["mt_anchor_upper"]]
  if len(dep)>=3:mode,est="depressed_anchors",1-sum(dep)/len(dep)
  elif fall:mode,est="fallback_anchors",1-sum(fall)/len(fall)
  else:mode,est="no_anchor_observed",None
  source_count=sum(bool(lowkeys&ks) for ks in otherhigh.values())
  if not others:level="singleton_species"
  elif not anchor_pool:level="no_loo_anchor_pool"
  elif not tested:level="no_anchor_observed_in_sample_A"
  elif len(dep)>=3:level=("strong_multi_source_anchor_support" if source_count>1 else "strong_single_source_anchor_support")
  elif len(tested)>=2:level=("moderate_multi_source_anchor_support" if source_count>1 else "moderate_single_source_anchor_support")
  else:level="weak_anchor_support"
  npair,nmir,mfrac,norm=mirror_stats(own,p);frac=best_overlap/len(lowkeys) if lowkeys else None;overlap_candidate=len(lowkeys)>=p["min_n_lowA"] and best_overlap>=p["min_overlap"] and frac is not None and frac>=p["min_frac_lowA_in_highB_candidate"];overlap_high=overlap_candidate and frac>=p["min_frac_lowA_in_highB_highconf"]
  mtcand=est is not None and est>=p["contam_threshold_candidate"];mthigh=est is not None and est>=p["contam_threshold_highconf"];candidate=overlap_candidate and mtcand;highconf=overlap_high and mthigh
  mcand=p95 is not None and norm>=p95;mhigh=p99 is not None and norm>=p99
  if not others:status="insufficient_singleton_species"
  elif not own:status="insufficient_variant_data"
  elif highconf:status="high_confidence_contaminated"
  elif candidate:status="candidate_contaminated"
  elif overlap_candidate and not mtcand:status="lowA_highB_overlap_only"
  elif mtcand and not overlap_candidate:status="mt_high_hets_only"
  elif est is None:status="insufficient_anchor_data"
  else:status="no_strong_evidence"
  qc="FAIL" if status=="high_confidence_contaminated" else "PASS" if status=="no_strong_evidence" else "WARN"
  out.append(dict(sample=sample,species=species,n_species_samples=counts[species],n_usable_variants=len(own),n_lowA=len(lowkeys),best_source_sample=best_source,best_overlap=best_overlap,best_frac_lowA_in_highB=frac,n_anchor_pool_excluding_A=len(anchor_pool),n_anchor_tested_in_A=len(tested),n_depressed_anchor=len(dep),mt_high_hets_contamination=est,mt_high_hets_mode=mode,anchor_evidence_level=level,anchor_source_count=source_count,n_mirror_pairs=npair,n_low_variants_with_mirror=nmir,mirror_low_fraction=mfrac,normalized_mirror_support=norm,mirror_p95_threshold=p95,mirror_p99_threshold=p99,mirror_calibration_status=cal_status,mirror_support_candidate=mcand,mirror_support_highconf=mhigh,contamination_status=status,contamination_flag_candidate=candidate,contamination_flag_highconf=highconf,qc_status=qc,qc_reason=status))
 return out

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("--config",type=Path);ap.add_argument("--variant-table",type=Path);ap.add_argument("--outdir",type=Path);ap.add_argument("--negative-control-pairs",type=Path);ap.add_argument("--overwrite",action="store_true");a=ap.parse_args();sec={}
 if a.config:
  cfg=read_simple_yaml(a.config)
  if "intraspecies_contamination" not in cfg:raise ValueError("missing 'intraspecies_contamination' section in configuration")
  sec=cfg.get("intraspecies_contamination") or {}
  if not isinstance(sec,dict):raise ValueError("'intraspecies_contamination' must be a YAML mapping")
  if not truth(sec.get("enabled"),False):
   print("[intraspecies] enabled=false")
   for k in ("build_variant_table","vcf_dir","metadata","variant_table","outdir"):
    v=sec.get(k);print(f"[intraspecies] {k}={str(v).lower() if isinstance(v,bool) else v if v not in (None,'') else '<not set>'}")
   print("[intraspecies] disabled; skipping.");return 0
 out=path(a.outdir or sec.get("outdir","results/qc/intraspecies_contamination"));report=out/"reports/intraspecies_contamination_report.tsv"
 if report.exists():print(f"[intraspecies] replacing existing report: {report}",file=sys.stderr)
 for d in (out/"logs",out/"reports"):d.mkdir(parents=True,exist_ok=True)
 table=a.variant_table or (path(sec["variant_table"]) if sec.get("variant_table") else None);sample_pairs=set()
 if table is None and truth(sec.get("build_variant_table"),True):
  vcf=path(sec.get("vcf_dir","results/qc/collected_variant_calling_results/collected_vcf"));meta=path(sec.get("sample_summary","results/qc/collected_variant_calling_results/reports/variant_calling_collection_summary.tsv"));table=out/".work/all_PASS_variants_core_table.tsv";table.parent.mkdir(parents=True,exist_ok=True);cmd=[sys.executable,str(ROOT/"qc_analysis/scripts/build_intraspecies_variant_table.py"),"--vcf-dir",str(vcf),"--metadata",str(meta),"--output",str(table),"--min-dp",str(sec.get("dp_min",100)),"--pass-only","--overwrite","--log-file",str(out/"logs/variant_table_build.log")];cmd += ["--snv-only"] if truth(sec.get("use_snv_only"),True) else [];subprocess.run(cmd,check=True)
  for r in load_rows(meta):
   s=r.get("sample") or r.get("Sample") or "";sp=r.get("species") or r.get("Species") or ""
   if s and sp:sample_pairs.add((sp,s))
 if table is None:raise ValueError("build_variant_table=false requires variant_table")
 p={**DEFAULTS,**{k:v for k,v in sec.items() if k in DEFAULTS},"use_snv_only":truth(sec.get("use_snv_only"),True),"pass_only":truth(sec.get("pass_only"),True)};nc=a.negative_control_pairs or sec.get("negative_control_pairs");findings=analyse(load_rows(path(table)),p,sample_pairs,nc)
 with report.open("w",newline="",encoding="utf-8") as h:w=csv.DictWriter(h,fieldnames=REPORT_COLUMNS,delimiter="\t",extrasaction="ignore");w.writeheader();w.writerows(findings)
 with (out/"run_parameters.tsv").open("w",newline="",encoding="utf-8") as h:w=csv.writer(h,delimiter="\t");w.writerow(("parameter","value"));w.writerows(sorted(p.items()));w.writerow(("variant_table",table));w.writerow(("timestamp",dt.datetime.now(dt.timezone.utc).isoformat()))
 (out/"logs/intraspecies_contamination.log").write_text(f"samples={len(findings)}\nreport={report}\n",encoding="utf-8");print(f"[intraspecies] report={report} samples={len(findings)}");return 0
if __name__=="__main__":
 try:raise SystemExit(main())
 except (ValueError,OSError,subprocess.CalledProcessError) as e:print(f"ERROR: {e}",file=sys.stderr);raise SystemExit(2)
