#!/usr/bin/env python3
"""Build the original-coordinate core table for intra-species QC.

This deliberately uses a small streaming VCF parser rather than INFO fields or a
coordinate conversion library: FORMAT/AF and FORMAT/DP are the sample-level
measurements used by the downstream analysis.
"""
from __future__ import annotations
import argparse, csv, gzip, logging, os, sys
from collections import defaultdict
from pathlib import Path

COLUMNS = "Sample Species CHROM POS REF ALT Type FILTER DP VAF AD_ref AD_alt".split()
SUMMARY = "sample species vcf_file n_total_records n_pass_records n_dp_ge_threshold n_snv_records n_output_variants n_missing_af n_missing_dp n_malformed_ad status notes".split()

def open_text(path): return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()
def clean(v): return "" if v is None else str(v).strip()
def metadata(path, sample_col, species_col):
    with Path(path).open(newline="") as h:
        first=h.readline(); h.seek(0); dialect=csv.excel_tab if "\t" in first else csv.excel
        rows=list(csv.reader(h,dialect=dialect))
    if not rows: raise ValueError("metadata is empty")
    header=[x.strip() for x in rows[0]]
    low=[x.lower() for x in header]
    if sample_col in header: si=header.index(sample_col)
    elif sample_col.lower() in low: si=low.index(sample_col.lower())
    else: si=0; rows=[header]+rows[1:] # headerless existing sample_ref_file format
    if species_col in header: pi=header.index(species_col)
    elif species_col.lower() in low: pi=low.index(species_col.lower())
    else: pi=1
    start=1 if (sample_col in header or sample_col.lower() in low or species_col in header or species_col.lower() in low) else 0
    result={}
    for row in rows[start:]:
        if len(row)<=max(si,pi): continue
        s,p=clean(row[si]),clean(row[pi])
        if not s or not p: raise ValueError("metadata contains an empty sample or species")
        if s in result and result[s]!=p: raise ValueError(f"sample maps to more than one species: {s}")
        result[s]=p
    return result
def expected_name(path):
    name=path.name
    for suffix in (".vcf.gz", ".vcf"): 
        if name.endswith(suffix): return name[:-len(suffix)].split(".")[0]
    return name
def main():
 p=argparse.ArgumentParser(); p.add_argument("--vcf-dir",required=True);p.add_argument("--metadata",required=True);p.add_argument("--output",required=True);p.add_argument("--metadata-sample-column",default="sample");p.add_argument("--metadata-species-column",default="species");p.add_argument("--min-dp",type=float,default=0);p.add_argument("--pass-only",action="store_true");p.add_argument("--snv-only",action="store_true");p.add_argument("--overwrite",action="store_true");p.add_argument("--log-file");a=p.parse_args()
 out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True)
 if out.exists() and not a.overwrite: p.error(f"output exists: {out}; use --overwrite")
 logpath=Path(a.log_file) if a.log_file else out.parent/"variant_table_build_warnings.log"; logging.basicConfig(filename=logpath,level=logging.WARNING,format="%(levelname)s: %(message)s")
 try: meta=metadata(a.metadata,a.metadata_sample_column,a.metadata_species_column)
 except ValueError as e: p.error(str(e))
 files=sorted({x.resolve() for x in Path(a.vcf_dir).rglob("*.vcf")}|{x.resolve() for x in Path(a.vcf_dir).rglob("*.vcf.gz")})
 if not files: p.error("no .vcf or .vcf.gz files found")
 seen={}; rows=[]; summaries=[]
 for f in files:
  st=defaultdict(int); st.update(status="ok",notes=""); sample=None
  try:
   with open_text(f) as h:
    for line in h:
     if line.startswith("#CHROM"):
      fields=line.rstrip("\n").split("\t")
      if len(fields)!=10: raise ValueError("VCF must have exactly one sample column")
      sample=fields[9]
      if sample != expected_name(f): raise ValueError(f"VCF sample {sample!r} does not match filename-derived sample {expected_name(f)!r}")
      if sample not in meta: raise ValueError(f"VCF sample missing from metadata: {sample}")
      if sample in seen and seen[sample] != f: raise ValueError(f"two VCF files resolve to sample {sample}: {seen[sample]} and {f}")
      seen[sample]=f; continue
     if not line.startswith("#"):
      if sample is None: raise ValueError("record before #CHROM header")
      st["n_total_records"]+=1; q=line.rstrip("\n").split("\t")
      if len(q)<10: raise ValueError("malformed VCF record")
      chrom,pos,ref,alts,flt,fmt,call=q[0],q[1],q[3],q[4],q[6],q[8],q[9]
      if flt=="PASS": st["n_pass_records"]+=1
      if a.pass_only and flt!="PASS": continue
      keys=fmt.split(":"); vals=call.split(":")
      if len(vals)!=len(keys): logging.warning("%s:%s malformed sample genotype fields",f,pos); continue
      d=dict(zip(keys,vals)); dp=d.get("DP")
      if dp in (None,".",""): st["n_missing_dp"]+=1; logging.warning("%s:%s missing FORMAT/DP",f,pos); continue
      try: dp=float(dp)
      except ValueError: st["n_missing_dp"]+=1; logging.warning("%s:%s malformed FORMAT/DP",f,pos); continue
      if dp>=a.min_dp: st["n_dp_ge_threshold"]+=1
      else: continue
      af=d.get("AF")
      if af in (None,".",""): st["n_missing_af"]+=1; logging.warning("%s:%s missing FORMAT/AF",f,pos); continue
      av=af.split(","); altv=alts.split(","); ad=d.get("AD"); adv=ad.split(",") if ad not in (None,".","") else None
      for i,alt in enumerate(altv):
       try: vaf=float(av[i])
       except (ValueError,IndexError): st["n_missing_af"]+=1; logging.warning("%s:%s missing/malformed AF for ALT %s",f,pos,alt); continue
       adref=adalt=""
       if adv:
        if len(adv)<i+2: st["n_malformed_ad"]+=1; logging.warning("%s:%s malformed AD",f,pos)
        else: adref,adalt=adv[0],adv[i+1]
       typ="SNV" if len(ref)==len(alt)==1 else "INDEL"
       if typ=="SNV": st["n_snv_records"]+=1
       if a.snv_only and typ!="SNV": continue
       rows.append([sample,meta[sample],chrom,pos,ref,alt,typ,flt,dp,vaf,adref,adalt]); st["n_output_variants"]+=1
   if sample is None: raise ValueError("no #CHROM header")
  except Exception as e:
   st["status"]="error";st["notes"]=str(e); logging.error("%s",e); summaries.append([sample or "",meta.get(sample or "",""),str(f)]+[st[x] for x in SUMMARY[3:-2]]+[st["status"],st["notes"]]); continue
  summaries.append([sample,meta[sample],str(f)]+[st[x] for x in SUMMARY[3:-2]]+[st["status"],st["notes"]])
 if any(x[-2]=="error" for x in summaries):
  with (out.parent/"variant_table_build_summary.tsv").open("w",newline="") as h: csv.writer(h,delimiter="\t").writerows([SUMMARY]+summaries)
  raise SystemExit("VCF table build failed; see summary and warnings log")
 rows.sort(key=lambda x:(x[1],x[0],x[2],int(x[3]),x[4],x[5]))
 with out.open("w",newline="") as h: csv.writer(h,delimiter="\t").writerows([COLUMNS]+rows)
 with (out.parent/"variant_table_build_summary.tsv").open("w",newline="") as h: csv.writer(h,delimiter="\t").writerows([SUMMARY]+summaries)
if __name__=="__main__": main()
