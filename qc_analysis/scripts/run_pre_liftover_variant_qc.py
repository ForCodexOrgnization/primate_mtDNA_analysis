#!/usr/bin/env python3
"""Annotate native-coordinate VCFs with immutable source-call provenance."""
from __future__ import annotations

import argparse, csv, fcntl, gzip, math, shutil, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from qc_analysis.lib.simple_yaml import read_simple_yaml

TAGS = {
    "SOURCE_CHROM": ("String", "Original species-reference chromosome"),
    "SOURCE_POS": ("Integer", "Original species-reference position"),
    "SOURCE_REF": ("String", "Original species-reference REF allele"),
    "SOURCE_ALT": ("String", "Original species-reference ALT allele"),
    "SOURCE_AF": ("Float", "Original species-reference ALT allele frequency"),
    "SOURCE_DP": ("Integer", "Original species-reference depth"),
    "SOURCE_CALL_CLASS": ("String", "Biological class assigned in native coordinates"),
    "SOURCE_VARIANT_QC": ("String", "Native-coordinate variant QC status/reasons"),
}
FIELDS = "sample source_chrom source_pos source_ref source_alt source_af source_dp source_filter source_call_class source_variant_qc".split()

def op(path): return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()
def num(x):
    try:
        v=float(str(x).split(",")[0]); return v if math.isfinite(v) else None
    except (ValueError, TypeError): return None
def evidence(p):
    info={x.split("=",1)[0]:x.split("=",1)[1] for x in p[7].split(";") if "=" in x}
    fmt=dict(zip(p[8].split(":"),p[9].split(":"))) if len(p)>9 else {}
    af=num(fmt.get("AF",info.get("AF"))); dp=num(fmt.get("DP",info.get("DP")))
    if af is None:
        ad=[num(x) for x in fmt.get("AD","").split(",")]
        if len(ad)==2 and None not in ad and sum(ad)>0: af=ad[1]/sum(ad)
    return info,af,dp
def fmt(v): return f"{v:.10g}" if isinstance(v,float) else str(v)

def rebuild(out):
    reports=out/"reports"; lock=reports/".merge.lock"
    with lock.open("w") as h:
        fcntl.flock(h,fcntl.LOCK_EX)
        rows=[]
        for path in sorted(reports.glob("*.source_variant_qc.tsv")):
            with path.open() as f: rows.extend(csv.DictReader(f,delimiter="\t"))
        with (reports/"source_variant_qc.tsv").open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=FIELDS,delimiter="\t");w.writeheader();w.writerows(rows)
        counts={}
        for r in rows: counts[(r["sample"],r["source_call_class"],r["source_variant_qc"])]=counts.get((r["sample"],r["source_call_class"],r["source_variant_qc"]),0)+1
        with (reports/"source_variant_qc_summary.tsv").open("w",newline="") as f:
            w=csv.writer(f,delimiter="\t");w.writerow(("sample","source_call_class","source_variant_qc","n_variants"))
            for key,n in sorted(counts.items()):w.writerow((*key,n))

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--config",type=Path,required=True);ap.add_argument("--sample");a=ap.parse_args()
    sec=read_simple_yaml(a.config).get("pre_liftover_variant_qc",{})
    if sec.get("enabled",True) is False:return 0
    inp=ROOT/str(sec.get("input_vcf_dir","results/qc/collected_variant_calling_results/collected_vcf"));out=ROOT/str(sec.get("output_dir","results/qc/pre_liftover_variant_qc"))
    (out/"vcf_source_qc").mkdir(parents=True,exist_ok=True);(out/"reports").mkdir(parents=True,exist_ok=True)
    samples=[a.sample] if a.sample else sorted({p.name.split(".round2",1)[0] for p in inp.glob("*.round2.original_coords.clean.final.split.vcf*")})
    for sample in samples:
        choices=[inp/f"{sample}.round2.original_coords.clean.final.split.vcf.gz",inp/f"{sample}.round2.original_coords.clean.final.split.vcf"]
        src=next((p for p in choices if p.is_file()),None)
        if not src: raise FileNotFoundError(f"source VCF not found for {sample}")
        dest=out/"vcf_source_qc"/f"{sample}.source_qc.vcf.gz"; rows=[]
        with op(src) as ih,gzip.open(dest,"wt") as oh:
            for line in ih:
                if line.startswith("#CHROM"):
                    for tag,(typ,desc) in TAGS.items():oh.write(f'##INFO=<ID={tag},Number=1,Type={typ},Description="{desc}">\n')
                    oh.write(line);continue
                if line.startswith("#"):oh.write(line);continue
                p=line.rstrip().split("\t");info,af,dp=evidence(p)
                canonical=("," not in p[4] and len(p[3])==len(p[4])==1 and p[3].upper() in "ACGT" and p[4].upper() in "ACGT")
                klass="UNCLASSIFIED" if af is None else "LOW_AF" if af<float(sec.get("heteroplasmy_af_min",.10)) else "HET" if af<float(sec.get("homoplasmy_af_min",.95)) else "HOM"
                reasons=[]
                if sec.get("pass_only",True) and p[6]!="PASS":reasons.append("FILTER_NOT_PASS")
                if dp is None or dp<float(sec.get("dp_min",100)):reasons.append("DP_BELOW_MIN")
                if sec.get("snv_only",True) and not canonical:reasons.append("NOT_CANONICAL_BIALLELIC_SNV")
                qc="PASS" if not reasons else ";".join(reasons)
                vals={"SOURCE_CHROM":p[0],"SOURCE_POS":p[1],"SOURCE_REF":p[3],"SOURCE_ALT":p[4],"SOURCE_AF":"." if af is None else fmt(af),"SOURCE_DP":"." if dp is None else fmt(dp),"SOURCE_CALL_CLASS":klass,"SOURCE_VARIANT_QC":qc}
                p[7]=(p[7]+";" if p[7]!="." else "")+";".join(f"{k}={v}" for k,v in vals.items())
                oh.write("\t".join(p)+"\n")
                rows.append(dict(zip(FIELDS,[sample,p[0],p[1],vals["SOURCE_REF"],vals["SOURCE_ALT"],vals["SOURCE_AF"],vals["SOURCE_DP"],p[6],klass,qc])))
        with (out/"reports"/f"{sample}.source_variant_qc.tsv").open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=FIELDS,delimiter="\t");w.writeheader();w.writerows(rows)
    rebuild(out);return 0
if __name__=="__main__":raise SystemExit(main())
