#!/usr/bin/env python3
"""Read-only audit of duplicate round2 native-coordinate mtDNA variants."""
from __future__ import annotations
import argparse
import csv
import gzip
import math
from collections import Counter, defaultdict
from pathlib import Path

SUFFIXES=(".round2.original_coords.clean.final.split.vcf.gz",".round2.original_coords.clean.final.split.vcf")
GROUP_FIELDS="sample chrom pos ref alt n_records n_excess_records n_vcf_pass n_pass_dp50 filter_pattern distinct_filters min_dp max_dp different_cons_pos n_distinct_cons_pos cons_positions all_records_low_dp any_pass_dp50".split()
RECORD_FIELDS="sample chrom pos ref alt record_number source_file line_number dp af filter cons_pos cons_ref cons_alt liftback_status qual id line".split()

def number(x):
    if x in (None,"",".","NA","nan"):
        return None
    try:
        n=float(str(x).split(",")[0])
        return n if math.isfinite(n) else None
    except (ValueError,TypeError):
        return None

def observe(f):
    info=dict(part.split("=",1) for part in f[7].split(";") if "=" in part)
    fmt=dict(zip(f[8].split(":"),f[9].split(":"))) if len(f)>9 else {}
    dp=number(fmt.get("DP",info.get("DP")))
    af=number(fmt.get("AF",info.get("AF")))
    if af is None and "AD" in fmt:
        ad=[number(v) for v in fmt["AD"].split(",")]
        if len(ad)==2 and None not in ad and sum(ad)>0:
            af=ad[1]/sum(ad)
    return dict(dp=dp,af=af,filter=f[6],cons_pos=info.get("CONS_POS",""),
                cons_ref=info.get("CONS_REF",""),cons_alt=info.get("CONS_ALT",""),
                liftback_status=info.get("LIFTBACK_STATUS",""),
                qual=number(f[5]),id=f[2],line="\t".join(f[:5]))

def discover(directory):
    files={}
    for suffix in SUFFIXES:
        for path in sorted(directory.glob("*"+suffix)):
            files.setdefault(path.name[:-len(suffix)],path)
    return files

def write_tsv(path,columns,rows):
    with path.open("w",newline="") as h:
        w=csv.DictWriter(h,delimiter="\t",fieldnames=columns,extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf-dir",type=Path,default=Path("results/qc_coverage_only/collected_variant_calling_results/collected_vcf"))
    p.add_argument("--outdir",type=Path,default=Path("results/qc_coverage_only_dp50/duplicate_variant_audit"))
    p.add_argument("--dp-min",type=float,default=50)
    p.add_argument("--sample")
    a=p.parse_args()
    if not a.vcf_dir.is_dir():
        p.error("Missing source VCF directory: "+str(a.vcf_dir))
    files=discover(a.vcf_dir)
    if a.sample:
        files={a.sample:files[a.sample]} if a.sample in files else {}
    if not files:
        p.error("No source VCF files to audit")
    a.outdir.mkdir(parents=True,exist_ok=True)
    samples=[]
    group_rows=[]
    record_rows=[]
    overall=Counter()
    for i,(sample,path) in enumerate(sorted(files.items()),1):
        by_key=defaultdict(list)
        opener=gzip.open if path.suffix==".gz" else open
        total=0
        with opener(path,"rt") as h:
            for line_no,line in enumerate(h,1):
                if not line or line.startswith("#"):
                    continue
                fields=line.rstrip("\r\n").split("\t")
                if len(fields)<8:
                    raise ValueError("Malformed VCF at %s:%d"%(path,line_no))
                key=tuple(fields[j] for j in (0,1,3,4))
                rec=observe(fields)
                rec["line_number"]=line_no
                by_key[key].append(rec)
                total+=1
        stats=Counter()
        for (chrom,pos,ref,alt),records in by_key.items():
            if len(records)<2:
                continue
            n=len(records)
            filters={r["filter"] for r in records}
            n_pass=sum(r["filter"]=="PASS" for r in records)
            n_pass_dp=sum(r["filter"]=="PASS" and r["dp"] is not None and r["dp"]>=a.dp_min for r in records)
            pattern="ALL_PASS" if n_pass==n else "MIXED_PASS_FAIL" if n_pass else "ALL_FAIL"
            cons={r["cons_pos"] for r in records if r["cons_pos"] not in ("",".")}
            known=all(r["cons_pos"] not in ("",".") for r in records)
            different="YES" if len(cons)>1 else "NO" if known else "UNKNOWN"
            depths=[r["dp"] for r in records if r["dp"] is not None]
            has_dp=any(d>=a.dp_min for d in depths)
            stats.update({"duplicate_groups":1,"duplicate_records":n,"excess_records":n-1,pattern:1,
                          "different_cons_pos_"+different:1,"any_dp_ge_threshold":int(has_dp),
                          "any_vcf_pass_dp_ge_threshold":int(n_pass_dp>0)})
            group_rows.append(dict(sample=sample,chrom=chrom,pos=pos,ref=ref,alt=alt,n_records=n,
                                   n_excess_records=n-1,n_vcf_pass=n_pass,n_pass_dp50=n_pass_dp,
                                   filter_pattern=pattern,distinct_filters=",".join(sorted(filters)),
                                   min_dp=min(depths,default=None),max_dp=max(depths,default=None),
                                   different_cons_pos=different,n_distinct_cons_pos=len(cons),
                                   cons_positions=",".join(sorted(cons,key=lambda v:(int(v) if v.isdigit() else 10**9,v))),
                                   all_records_low_dp="YES" if not has_dp else "NO",
                                   any_pass_dp50="YES" if n_pass_dp else "NO"))
            for j,rec in enumerate(records,1):
                record_rows.append(dict(sample=sample,chrom=chrom,pos=pos,ref=ref,alt=alt,
                                        record_number=j,source_file=str(path),**rec))
        samples.append(dict(sample=sample,source_vcf=str(path),n_total_records=total,
                            n_unique_keys=len(by_key),n_duplicate_groups=stats["duplicate_groups"],
                            n_duplicate_records=stats["duplicate_records"],n_excess_records=stats["excess_records"],
                            duplicate_rate_records=stats["duplicate_records"]/total if total else 0,
                            groups_all_pass=stats["ALL_PASS"],groups_mixed_pass_fail=stats["MIXED_PASS_FAIL"],
                            groups_all_fail=stats["ALL_FAIL"],
                            groups_different_cons_pos=stats["different_cons_pos_YES"],
                            groups_same_cons_pos=stats["different_cons_pos_NO"],
                            groups_cons_pos_unknown=stats["different_cons_pos_UNKNOWN"],
                            groups_any_dp_ge_threshold=stats["any_dp_ge_threshold"],
                            groups_any_pass_dp50=stats["any_vcf_pass_dp_ge_threshold"]))
        overall.update(stats)
        overall["vcf_records"]+=total
        overall["unique_keys"]+=len(by_key)
        overall["samples_with_duplicates"]+=int(stats["duplicate_groups"]>0)
        if i%100==0:
            print("[audit] %d/%d samples"%(i,len(files)),flush=True)
    write_tsv(a.outdir/"sample_duplicate_summary.tsv",list(samples[0]),samples)
    write_tsv(a.outdir/"duplicate_variant_groups.tsv",GROUP_FIELDS,group_rows)
    write_tsv(a.outdir/"duplicate_variant_records.tsv",RECORD_FIELDS,record_rows)
    metrics=dict(n_samples_scanned=len(files),n_samples_with_duplicates=overall["samples_with_duplicates"],
                 n_total_vcf_records=overall["vcf_records"],n_unique_variant_keys=overall["unique_keys"],
                 n_duplicate_groups=overall["duplicate_groups"],n_duplicate_records=overall["duplicate_records"],
                 n_excess_records=overall["excess_records"],n_groups_all_vcf_pass=overall["ALL_PASS"],
                 n_groups_mixed_vcf_pass_fail=overall["MIXED_PASS_FAIL"],
                 n_groups_all_vcf_fail=overall["ALL_FAIL"],
                 n_groups_different_cons_pos=overall["different_cons_pos_YES"],
                 n_groups_same_cons_pos=overall["different_cons_pos_NO"],
                 n_groups_cons_pos_unknown=overall["different_cons_pos_UNKNOWN"],
                 n_groups_with_any_dp_ge_threshold=overall["any_dp_ge_threshold"],
                 n_groups_with_any_vcf_pass_dp_ge_threshold=overall["any_vcf_pass_dp_ge_threshold"],
                 dp_threshold=a.dp_min)
    with (a.outdir/"cohort_duplicate_summary.tsv").open("w",newline="") as h:
        writer=csv.writer(h,delimiter="\t")
        writer.writerow(("metric","value"))
        writer.writerows(metrics.items())
    print("\n".join("%s: %s"%item for item in metrics.items()))
    print("Reports:",a.outdir.resolve())

if __name__=="__main__":
    main()
