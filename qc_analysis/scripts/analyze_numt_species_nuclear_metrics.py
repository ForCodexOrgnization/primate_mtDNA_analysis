#!/usr/bin/env python3
"""Report species/reference-level high-confidence nuclear NUMT locus counts and bp.

Uses precisely the samples included in sample_numt_region_base_coverage.tsv
(coverage tier highconf_any_overlap), and their original per-sample highconf
NUMT BED files. Merges BED (0-based, half-open) intervals across samples
for each cohort/species/reference, yielding unique nuclear NUMT loci/bases.
The comparison to published T2T values is descriptive only: the published
reference catalog uses identity >=90% while ours has additional read filters.
"""
import argparse
import csv
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
DEFAULT_BED=Path("/vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis/data/numt_besthit")

def read_tsv(path):
    with path.open(newline="",encoding="utf-8") as f:
        yield from csv.DictReader(f,delimiter="\t")

def merge(spans):
    spans.sort()
    out=[]
    for a,b in spans:
        if not out or a>out[-1][1]:
            out.append([a,b])
        else:
            out[-1][1]=max(out[-1][1],b)
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--sample-coverage",type=Path,default=ROOT/"results/numt_region_base_coverage/sample_numt_region_base_coverage.tsv")
    ap.add_argument("--species-coverage",type=Path,default=ROOT/"results/numt_region_base_coverage/species_union_numt_region_base_coverage.tsv")
    ap.add_argument("--highconf-dir",type=Path,default=DEFAULT_BED)
    ap.add_argument("--output",type=Path,default=ROOT/"results/numt_region_base_coverage/species_union_numt_nuclear_metrics.tsv")
    args=ap.parse_args()
    included=defaultdict(set)
    for r in read_tsv(args.sample_coverage):
        if r["tier"]!="highconf_any_overlap":continue
        included[(r["cohort"],r["species"],r["reference_key"])].add(r["sample"])
    coverage={}
    for r in read_tsv(args.species_coverage):
        if r["tier"]!="highconf_any_overlap":continue
        key=(r["cohort"],r["species"],r["reference_key"])
        if key in coverage:raise ValueError(f"Duplicate species coverage: {key}")
        coverage[key]=r
    fields=["cohort","species","reference_key","n_samples","n_numt_loci",
            "nuclear_numt_union_bp","nuclear_numt_sum_bp","chrM_coverage_pct","mt_length","n_empty_bed"]
    rows=[]
    for key,samples in sorted(included.items()):
        if key not in coverage:raise ValueError(f"No matching species coverage for {key}")
        bychrom=defaultdict(list)
        n_empty=0
        for sample in samples:
            bed=args.highconf_dir/f"{sample}.highconf_numt.bed"
            if not bed.is_file():raise FileNotFoundError(bed)
            found=0
            with bed.open(encoding="utf-8") as stream:
                for line in stream:
                    if not line.strip() or line.startswith("#"):continue
                    f=line.split()
                    if len(f)<3:continue
                    chrom,a,b=f[0],int(f[1]),int(f[2])
                    if b<=a:raise ValueError(f"Invalid BED interval in {bed}: {line}")
                    bychrom[chrom].append((a,b))
                    found+=1
            if not found:n_empty+=1
        merged={chrom:merge(v) for chrom,v in bychrom.items()}
        n=sum(len(v) for v in merged.values())
        unionbp=sum(b-a for v in merged.values() for a,b in v)
        sumbp=sum(b-a for v in bychrom.values() for a,b in v)
        c=coverage[key]
        rows.append(dict(cohort=key[0],species=key[1],reference_key=key[2],
                         n_samples=len(samples),n_numt_loci=n,
                         nuclear_numt_union_bp=unionbp,nuclear_numt_sum_bp=sumbp,
                         chrM_coverage_pct=c["coverage_pct"],mt_length=c["mt_length"],
                         n_empty_bed=n_empty))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,delimiter="\t",fieldnames=fields)
        w.writeheader();w.writerows(rows)
    print(f"Wrote {len(rows)} cohort/species/reference rows to {args.output}")
    for r in rows:
        if r["species"] in ("Pan_paniscus","Pan_troglodytes","Gorilla_gorilla","Pongo_pygmaeus","Pongo_abelii"):
            print(r["cohort"],r["species"],"n_loci",r["n_numt_loci"],"nuclear_bp",r["nuclear_numt_union_bp"],"chrM_cov",r["chrM_coverage_pct"])

if __name__=="__main__":main()
