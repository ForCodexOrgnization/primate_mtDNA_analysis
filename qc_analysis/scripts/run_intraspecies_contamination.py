#!/usr/bin/env python3
"""Produce the cohort-level, sample-only intra-species contamination report.

This is the production implementation of the three-evidence algorithm formerly
implemented in ``validation/contamination_reference.R``. It intentionally uses no
third-party Python packages and never modifies an input VCF.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import math
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from qc_analysis.lib.simple_yaml import read_simple_yaml

REPORT_COLUMNS = """sample species n_species_samples n_usable_variants n_lowA
best_source_sample best_overlap best_frac_lowA_in_highB
n_anchor_pool_excluding_A n_anchor_tested_in_A n_depressed_anchor
mt_high_hets_contamination mt_high_hets_mode anchor_evidence_level
n_mirror_pairs n_low_variants_with_mirror mirror_low_fraction
mirror_support_candidate mirror_support_highconf contamination_status
contamination_flag_candidate contamination_flag_highconf qc_status qc_reason""".split()

DEFAULTS = dict(dp_min=100, low_vaf_min=.01, low_vaf_max=.20,
 high_vaf_min=.99, mt_lower=.80, mt_depressed_upper=.998, mt_anchor_upper=1.0,
 min_n_lowA=5, min_overlap=3, min_frac_lowA_in_highB_candidate=.50,
 min_frac_lowA_in_highB_highconf=.6213636363636358,
 contam_threshold_candidate=.036420574377757434,
 contam_threshold_highconf=.07103935483870959, mirror_low_vaf_min=.01,
 mirror_low_vaf_max=.20, mirror_high_vaf_min=.80,
 mirror_high_vaf_max=.998, mirror_tolerance=0,
 min_mirror_pairs_for_raw_flag=3, min_low_variants_with_mirror_for_flag=3)

def path(value: Any) -> Path:
    p = Path(str(value)).expanduser()
    return p if p.is_absolute() else ROOT / p

def truth(value: Any, default: bool = False) -> bool:
    if value is None: return default
    if isinstance(value, bool): return value
    raise ValueError(f"expected boolean, got {value!r}")

def load_rows(p: Path) -> list[dict[str, str]]:
    with p.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def number(row: dict[str, str], key: str) -> float:
    return float(row[key])

def variant_key(row: dict[str, str]) -> tuple[str, str, str, str]:
    return row["CHROM"], row["POS"], row["REF"], row["ALT"]

def analyse(rows: list[dict[str, str]], parameters: dict[str, Any], sample_pairs: set[tuple[str,str]] | None = None) -> list[dict[str, Any]]:
    """Apply lowA/highB, leave-one-out anchor, and mirror evidence rules."""
    samples = sorted({(r["Species"], r["Sample"]) for r in rows} | (sample_pairs or set()))
    usable = [r for r in rows if number(r, "DP") >= float(parameters["dp_min"])
              and (not parameters["use_snv_only"] or r["Type"] == "SNV")
              and (not parameters["pass_only"] or r["FILTER"] == "PASS")]
    by_sample: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in usable: by_sample[row["Species"], row["Sample"]].append(row)
    species_counts = Counter(species for species, _ in samples)
    result = []
    for species, sample in samples:
        own = by_sample[species, sample]
        others = [s for sp, s in samples if sp == species and s != sample]
        low = [r for r in own if parameters["low_vaf_min"] <= number(r,"VAF") <= parameters["low_vaf_max"]]
        low_keys = {variant_key(r) for r in low}
        best_source, best_overlap = "", 0
        other_high: dict[str, set[tuple[str,str,str,str]]] = {}
        for source in others:
            keys = {variant_key(r) for r in by_sample[species,source] if number(r,"VAF") >= parameters["high_vaf_min"]}
            other_high[source] = keys
            overlap = len(low_keys & keys)
            if overlap > best_overlap: best_source, best_overlap = source, overlap
        anchor_pool = set().union(*other_high.values()) if other_high else set()
        own_map = defaultdict(list)
        for r in own: own_map[variant_key(r)].append(number(r,"VAF"))
        tested = [v for key in anchor_pool for v in own_map.get(key, [])]
        depressed = [v for v in tested if parameters["mt_lower"] <= v <= parameters["mt_depressed_upper"]]
        fallback = [v for v in tested if parameters["mt_lower"] <= v <= parameters["mt_anchor_upper"]]
        if len(depressed) >= 3: mode, estimate = "depressed_anchors", 1-sum(depressed)/len(depressed)
        elif fallback: mode, estimate = "fallback_anchors", 1-sum(fallback)/len(fallback)
        else: mode, estimate = "no_anchor_observed", None
        # Preserve the reference algorithm's within-sample, same-allele mirror
        # definition; duplicate low/high calls are paired only when AFs sum to 1.
        lows = [r for r in own if parameters["mirror_low_vaf_min"] <= number(r,"VAF") <= parameters["mirror_low_vaf_max"]]
        highs = [r for r in own if parameters["mirror_high_vaf_min"] <= number(r,"VAF") <= parameters["mirror_high_vaf_max"]]
        pairs = [(a,b) for a in lows for b in highs if variant_key(a)==variant_key(b)
                 and abs(number(a,"VAF")+number(b,"VAF")-1) <= parameters["mirror_tolerance"]]
        mirrored = len({variant_key(a) for a,_ in pairs})
        frac = mirrored/len({variant_key(x) for x in lows}) if lows else 0.0
        overlap_frac = best_overlap/len(low_keys) if low_keys else None
        common = len(low_keys)>=parameters["min_n_lowA"] and best_overlap>=parameters["min_overlap"]
        candidate = bool(common and overlap_frac is not None and overlap_frac>=parameters["min_frac_lowA_in_highB_candidate"] and estimate is not None and estimate>=parameters["contam_threshold_candidate"])
        highconf = bool(common and overlap_frac is not None and overlap_frac>=parameters["min_frac_lowA_in_highB_highconf"] and estimate is not None and estimate>=parameters["contam_threshold_highconf"])
        mirror_candidate = len(pairs)>=parameters["min_mirror_pairs_for_raw_flag"] and mirrored>=parameters["min_low_variants_with_mirror_for_flag"]
        mirror_highconf = mirror_candidate and highconf
        if not others: status="insufficient_singleton_species"
        elif not own: status="insufficient_variant_data"
        elif highconf: status="high_confidence_contaminated"
        elif candidate: status="candidate_contaminated"
        elif estimate is None: status="insufficient_anchor_data"
        elif best_overlap>=parameters["min_overlap"]: status="lowA_highB_overlap_only"
        elif estimate>=parameters["contam_threshold_candidate"]: status="mt_high_hets_only"
        else: status="no_strong_evidence"
        qc = "FAIL" if status=="high_confidence_contaminated" else "PASS" if status=="no_strong_evidence" else "WARN"
        result.append(dict(sample=sample,species=species,n_species_samples=species_counts[species],n_usable_variants=len(own),n_lowA=len(low_keys),
          best_source_sample=best_source,best_overlap=best_overlap,best_frac_lowA_in_highB=overlap_frac,
          n_anchor_pool_excluding_A=len(anchor_pool),n_anchor_tested_in_A=len(tested),n_depressed_anchor=len(depressed),
          mt_high_hets_contamination=estimate,mt_high_hets_mode=mode,anchor_evidence_level="strong" if len(depressed)>=3 else "limited" if tested else "none",
          n_mirror_pairs=len(pairs),n_low_variants_with_mirror=mirrored,mirror_low_fraction=frac,
          mirror_support_candidate=mirror_candidate,mirror_support_highconf=mirror_highconf,
          contamination_status=status,contamination_flag_candidate=candidate,contamination_flag_highconf=highconf,qc_status=qc,qc_reason=status))
    return result

def main() -> int:
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument("--config",type=Path); ap.add_argument("--variant-table",type=Path); ap.add_argument("--outdir",type=Path); ap.add_argument("--overwrite",action="store_true"); args=ap.parse_args()
    cfg={}; section={}
    if args.config:
        cfg=read_simple_yaml(args.config)
        if "intraspecies_contamination" not in cfg: raise ValueError("missing 'intraspecies_contamination' section in configuration")
        section=cfg.get("intraspecies_contamination") or {}
        if not isinstance(section,dict): raise ValueError("'intraspecies_contamination' must be a YAML mapping")
        if not truth(section.get("enabled"),False):
            print("[intraspecies] enabled=false")
            for key in ("build_variant_table","vcf_dir","metadata","variant_table","outdir"):
                value=section.get(key);print(f"[intraspecies] {key}={str(value).lower() if isinstance(value,bool) else value if value not in (None,'') else '<not set>'}")
            print("[intraspecies] disabled; skipping."); return 0
    out=args.outdir or path(section.get("outdir","results/qc/intraspecies_contamination")); out=path(out)
    report=out/"reports/intraspecies_contamination_report.tsv"
    if report.exists() and not (args.overwrite or truth(section.get("overwrite"),False)): raise ValueError(f"output exists: {report}; use --overwrite")
    for d in (out/"logs",out/"reports"): d.mkdir(parents=True,exist_ok=True)
    table=args.variant_table or (path(section["variant_table"]) if section.get("variant_table") else None)
    sample_pairs: set[tuple[str,str]] = set()
    if table is None and truth(section.get("build_variant_table"),True):
        vcf=path(section.get("vcf_dir","results/qc/collected_variant_calling_results/collected_vcf"))
        metadata=path(section.get("sample_summary","results/qc/collected_variant_calling_results/reports/variant_calling_collection_summary.tsv"))
        table=out/".work/all_PASS_variants_core_table.tsv"; table.parent.mkdir(parents=True,exist_ok=True)
        cmd=[sys.executable,str(ROOT/"qc_analysis/scripts/build_intraspecies_variant_table.py"),"--vcf-dir",str(vcf),"--metadata",str(metadata),"--output",str(table),"--min-dp",str(section.get("dp_min",100)),"--pass-only","--overwrite","--log-file",str(out/"logs/variant_table_build.log")]
        if truth(section.get("use_snv_only"),True): cmd.append("--snv-only")
        subprocess.run(cmd,check=True)
        for row in load_rows(metadata):
            sample=row.get("sample") or row.get("Sample") or ""; species=row.get("species") or row.get("Species") or ""
            if sample and species: sample_pairs.add((species,sample))
    if table is None: raise ValueError("build_variant_table=false requires variant_table")
    parameters={**DEFAULTS,**{k:v for k,v in section.items() if k in DEFAULTS}}
    parameters.update(use_snv_only=truth(section.get("use_snv_only"),True),pass_only=truth(section.get("pass_only"),True))
    rows=load_rows(path(table)); findings=analyse(rows,parameters,sample_pairs)
    with report.open("w",newline="",encoding="utf-8") as h:
        w=csv.DictWriter(h,fieldnames=REPORT_COLUMNS,delimiter="\t",extrasaction="ignore");w.writeheader();w.writerows(findings)
    with (out/"run_parameters.tsv").open("w",newline="",encoding="utf-8") as h:
        w=csv.writer(h,delimiter="\t");w.writerow(("parameter","value"));w.writerows(sorted(parameters.items()));w.writerow(("variant_table",table));w.writerow(("timestamp",dt.datetime.now(dt.timezone.utc).isoformat()))
    (out/"logs/intraspecies_contamination.log").write_text(f"samples={len(findings)}\nreport={report}\n",encoding="utf-8")
    print(f"[intraspecies] report={report} samples={len(findings)}")
    return 0

if __name__=="__main__":
    try: raise SystemExit(main())
    except (ValueError,OSError,subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}",file=sys.stderr);raise SystemExit(2)
