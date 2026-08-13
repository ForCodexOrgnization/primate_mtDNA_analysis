#!/usr/bin/env python3
"""Screen post-liftover VCFs for low-frequency human mtDNA marker alleles.

This module deliberately uses the VCF POS/REF/ALT and sample FORMAT values:
SRC_* annotations describe the pre-liftover allele and are audit metadata only.
HaploGrep is supplementary and receives a synthetic marker-only profile.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterable

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from qc_analysis.lib.simple_yaml import read_simple_yaml

MARKER_RE = re.compile(r"^([0-9]+)([ACGT])(!?)$")
REPORT_COLUMNS = """sample species human_contamination_status human_contamination_evidence
n_usable_variants n_low_variants n_variants_missing_af n_human_marker_hits
frac_low_variants_human_marker baseline_marker_screen_pass
n_human_marker_hits_control_region n_human_marker_hits_non_control_region
frac_human_marker_hits_control_region non_control_marker_pass median_human_marker_af
mean_human_marker_af min_human_marker_af max_human_marker_af human_marker_af_mad
human_marker_af_iqr fraction_human_markers_near_median_af vaf_coherence_pass
estimated_human_fraction_median_af estimated_human_mtDNA_fraction
n_human_markers_with_ad human_marker_total_depth human_marker_total_alt_depth
human_marker_list n_back_mutation_marker_hits haplogrep_status haplogrep_best_haplogroup
haplogrep_quality haplogrep_n_input_markers haplogrep_sparse_input
haplogrep_support_status haplogrep_missing_markers haplogrep_private_markers
haplogrep_n_missing_markers haplogrep_n_private_markers""".split()
AUDIT_COLUMNS = """sample species human_chrom human_pos human_ref human_alt dp af
af_source ref_depth alt_depth liftover_allele_status human_marker is_back_mutation
in_control_region used_in_candidate_screen used_in_fail_screen used_in_haplogrep""".split()


def resolve(value: Any) -> Path:
    p = Path(str(value)).expanduser()
    return p if p.is_absolute() else ROOT / p


def boolean(value: bool) -> str:
    return "true" if value else "false"


def load_samples(path: Path) -> list[tuple[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows:
        return []
    header = [x.lower() for x in rows[0]]
    start = 1 if "sample" in header else 0
    si = header.index("sample") if start else 0
    pi = header.index("species") if start and "species" in header else 1
    result = [(r[si].strip(), r[pi].strip() if len(r) > pi else "") for r in rows[start:] if len(r) > si and r[si].strip()]
    if len({x[0] for x in result}) != len(result):
        raise ValueError(f"duplicate samples in {path}")
    return result


def load_markers(path: Path) -> tuple[dict[tuple[int, str], dict[str, Any]], dict[str, int]]:
    """Parse notation rather than trusting normalized columns; first POS/ALT wins."""
    markers: dict[tuple[int, str], dict[str, Any]] = {}
    qc = {"raw": 0, "simple": 0, "back": 0, "duplicates": 0, "excluded": 0}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "marker" not in {x.lower() for x in reader.fieldnames}:
            raise ValueError(f"marker table requires a marker column: {path}")
        for row in reader:
            qc["raw"] += 1
            low = {str(k).lower(): str(v or "").strip() for k, v in row.items()}
            notation = low.get("marker", "")
            match = MARKER_RE.fullmatch(notation)
            if not match:
                qc["excluded"] += 1
                continue
            pos, alt, bang = int(match.group(1)), match.group(2), bool(match.group(3))
            # A bang in notation is authoritative; normalized true is also retained.
            back = bang or low.get("is_back_mutation", "").lower() in {"true", "1", "yes"}
            qc["simple"] += 1
            qc["back"] += int(back)
            key = (pos, alt)
            if key in markers:
                qc["duplicates"] += 1
                # Preserve back-mutation distinction conservatively across duplicates.
                markers[key]["is_back_mutation"] |= back
                continue
            markers[key] = {"marker": notation, "pos": pos, "alt": alt, "is_back_mutation": back}
    if not markers:
        raise ValueError(f"marker table contains no usable simple SNVs: {path}")
    return markers, qc


def find_vcf(directory: Path, pattern: str, sample: str) -> Path | None:
    base = directory / pattern.format(sample=sample)
    candidates = [base]
    if base.suffix == ".gz":
        candidates.append(Path(str(base)[:-3]))
    else:
        candidates.append(Path(str(base) + ".gz"))
    found = [p for p in candidates if p.is_file()]
    if len(found) > 1:
        raise ValueError(f"ambiguous lifted VCF for {sample}: {', '.join(map(str, found))}")
    return found[0] if found else None


def parse_number(value: str) -> float | None:
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def parse_vcf(path: Path) -> Iterable[dict[str, Any]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            fmt = fields[8].split(":") if len(fields) > 8 else []
            values = fields[9].split(":") if len(fields) > 9 else []
            sample = dict(zip(fmt, values))
            info = {part.split("=", 1)[0]: part.split("=", 1)[1] if "=" in part else True for part in fields[7].split(";") if part}
            dp = parse_number(sample.get("DP", ""))
            if dp is None:
                dp = parse_number(str(info.get("DP", "")))
            ref_depth = alt_depth = None
            ad = sample.get("AD", "").split(",")
            if len(ad) == 2:
                ref_depth, alt_depth = parse_number(ad[0]), parse_number(ad[1])
                if ref_depth is None or alt_depth is None or ref_depth < 0 or alt_depth < 0:
                    ref_depth = alt_depth = None
            af = parse_number(sample.get("AF", "").split(",")[0])
            source = "FORMAT_AF"
            if af is None and ref_depth is not None and alt_depth is not None and ref_depth + alt_depth > 0:
                af = alt_depth / (ref_depth + alt_depth)
                source = "CALCULATED_FROM_AD"
            elif af is None:
                source = "MISSING"
            try:
                pos = int(fields[1])
            except ValueError:
                pos = None
            yield dict(chrom=fields[0], pos=pos, ref=fields[3].upper(), alt=fields[4].upper(),
                       filter=fields[6], dp=dp, af=af, af_source=source,
                       ref_depth=ref_depth, alt_depth=alt_depth,
                       liftover_allele_status=str(info.get("LIFTOVER_ALLELE_STATUS", "")))


def quantile(values: list[float], q: float) -> float:
    if len(values) == 1:
        return values[0]
    ordered = sorted(values); index = (len(ordered) - 1) * q
    lo, hi = math.floor(index), math.ceil(index)
    return ordered[lo] + (ordered[hi] - ordered[lo]) * (index - lo)


def haplogrep_tool(cfg: dict[str, Any]) -> tuple[list[str] | None, str]:
    executable = str(cfg.get("executable") or "").strip()
    jar = str(cfg.get("jar") or "").strip()
    if executable:  # documented deterministic precedence
        p = resolve(executable)
        return ([str(p)], "executable") if p.is_file() and os.access(p, os.X_OK) else (None, "executable")
    if jar:
        p = resolve(jar); java = str(cfg.get("java") or "java")
        java_path = shutil.which(java) if not Path(java).is_absolute() else (java if Path(java).is_file() else None)
        xmx = str(cfg.get("java_xmx") or "").strip()
        return ([str(java_path)] + ([f"-Xmx{xmx}"] if xmx else []) + ["-jar", str(p)], "jar") if p.is_file() and java_path else (None, "jar")
    return None, "unconfigured"


def parse_haplogrep_output(path: Path) -> dict[str, Any]:
    rows = []
    for delimiter in ("\t", ","):
        with path.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter=delimiter))
        if rows and len(rows[0]) > 1: break
    if not rows: return {}
    row = rows[0]; norm = {re.sub(r"[^a-z0-9]", "", k.lower()): v for k, v in row.items()}
    def get(*keys: str) -> str:
        return next((norm[k] for k in keys if norm.get(k) not in (None, "")), "")
    missing = get("missingpolys", "missingmarkers", "missing")
    private = get("privatepolys", "privatemarkers", "private")
    return {"haplogrep_best_haplogroup": get("haplogroup", "besthaplogroup", "hg"),
            "haplogrep_quality": parse_number(get("quality", "rankedquality", "score")),
            "haplogrep_missing_markers": missing, "haplogrep_private_markers": private,
            "haplogrep_n_missing_markers": len([x for x in re.split(r"[,; ]+", missing) if x]),
            "haplogrep_n_private_markers": len([x for x in re.split(r"[,; ]+", private) if x])}


def analyze_sample(sample: str, species: str, path: Path | None, markers: dict, cfg: dict, audit: list) -> tuple[dict, list]:
    vf, ms, vc, cr, hg = (cfg[x] for x in ("variant_filters", "marker_screen", "vaf_coherence", "control_region", "haplogrep"))
    usable, missing_af, low = [], 0, []
    if path:
        for v in parse_vcf(path):
            structural = v["pos"] is not None and len(v["ref"]) == len(v["alt"]) == 1 and v["ref"] in "ACGT" and v["alt"] in "ACGT" and "," not in v["alt"]
            base_ok = structural and (not vf["pass_only"] or v["filter"] in {"PASS", "."}) and v["dp"] is not None and v["dp"] >= vf["dp_min"]
            if base_ok and v["af"] is None: missing_af += 1
            if base_ok and v["af"] is not None and 0 <= v["af"] <= 1:
                usable.append(v)
                if vf["low_vaf_min"] <= v["af"] <= vf["low_vaf_max"]: low.append(v)
    # Distinct final human POS+ALT variants are the statistical unit.
    distinct: dict[tuple[int, str], dict] = {}
    for v in low:
        key = (v["pos"], v["alt"])
        if key in markers and key not in distinct: distinct[key] = v
    candidate_hits = [(k, v) for k, v in distinct.items() if ms["include_back_mutations_in_candidate_screen"] or not markers[k]["is_back_mutation"]]
    fail_hits = [(k, v) for k, v in distinct.items() if ms["include_back_mutations_in_fail_screen"] or not markers[k]["is_back_mutation"]]
    fraction = len(candidate_hits) / len(low) if low else None
    baseline = len(low) >= ms["min_low_variants_for_screen"] and len(candidate_hits) >= ms["min_human_marker_hits"] and fraction >= ms["min_fraction_low_variants_human_marker"]
    afs = [v["af"] for _, v in candidate_hits]; median = statistics.median(afs) if afs else None
    coherent_fraction = sum(abs(x - median) <= vc["tolerance"] for x in afs) / len(afs) if afs else None
    coherence = (not vc["enabled"] or (coherent_fraction is not None and coherent_fraction >= vc["min_fraction_markers_coherent"]))
    def in_cr(pos: int) -> bool: return pos >= cr["start"] or pos <= cr["end"]
    n_control = sum(in_cr(k[0]) for k, _ in candidate_hits); n_noncontrol = len(candidate_hits) - n_control
    fail_noncontrol = sum(not in_cr(k[0]) for k, _ in fail_hits)
    noncontrol_pass = fail_noncontrol >= cr["min_non_control_region_hits_for_fail"]
    classification = cfg["classification"]
    strict = ((baseline or not classification["require_baseline_marker_screen_for_fail"]) and
              (coherence or not classification["require_vaf_coherence_for_fail"]) and
              (noncontrol_pass or not classification["require_non_control_markers_for_fail"]))
    if len(low) < ms["min_low_variants_for_screen"]: status, evidence = "INSUFFICIENT_DATA", "insufficient_low_vaf_variants"
    elif baseline and strict: status, evidence = "FAIL", "baseline_marker_enrichment_with_coherent_vaf_and_non_control_signal"
    elif baseline: status, evidence = "CANDIDATE", "baseline_marker_enrichment_without_complete_strict_support"
    else: status, evidence = "PASS", "no_human_marker_enrichment"
    ad_hits = [v for _, v in candidate_hits if v["ref_depth"] is not None and v["alt_depth"] is not None]
    total_depth = sum(v["ref_depth"] + v["alt_depth"] for v in ad_hits); alt_depth = sum(v["alt_depth"] for v in ad_hits)
    row = {key: None for key in REPORT_COLUMNS}
    row.update(sample=sample, species=species, human_contamination_status=status, human_contamination_evidence=evidence,
               n_usable_variants=len(usable), n_low_variants=len(low), n_variants_missing_af=missing_af,
               n_human_marker_hits=len(candidate_hits), frac_low_variants_human_marker=fraction,
               baseline_marker_screen_pass=baseline, n_human_marker_hits_control_region=n_control,
               n_human_marker_hits_non_control_region=n_noncontrol,
               frac_human_marker_hits_control_region=n_control/len(candidate_hits) if candidate_hits else None,
               non_control_marker_pass=noncontrol_pass, median_human_marker_af=median,
               mean_human_marker_af=statistics.mean(afs) if afs else None, min_human_marker_af=min(afs) if afs else None,
               max_human_marker_af=max(afs) if afs else None,
               human_marker_af_mad=statistics.median(abs(x-median) for x in afs) if afs else None,
               human_marker_af_iqr=quantile(afs,.75)-quantile(afs,.25) if afs else None,
               fraction_human_markers_near_median_af=coherent_fraction, vaf_coherence_pass=coherence,
               estimated_human_fraction_median_af=median,
               estimated_human_mtDNA_fraction=alt_depth/total_depth if total_depth else None,
               n_human_markers_with_ad=len(ad_hits), human_marker_total_depth=total_depth if ad_hits else None,
               human_marker_total_alt_depth=alt_depth if ad_hits else None,
               human_marker_list=";".join(markers[k]["marker"] for k,_ in candidate_hits),
               n_back_mutation_marker_hits=sum(markers[k]["is_back_mutation"] for k,_ in distinct.items()))
    # HaploGrep-selected markers remain marker-only and have their own VAF/back-mutation policy.
    selected = [(k,v) for k,v in distinct.items() if hg["input_vaf_min"] <= v["af"] <= hg["input_vaf_max"] and (not hg["exclude_back_mutations"] or not markers[k]["is_back_mutation"])]
    if not hg["require_phylotree_marker"]:  # all low variants is supported, but never the complete VCF
        selected = [((v["pos"],v["alt"]),v) for v in low if hg["input_vaf_min"] <= v["af"] <= hg["input_vaf_max"]]
    selected_keys = {k for k,_ in selected}; candidate_keys={k for k,_ in candidate_hits}; fail_keys={k for k,_ in fail_hits}
    for k, v in distinct.items():
        m=markers[k]; audit.append(dict(sample=sample,species=species,human_chrom=v["chrom"],human_pos=v["pos"],human_ref=v["ref"],human_alt=v["alt"],dp=v["dp"],af=v["af"],af_source=v["af_source"],ref_depth=v["ref_depth"],alt_depth=v["alt_depth"],liftover_allele_status=v["liftover_allele_status"],human_marker=m["marker"],is_back_mutation=m["is_back_mutation"],in_control_region=in_cr(v["pos"]),used_in_candidate_screen=k in candidate_keys,used_in_fail_screen=k in fail_keys,used_in_haplogrep=k in selected_keys))
    return row, selected


def write_tsv(path: Path, columns: list[str], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=columns,delimiter="\t",extrasaction="ignore");writer.writeheader()
        for row in rows: writer.writerow({k: ("NA" if v is None else boolean(v) if isinstance(v,bool) else v) for k,v in row.items()})


def main() -> int:
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("--config",type=Path,required=True);ap.add_argument("--validate-inputs",action="store_true");ap.add_argument("--overwrite",action="store_true");args=ap.parse_args()
    config=read_simple_yaml(args.config); sec=config.get("human_contamination") or {}
    if sec.get("enabled",True) is False: print("[human_contamination] disabled; skipping"); return 0
    paths=sec["paths"]; sample_file=resolve(paths["sample_ref_file"]); marker_file=resolve(paths["phylotree_marker_file"]); vcf_dir=resolve(paths["input_vcf_dir"]); out=resolve(paths["output_dir"])
    if not sample_file.is_file(): raise ValueError(f"sample metadata missing: {sample_file}")
    if not marker_file.is_file(): raise ValueError(f"PhyloTree marker file missing: {marker_file}")
    samples=load_samples(sample_file); markers,marker_qc=load_markers(marker_file)
    vcfs={s:find_vcf(vcf_dir,paths["input_vcf_pattern"],s) if vcf_dir.is_dir() else None for s,_ in samples}
    hg=sec["haplogrep"]; tool,mode=haplogrep_tool(hg); available=tool is not None
    validation={"samples expected":len(samples),"lifted VCFs found":sum(v is not None for v in vcfs.values()),"lifted VCFs missing":sum(v is None for v in vcfs.values()),"PhyloTree markers raw":marker_qc["raw"],"PhyloTree simple SNVs":marker_qc["simple"],"PhyloTree back-mutation SNVs":marker_qc["back"],"PhyloTree duplicate POS/ALT removed":marker_qc["duplicates"],"HaploGrep enabled":boolean(hg["enabled"]),"HaploGrep executable mode":mode,"HaploGrep tool available":boolean(available)}
    if args.validate_inputs:
        for k,v in validation.items(): print(f"{k}: {v}")
        if not vcf_dir.is_dir(): raise ValueError(f"input VCF directory missing: {vcf_dir}")
        if hg["enabled"] and hg["require_tool_when_enabled"] and not available: raise ValueError("HaploGrep enabled but configured tool is unavailable")
        out.mkdir(parents=True,exist_ok=True); test=out/".write_test";test.write_text("");test.unlink(); return 0 if all(vcfs.values()) else 1
    reports=out/"reports"
    if reports.exists() and any(reports.iterdir()) and not args.overwrite: raise ValueError(f"outputs already exist; use --overwrite: {reports}")
    if args.overwrite and out.exists(): shutil.rmtree(out)
    reports.mkdir(parents=True); input_dir=resolve(hg["input_dir"]); output_dir=resolve(hg["output_dir"]);input_dir.mkdir(parents=True,exist_ok=True);output_dir.mkdir(parents=True,exist_ok=True)
    rows=[];audit=[]
    for sample,species in samples:
        row,selected=analyze_sample(sample,species,vcfs[sample],markers,sec,audit); n=len(selected)
        row.update(haplogrep_status="NOT_RUN",haplogrep_support_status="NOT_RUN",haplogrep_n_input_markers=n,haplogrep_sparse_input=n<hg["sparse_input_marker_threshold"])
        mode_run=hg["run_mode"]; eligible=hg["enabled"] and mode_run!="disabled" and (mode_run=="all" or mode_run=="all_with_min_markers" and n>=hg["min_input_markers"] or mode_run=="candidate_only" and row["baseline_marker_screen_pass"])
        if eligible and n<hg["min_input_markers"]: row.update(haplogrep_status="INSUFFICIENT_MARKERS",haplogrep_support_status="INSUFFICIENT_MARKERS")
        elif eligible and not available:
            row.update(haplogrep_status="TOOL_UNAVAILABLE",haplogrep_support_status="TOOL_UNAVAILABLE")
            if hg["require_tool_when_enabled"]: raise ValueError("HaploGrep tool unavailable")
        elif eligible:
            profile=input_dir/f"{sample}.human_contaminant.hsd"; audit_path=input_dir/f"{sample}.human_contaminant.audit.tsv"
            profile.write_text("SampleID\tRange\tHaplogroup\tPolymorphisms\n"+f"{sample}\t1-16569\t?\t"+" ".join(f"{k[0]}{k[1]}" for k,_ in selected)+"\n")
            write_tsv(audit_path,["sample","human_pos","human_ref","human_alt","af","dp","ad","marker","is_back_mutation"],[dict(sample=sample,human_pos=k[0],human_ref=v["ref"],human_alt=k[1],af=v["af"],dp=v["dp"],ad=f'{v["ref_depth"]},{v["alt_depth"]}' if v["ref_depth"] is not None else "NA",marker=markers.get(k,{}).get("marker",f"{k[0]}{k[1]}"),is_back_mutation=markers.get(k,{}).get("is_back_mutation",False)) for k,v in selected])
            raw=output_dir/f"{sample}.haplogrep.tsv"; command=tool+["classify","--in",str(profile),"--out",str(raw),"--tree",str(hg["tree"])]
            try:
                subprocess.run(command,check=True,capture_output=True,text=True); parsed=parse_haplogrep_output(raw);row.update(parsed);quality=parsed.get("haplogrep_quality");support="SUPPORTIVE" if quality is not None and quality>=hg["min_quality_for_support"] else "LOW_INFORMATION";row.update(haplogrep_status="COMPLETED",haplogrep_support_status=support)
            except subprocess.CalledProcessError as exc:
                (output_dir/f"{sample}.stderr.txt").write_text(exc.stderr or "");row.update(haplogrep_status="FAILED",haplogrep_support_status="LOW_INFORMATION")
            if not hg["keep_input_files"]: profile.unlink(missing_ok=True);audit_path.unlink(missing_ok=True)
            if not hg["keep_raw_output"]: raw.unlink(missing_ok=True)
        rows.append(row)
    write_tsv(reports/"human_contamination_report.tsv",REPORT_COLUMNS,rows)
    with gzip.open(reports/"human_marker_overlap_variants.tsv.gz","wt",encoding="utf-8",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=AUDIT_COLUMNS,delimiter="\t");writer.writeheader()
        for r in audit:writer.writerow({k:boolean(v) if isinstance(v,bool) else "NA" if v is None else v for k,v in r.items()})
    metrics={"n_samples_total":len(rows),**{f"n_{s}":sum(r["human_contamination_status"]==s for r in rows) for s in ("PASS","CANDIDATE","FAIL","INSUFFICIENT_DATA")},"n_baseline_marker_screen_pass":sum(r["baseline_marker_screen_pass"] for r in rows),"n_vaf_coherence_pass":sum(r["vaf_coherence_pass"] for r in rows),"n_haplogrep_run":sum(r["haplogrep_status"]=="COMPLETED" for r in rows),"n_haplogrep_supportive":sum(r["haplogrep_support_status"] in {"SUPPORTIVE","HIGH_SUPPORT"} for r in rows),"n_haplogrep_low_information":sum(r["haplogrep_support_status"]=="LOW_INFORMATION" for r in rows)}
    with (reports/"human_contamination_summary.tsv").open("w",newline="") as h:w=csv.writer(h,delimiter="\t");w.writerow(("metric","value"));w.writerows(metrics.items());w.writerows((f"threshold.{section}.{k}",v) for section in ("variant_filters","marker_screen","vaf_coherence","control_region") for k,v in sec[section].items())
    with (reports/"human_contamination_run_parameters.tsv").open("w",newline="") as h:w=csv.writer(h,delimiter="\t");w.writerow(("parameter","value"));w.writerow(("effective_config_json",json.dumps(sec,sort_keys=True,separators=(",",":"))));w.writerows((f"marker_qc.{k}",v) for k,v in marker_qc.items())
    print(f"[human_contamination] samples={len(rows)} output={reports}");return 0

if __name__=="__main__":
    try: raise SystemExit(main())
    except (KeyError,ValueError,OSError,subprocess.SubprocessError) as exc: print(f"ERROR: {exc}",file=sys.stderr);raise SystemExit(2)
