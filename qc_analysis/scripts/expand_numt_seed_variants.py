#!/usr/bin/env python3
"""Limited AF expansion around REMOVE-level NUMT-supported HET clusters."""
from __future__ import annotations

import argparse, csv, math, sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from qc_analysis.lib.simple_yaml import read_simple_yaml


def resolve(x):
    p = Path(str(x)).expanduser()
    return p if p.is_absolute() else ROOT / p


def read_tsv(p):
    with Path(p).open(newline="", encoding="utf-8") as h:
        return list(csv.DictReader(h, delimiter="\t"))


def write_tsv(p, rows, fields):
    Path(p).parent.mkdir(parents=True, exist_ok=True)
    with Path(p).open("w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)


def fnum(x):
    try:
        v = float(x); return v if math.isfinite(v) else None
    except (TypeError, ValueError): return None


def inum(x):
    try: return int(float(x))
    except (TypeError, ValueError): return None


def yes(x): return str(x).strip().upper() in {"YES", "TRUE", "T", "1"}


def circ_dist(a, b, L):
    d = abs(a - b); return min(d, L - d)


def inside_cluster(pos, start, end, wraps):
    return (pos >= start or pos <= end) if wraps else start <= pos <= end


def cluster_dist(pos, start, end, wraps, L):
    if inside_cluster(pos, start, end, wraps): return 0
    return min(circ_dist(pos, start, L), circ_dist(pos, end, L))


def inside_interval(pos, start, end):
    return start <= pos <= end if start <= end else (pos >= start or pos <= end)


def strong_seeds(cluster_rows):
    out = []
    for r in cluster_rows:
        if str(r.get("filter_action", "")).upper() != "REMOVE": continue
        sample_ok = str(r.get("sample_numt_evidence", "")).upper() == "YES"
        species_rec = str(r.get("species_numt_evidence", "")).upper() == "YES" and str(r.get("recurrence", "")).upper() == "YES"
        if not (sample_ok or species_rec): continue
        vals = [inum(r.get(k)) for k in ("cluster_start", "cluster_end", "numt_chrm_start", "numt_chrm_end")]
        afs = [fnum(r.get(k)) for k in ("median_af", "af_min", "af_max")]
        if any(v is None for v in vals + afs): continue
        out.append({
            "sample": str(r.get("sample", "")), "id": str(r.get("cluster_id", "")),
            "scope": "SAMPLE" if sample_ok else "SPECIES", "tier": str(r.get("numt_tier", "NONE")),
            "start": vals[0], "end": vals[1], "wraps": yes(r.get("wraps_origin")),
            "numt_start": vals[2], "numt_end": vals[3],
            "median": afs[0], "af_min": afs[1], "af_max": afs[2], "row": r,
        })
    return out


def reset_previous(variants, clusters):
    for r in variants:
        if str(r.get("numt_expanded", "")).upper() != "YES": continue
        r["filter_action"] = r.get("pre_expansion_filter_action", "KEEP") or "KEEP"
        r["filter_reason"] = r.get("pre_expansion_filter_reason", "")
        r["numt_scope"] = r.get("pre_expansion_numt_scope", "NONE") or "NONE"
        r["numt_tier"] = r.get("pre_expansion_numt_tier", "NONE") or "NONE"
        for k in ("numt_expansion_seed_id", "numt_expansion_delta_af", "numt_expansion_distance_bp", "numt_expansion_final_af_span"):
            r[k] = ""
        r["numt_expanded"] = "NO"
    for r in clusters:
        r["numt_expansion_n_variants"] = 0
        r["numt_expansion_final_af_min"] = r.get("af_min", "")
        r["numt_expansion_final_af_max"] = r.get("af_max", "")
        r["numt_expansion_final_af_span"] = r.get("af_span", "")


def apply_expansion(clusters, variants, mt_length=16569, max_distance_bp=250, max_delta_af=0.07, max_total_af_span=0.10):
    """Expand strict REMOVE seeds without changing the original seed centre."""
    reset_previous(variants, clusters)
    seeds = strong_seeds(clusters)
    by_sample = defaultdict(list)
    for s in seeds: by_sample[s["sample"]].append(s)

    assigned = defaultdict(list); diagnostics = []
    for i, r in enumerate(variants):
        if str(r.get("filter_action", "")).upper() == "REMOVE": continue
        sample, pos, af = str(r.get("sample", "")), inum(r.get("source_pos")), fnum(r.get("source_af"))
        if pos is None or af is None: continue
        nearby = []
        for s in by_sample.get(sample, []):
            dist = cluster_dist(pos, s["start"], s["end"], s["wraps"], mt_length)
            if dist > max_distance_bp: continue
            same = inside_interval(pos, s["numt_start"], s["numt_end"])
            delta = abs(af - s["median"])
            nearby.append((0 if same else 1, dist, delta, s["id"], same, s))
        if not nearby: continue
        _, dist, delta, _, same, s = min(nearby, key=lambda x: (x[0], x[1], x[2], x[3]))
        if not same or delta > max_delta_af + 1e-12: continue
        d = {"i": i, "key": "|".join(str(r.get(k, "")) for k in ("sample", "source_pos", "source_ref", "source_alt")),
             "sample": sample, "pos": pos, "ref": str(r.get("source_ref", "")), "alt": str(r.get("source_alt", "")),
             "af": af, "seed": s, "distance": dist, "delta": delta}
        assigned[s["id"]].append(d); diagnostics.append(d)

    accepted = set(); seed_map = {s["id"]: s for s in seeds}
    for sid in sorted(assigned):
        s = seed_map[sid]; lo, hi = s["af_min"], s["af_max"]; kept = []
        for d in sorted(assigned[sid], key=lambda x: (x["delta"], x["distance"], x["pos"], x["ref"], x["alt"])):
            nlo, nhi = min(lo, d["af"]), max(hi, d["af"])
            if nhi - nlo > max_total_af_span + 1e-12: continue
            lo, hi = nlo, nhi; kept.append(d); accepted.add(d["key"])
            r = variants[d["i"]]
            if not r.get("pre_expansion_filter_action"):
                r["pre_expansion_filter_action"] = r.get("filter_action", "KEEP") or "KEEP"
                r["pre_expansion_filter_reason"] = r.get("filter_reason", "")
                r["pre_expansion_numt_scope"] = r.get("numt_scope", "NONE") or "NONE"
                r["pre_expansion_numt_tier"] = r.get("numt_tier", "NONE") or "NONE"
            r.update({"numt_expanded": "YES", "numt_expansion_seed_id": sid,
                      "numt_expansion_delta_af": f"{d['delta']:.6g}", "numt_expansion_distance_bp": d["distance"],
                      "numt_expansion_final_af_span": f"{nhi-nlo:.6g}", "filter_action": "REMOVE",
                      "filter_reason": "NUMT_SEED_EXPANSION", "numt_scope": s["scope"], "numt_tier": s["tier"]})
            if not str(r.get("cluster_class", "")): r["cluster_class"] = "NUMT_EXPANDED"
        s["row"].update({"numt_expansion_n_variants": len(kept), "numt_expansion_final_af_min": f"{lo:.6g}",
                         "numt_expansion_final_af_max": f"{hi:.6g}", "numt_expansion_final_af_span": f"{hi-lo:.6g}"})

    out = []
    for d in diagnostics:
        s = d["seed"]; ok = d["key"] in accepted
        out.append({"variant_key": d["key"], "sample": d["sample"], "source_pos": d["pos"], "source_ref": d["ref"],
                    "source_alt": d["alt"], "source_af": f"{d['af']:.6g}", "seed_id": s["id"], "seed_scope": s["scope"],
                    "seed_median_af": f"{s['median']:.6g}", "seed_af_min": f"{s['af_min']:.6g}", "seed_af_max": f"{s['af_max']:.6g}",
                    "numt_chrm_start": s["numt_start"], "numt_chrm_end": s["numt_end"], "distance_bp": d["distance"],
                    "delta_af": f"{d['delta']:.6g}", "accepted": "YES" if ok else "NO",
                    "rejection_reason": "" if ok else "EXPANDED_AF_SPAN_CAP"})
    return out


def rebuild_removal(variants):
    out, seen = [], set()
    for r in variants:
        if str(r.get("filter_action", "")).upper() != "REMOVE": continue
        k = tuple(r.get(x, "") for x in ("sample", "source_chrom", "source_pos", "source_ref", "source_alt"))
        if k in seen: continue
        seen.add(k); out.append(dict(r))
    return out


def update_samples(samples, variants):
    rem, exp = Counter(), Counter()
    for r in variants:
        s = str(r.get("sample", ""))
        if str(r.get("filter_action", "")).upper() == "REMOVE": rem[s] += 1
        if str(r.get("numt_expanded", "")).upper() == "YES": exp[s] += 1
    out = []
    for raw in samples:
        r = dict(raw); s = str(r.get("sample", "")); n = inum(r.get("n_het")) or 0
        r["n_numt_variants_to_remove"] = rem[s]; r["n_numt_expanded_variants"] = exp[s]
        r["fraction_het_removed"] = f"{rem[s]/n:.6g}" if n else "0"
        if rem[s] > 0: r["numt_sample"] = "YES"
        out.append(r)
    return out


def add_fields(fields, extras):
    out = list(fields)
    for x in extras:
        if x not in out: out.append(x)
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--config", type=Path, required=True); args = ap.parse_args()
    sec = read_simple_yaml(args.config).get("local_heteroplasmy_qc", {})
    if sec.get("numt_seed_expansion_enabled", True) is False:
        print("[numt_seed_expansion] disabled", file=sys.stderr); return 0
    report = resolve(sec.get("output_dir", "results/qc/local_heteroplasmy_qc")) / "reports"
    cp, vp, sp = report/"local_heteroplasmy_cluster_summary.tsv", report/"local_heteroplasmy_variant_detail.tsv", report/"local_heteroplasmy_sample_summary.tsv"
    clusters, variants, samples = read_tsv(cp), read_tsv(vp), read_tsv(sp)
    L = int(sec.get("mt_length", 16569)); dist = int(sec.get("numt_seed_expansion_distance_bp", 250))
    delta = float(sec.get("numt_seed_expansion_max_delta_af", 0.07)); cap = float(sec.get("numt_seed_expansion_max_total_af_span", 0.10))
    diag = apply_expansion(clusters, variants, L, dist, delta, cap)
    removals, samples = rebuild_removal(variants), update_samples(samples, variants)
    numt_samples = [r for r in samples if str(r.get("numt_sample", "")).upper() == "YES"]

    vf = add_fields(list(variants[0]) if variants else [], ["pre_expansion_filter_action", "pre_expansion_filter_reason", "pre_expansion_numt_scope", "pre_expansion_numt_tier",
                    "numt_expanded", "numt_expansion_seed_id", "numt_expansion_delta_af", "numt_expansion_distance_bp", "numt_expansion_final_af_span"])
    cf = add_fields(list(clusters[0]) if clusters else [], ["numt_expansion_n_variants", "numt_expansion_final_af_min", "numt_expansion_final_af_max", "numt_expansion_final_af_span"])
    sf = add_fields(list(samples[0]) if samples else [], ["n_numt_expanded_variants"])
    rf = ["sample", "species", "source_chrom", "source_pos", "source_ref", "source_alt", "source_af", "source_dp", "clustered", "cluster_id", "cluster_class",
          "numt_scope", "numt_tier", "recurrence", "numt_expanded", "numt_expansion_seed_id", "numt_expansion_delta_af", "numt_expansion_distance_bp",
          "numt_expansion_final_af_span", "filter_action", "filter_reason"]
    write_tsv(cp, clusters, cf); write_tsv(vp, variants, vf); write_tsv(sp, samples, sf)
    write_tsv(report/"numt_samples.tsv", numt_samples, sf); write_tsv(report/"numt_variants_to_remove.tsv", removals, rf)
    df = ["variant_key", "sample", "source_pos", "source_ref", "source_alt", "source_af", "seed_id", "seed_scope", "seed_median_af", "seed_af_min", "seed_af_max",
          "numt_chrm_start", "numt_chrm_end", "distance_bp", "delta_af", "accepted", "rejection_reason"]
    write_tsv(report/"numt_seed_expansion_variants.tsv", diag, df)
    acc = [r for r in diag if r["accepted"] == "YES"]; scopes = Counter(r["seed_scope"] for r in acc)
    summary = [{"max_distance_bp": dist, "max_delta_af": f"{delta:.6g}", "max_total_af_span": f"{cap:.6g}", "n_remove_level_seeds": len(strong_seeds(clusters)),
                "n_expansion_candidates": len(diag), "n_expanded_variants": len(acc), "n_expanded_samples": len({r['sample'] for r in acc}),
                "n_sample_scope_expanded_variants": scopes["SAMPLE"], "n_species_scope_expanded_variants": scopes["SPECIES"], "n_total_variants_to_remove": len(removals)}]
    write_tsv(report/"numt_seed_expansion_summary.tsv", summary, list(summary[0]))
    print(f"[numt_seed_expansion] seeds={summary[0]['n_remove_level_seeds']} candidates={len(diag)} expanded={len(acc)} samples={summary[0]['n_expanded_samples']} total_remove={len(removals)}", file=sys.stderr)
    return 0


if __name__ == "__main__": raise SystemExit(main())
