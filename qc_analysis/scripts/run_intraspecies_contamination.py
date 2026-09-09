#!/usr/bin/env python3
"""Validated original-coordinate intra-species contamination analysis.

Target and source cohorts are deliberately asymmetric:
- TARGETS are deduplicated samples; contamination is tested only when sample QC
  is PASS and the sample has at least ``min_target_het`` strict HET calls.
- SOURCES are *all* deduplicated same-species samples, irrespective of sample QC.

Low-A target evidence uses the original VCF calls (before local-cluster removal)
so local artifact filtering cannot erase a real contamination signal.

Genome-wide dispersion metrics are diagnostic only. If a downstream local-
heteroplasmy report already exists, the script also recomputes source matching
after excluding detected local-cluster variants. This ALL-vs-NONCLUSTER
comparison is a sensitivity analysis only and does not yet alter the primary
contamination classification.
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

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from qc_analysis.lib.simple_yaml import read_simple_yaml

REPORT_COLUMNS = """sample species sample_qc_status n_strict_het target_eligible target_ineligible_reason n_species_samples n_source_candidates n_usable_variants n_lowA best_source_sample best_source_qc_status best_overlap best_frac_lowA_in_highB best_overlap_positions best_overlap_occupied_bins best_overlap_linear_span_bp best_overlap_circular_span_bp best_overlap_max_in_window best_overlap_max_local_fraction local_cluster_sensitivity_available n_lowA_clustered n_lowA_noncluster best_source_sample_noncluster best_overlap_noncluster best_frac_lowA_in_highB_noncluster best_overlap_positions_noncluster best_overlap_occupied_bins_noncluster best_overlap_circular_span_bp_noncluster best_overlap_max_local_fraction_noncluster n_overlap_from_cluster overlap_retention_after_cluster_removal source_stable_after_cluster_removal n_anchor_pool_excluding_A n_anchor_tested_in_A n_depressed_anchor mt_high_hets_contamination mt_high_hets_mode anchor_evidence_level anchor_source_count n_mirror_pairs n_low_variants_with_mirror mirror_low_fraction normalized_mirror_support mirror_p95_threshold mirror_p99_threshold mirror_calibration_status mirror_support_candidate mirror_support_highconf contamination_status contamination_flag_candidate contamination_flag_highconf qc_status qc_reason""".split()

ELIGIBILITY_COLUMNS = """sample species sample_qc_status n_strict_het min_target_het target_eligible target_ineligible_reason""".split()

DEFAULTS = dict(
    dp_min=100,
    low_vaf_min=.01,
    low_vaf_max=.20,
    min_alt_reads=3,
    source_high_vaf_min=.90,
    target_het_af_min=.10,
    target_het_af_max=.95,
    min_target_het=4,
    mt_length=16569,
    dispersion_bin_bp=1000,
    dispersion_window_bp=1000,
    mt_lower=.80,
    mt_depressed_upper=.998,
    mt_anchor_upper=1.,
    min_n_lowA=5,
    min_overlap=3,
    min_frac_lowA_in_highB_candidate=.50,
    min_frac_lowA_in_highB_highconf=.6213636363636358,
    contam_threshold_candidate=.036420574377757434,
    contam_threshold_highconf=.07103935483870959,
    mirror_low_vaf_min=.01,
    mirror_low_vaf_max=.20,
    mirror_high_vaf_min=.80,
    mirror_high_vaf_max=.998,
    mirror_tolerance=0.,
    min_negative_control_values=3,
    target_negative_control_tier="tier2_location_and_batch_different",
)


def path(v):
    p = Path(str(v)).expanduser()
    return p if p.is_absolute() else ROOT / p


def truth(v, d=False):
    if v is None:
        return d
    if isinstance(v, bool):
        return v
    raise ValueError(f"expected boolean, got {v!r}")


def load_rows(p):
    with p.open(newline="", encoding="utf-8") as h:
        return list(csv.DictReader(h, delimiter="\t"))


def write_rows(p, rows, fields):
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def as_float(value, default=None):
    try:
        x = float(value)
        return x if math.isfinite(x) else default
    except (TypeError, ValueError):
        return default


def num(r, k):
    return float(r[k])


def key(r):
    return r["CHROM"], r["POS"], r["REF"], r["ALT"]


def read_sample_pairs(p: Path) -> set[tuple[str, str]]:
    if not p.is_file():
        raise FileNotFoundError(f"sample list not found: {p}")
    with p.open(encoding="utf-8") as h:
        lines = [x.rstrip("\n") for x in h if x.strip() and not x.lstrip().startswith("#")]
    if not lines:
        raise ValueError(f"sample list is empty: {p}")
    first = lines[0].split("\t") if "\t" in lines[0] else lines[0].split()
    lower = [x.strip().lower() for x in first]
    if "sample" in lower and "species" in lower:
        if "\t" in lines[0]:
            rows = list(csv.DictReader(lines, delimiter="\t"))
        else:
            names = first
            rows = [dict(zip(names, x.split())) for x in lines[1:]]
        out = set()
        for r in rows:
            s = str(r.get("sample") or r.get("Sample") or "").strip()
            sp = str(r.get("species") or r.get("Species") or "").strip()
            if s and sp:
                out.add((sp, s))
        return out
    out = set()
    for line in lines:
        f = line.split()
        if len(f) >= 2:
            out.add((f[1], f[0]))
    if not out:
        raise ValueError(f"could not parse sample/species pairs from: {p}")
    return out


def read_sample_qc(p: Path) -> dict[str, str]:
    if not p.is_file():
        raise FileNotFoundError(f"sample QC report not found: {p}")
    rows = load_rows(p)
    if rows and ("sample" not in rows[0] or "qc_status" not in rows[0]):
        raise ValueError(f"sample QC report must contain sample and qc_status: {p}")
    return {
        str(r.get("sample", "")).strip(): str(r.get("qc_status", "")).strip().upper()
        for r in rows if str(r.get("sample", "")).strip()
    }


def count_strict_hets(source_qc_path: Path, allowed_samples: set[str], p: dict) -> Counter:
    if not source_qc_path.is_file():
        raise FileNotFoundError(
            f"source variant QC report not found: {source_qc_path}. "
            "Run pre_liftover_variant_qc before intraspecies_contamination."
        )
    counts = Counter()
    lo = float(p["target_het_af_min"])
    hi = float(p["target_het_af_max"])
    dp_min = float(p["dp_min"])
    for r in load_rows(source_qc_path):
        sample = str(r.get("sample", "")).strip()
        if sample not in allowed_samples:
            continue
        af = as_float(r.get("source_af"))
        dp = as_float(r.get("source_dp"))
        if (
            str(r.get("source_call_class", "")).upper() == "HET"
            and str(r.get("source_filter", "")).upper() == "PASS"
            and str(r.get("source_variant_qc", "")).upper() == "PASS"
            and af is not None and lo <= af < hi
            and dp is not None and dp >= dp_min
        ):
            counts[sample] += 1
    return counts


def load_clustered_variant_keys(cluster_report: Path | None):
    """Return sample -> exact SOURCE allele keys for all detected local clusters.

    The report contains only the strict HET range evaluated by local clustering,
    so low-A calls below 0.10 are necessarily retained in the NONCLUSTER branch.
    """
    if cluster_report is None or not cluster_report.is_file():
        return {}, False
    rows = load_rows(cluster_report)
    if rows and not {"sample", "source_chrom", "source_pos", "source_ref", "source_alt", "clustered"}.issubset(rows[0]):
        raise ValueError(f"invalid local heteroplasmy variant report: {cluster_report}")
    out = defaultdict(set)
    for r in rows:
        if str(r.get("clustered", "")).upper() != "YES":
            continue
        sample = str(r.get("sample", "")).strip()
        if not sample:
            continue
        out[sample].add((
            str(r.get("source_chrom", "")),
            str(r.get("source_pos", "")),
            str(r.get("source_ref", "")),
            str(r.get("source_alt", "")),
        ))
    return dict(out), True


def quantile7(values, p):
    x = sorted(values)
    h = (len(x) - 1) * p
    i = math.floor(h)
    return x[i] + (x[min(i + 1, len(x) - 1)] - x[i]) * (h - i)


def low_a_rows(rows, p):
    out = []
    for r in rows:
        af = as_float(r.get("VAF"))
        ad_alt = as_float(r.get("AD_alt"))
        if af is None or not (float(p["low_vaf_min"]) <= af <= float(p["low_vaf_max"])):
            continue
        if ad_alt is None or ad_alt < int(p["min_alt_reads"]):
            continue
        out.append(r)
    return out


def overlap_dispersion(overlap_keys, p, suffix=""):
    positions = sorted({int(k[1]) for k in overlap_keys})
    n = len(positions)
    prefix = "best_overlap_"
    def name(x):
        return f"{prefix}{x}{suffix}"
    if not positions:
        return {
            name("positions"): "",
            name("occupied_bins"): 0,
            name("linear_span_bp"): 0,
            name("circular_span_bp"): 0,
            name("max_in_window"): 0,
            name("max_local_fraction"): None,
        }
    mt_length = int(p["mt_length"])
    bin_bp = int(p["dispersion_bin_bp"])
    window_bp = int(p["dispersion_window_bp"])
    if mt_length <= 0 or bin_bp <= 0 or window_bp <= 0:
        raise ValueError("mt_length, dispersion_bin_bp and dispersion_window_bp must be positive")
    occupied_bins = len({(pos - 1) // bin_bp for pos in positions})
    linear_span = positions[-1] - positions[0] if n > 1 else 0
    if n > 1:
        gaps = [positions[i + 1] - positions[i] for i in range(n - 1)]
        gaps.append((positions[0] + mt_length) - positions[-1])
        circular_span = mt_length - max(gaps)
    else:
        circular_span = 0
    doubled = positions + [x + mt_length for x in positions]
    max_in_window = 0
    j = 0
    for i in range(n):
        if j < i:
            j = i
        limit = doubled[i] + window_bp
        while j < i + n and doubled[j] <= limit:
            j += 1
        max_in_window = max(max_in_window, j - i)
    return {
        name("positions"): ",".join(str(x) for x in positions),
        name("occupied_bins"): occupied_bins,
        name("linear_span_bp"): linear_span,
        name("circular_span_bp"): circular_span,
        name("max_in_window"): max_in_window,
        name("max_local_fraction"): max_in_window / n,
    }


def best_source_match(lowkeys, otherhigh):
    ranked = [(len(lowkeys & ks), s) for s, ks in otherhigh.items()]
    best_overlap, best_source = max(ranked, key=lambda x: (x[0], x[1])) if ranked else (0, "")
    best_keys = lowkeys & otherhigh.get(best_source, set()) if best_source else set()
    frac = best_overlap / len(lowkeys) if lowkeys else None
    return best_source, best_overlap, frac, best_keys


def mirror_stats(rows, p):
    lows = [r for r in rows if p["mirror_low_vaf_min"] <= num(r, "VAF") <= p["mirror_low_vaf_max"]]
    highs = [r for r in rows if p["mirror_high_vaf_min"] <= num(r, "VAF") <= p["mirror_high_vaf_max"]]
    pairs = [(i, j) for i, a in enumerate(lows) for j, b in enumerate(highs)
             if abs(num(a, "VAF") + num(b, "VAF") - 1) <= p["mirror_tolerance"] + 1e-12]
    mirrored = len({i for i, _ in pairs})
    frac = mirrored / len(lows) if lows else 0.
    normalized = len(pairs) / (len(lows) * len(highs)) if lows and highs else 0.
    return len(pairs), mirrored, frac, normalized


def calibration(rows, p, nc_path):
    if not nc_path:
        return None, None, "not_calibrated_no_file", 0
    q = path(nc_path)
    if not q.is_file():
        return None, None, "not_calibrated_missing_file", 0
    table = load_rows(q)
    tier = [r for r in table if r.get("negative_control_tier") == str(p["target_negative_control_tier"])]
    if not tier:
        return None, None, "not_calibrated_no_tier2_pairs", 0
    by = defaultdict(list)
    for r in rows:
        by[r["Sample"]].append(r)
    vals = []
    for pair in tier:
        a = by.get(pair.get("Sample_A", ""), [])
        b = by.get(pair.get("Sample_B", ""), [])
        if a and b:
            vals.append(mirror_stats(a + b, p)[3])
    if len(vals) < int(p["min_negative_control_values"]):
        return None, None, "not_calibrated_insufficient_values", len(vals)
    return quantile7(vals, .95), quantile7(vals, .99), "calibrated", len(vals)


def analyse(rows, p, source_pairs, qc_status, het_counts, clustered_keys=None,
            cluster_sensitivity_available=False, negative_control_pairs=None):
    source_pairs = set(source_pairs)
    clustered_keys = clustered_keys or {}
    source_samples = {s for _, s in source_pairs}
    usable = [
        r for r in rows
        if r.get("Sample") in source_samples
        and num(r, "DP") >= float(p["dp_min"])
        and (not p["use_snv_only"] or r["Type"] == "SNV")
        and (not p["pass_only"] or r["FILTER"] == "PASS")
    ]
    by = defaultdict(list)
    for r in usable:
        by[r["Species"], r["Sample"]].append(r)

    p95, p99, cal_status, _ = calibration(usable, p, negative_control_pairs)
    counts = Counter(sp for sp, _ in source_pairs)
    source_high_min = float(p["source_high_vaf_min"])
    min_target_het = int(p["min_target_het"])
    out = []
    eligibility = []

    for species, sample in sorted(source_pairs):
        sample_qc = qc_status.get(sample, "MISSING")
        n_het = int(het_counts.get(sample, 0))
        eligible = sample_qc == "PASS" and n_het >= min_target_het
        if sample_qc != "PASS":
            ineligible_reason = "sample_qc_not_pass"
        elif n_het < min_target_het:
            ineligible_reason = "insufficient_strict_heteroplasmy"
        else:
            ineligible_reason = ""
        eligibility.append(dict(
            sample=sample, species=species, sample_qc_status=sample_qc,
            n_strict_het=n_het, min_target_het=min_target_het,
            target_eligible=eligible, target_ineligible_reason=ineligible_reason,
        ))

        own = by[species, sample]
        others = [s for sp, s in source_pairs if sp == species and s != sample]
        base = dict(
            sample=sample, species=species, sample_qc_status=sample_qc,
            n_strict_het=n_het, target_eligible=eligible,
            target_ineligible_reason=ineligible_reason,
            n_species_samples=counts[species], n_source_candidates=len(others),
            n_usable_variants=len(own), n_lowA=0, best_source_sample="",
            best_source_qc_status="", best_overlap=0,
            best_frac_lowA_in_highB=None,
            best_overlap_positions="", best_overlap_occupied_bins=0,
            best_overlap_linear_span_bp=0, best_overlap_circular_span_bp=0,
            best_overlap_max_in_window=0, best_overlap_max_local_fraction=None,
            local_cluster_sensitivity_available=cluster_sensitivity_available,
            n_lowA_clustered=0, n_lowA_noncluster=0,
            best_source_sample_noncluster="", best_overlap_noncluster=0,
            best_frac_lowA_in_highB_noncluster=None,
            best_overlap_positions_noncluster="", best_overlap_occupied_bins_noncluster=0,
            best_overlap_circular_span_bp_noncluster=0,
            best_overlap_max_local_fraction_noncluster=None,
            n_overlap_from_cluster=0, overlap_retention_after_cluster_removal=None,
            source_stable_after_cluster_removal="NA",
            n_anchor_pool_excluding_A=0, n_anchor_tested_in_A=0,
            n_depressed_anchor=0, mt_high_hets_contamination=None,
            mt_high_hets_mode="not_tested", anchor_evidence_level="not_tested",
            anchor_source_count=0, n_mirror_pairs=0,
            n_low_variants_with_mirror=0, mirror_low_fraction=0.,
            normalized_mirror_support=0., mirror_p95_threshold=p95,
            mirror_p99_threshold=p99, mirror_calibration_status=cal_status,
            mirror_support_candidate=False, mirror_support_highconf=False,
            contamination_flag_candidate=False, contamination_flag_highconf=False,
        )
        if not eligible:
            status = "not_tested_" + ineligible_reason
            base.update(contamination_status=status, qc_status="PASS", qc_reason=status)
            out.append(base)
            continue

        low = low_a_rows(own, p)
        lowkeys = {key(r) for r in low}
        otherhigh = {
            s: {key(r) for r in by[species, s] if num(r, "VAF") >= source_high_min}
            for s in others
        }
        best_source, best_overlap, frac, best_overlap_keys = best_source_match(lowkeys, otherhigh)
        dispersion = overlap_dispersion(best_overlap_keys, p)

        sensitivity = {}
        if cluster_sensitivity_available:
            sample_clustered = clustered_keys.get(sample, set())
            low_clustered_keys = lowkeys & sample_clustered
            lowkeys_noncluster = lowkeys - sample_clustered
            source_nc, overlap_nc, frac_nc, overlap_keys_nc = best_source_match(lowkeys_noncluster, otherhigh)
            dispersion_nc = overlap_dispersion(overlap_keys_nc, p, suffix="_noncluster")
            n_overlap_from_cluster = len(best_overlap_keys & sample_clustered)
            retention = overlap_nc / best_overlap if best_overlap else None
            stable = "YES" if best_source and source_nc == best_source else "NO" if best_source or source_nc else "NA"
            sensitivity.update(
                n_lowA_clustered=len(low_clustered_keys),
                n_lowA_noncluster=len(lowkeys_noncluster),
                best_source_sample_noncluster=source_nc,
                best_overlap_noncluster=overlap_nc,
                best_frac_lowA_in_highB_noncluster=frac_nc,
                n_overlap_from_cluster=n_overlap_from_cluster,
                overlap_retention_after_cluster_removal=retention,
                source_stable_after_cluster_removal=stable,
                **dispersion_nc,
            )

        anchor_pool = set().union(*otherhigh.values()) if otherhigh else set()
        ownmap = defaultdict(list)
        for r in own:
            ownmap[key(r)].append(num(r, "VAF"))
        tested = [v for k in anchor_pool for v in ownmap.get(k, [])]
        dep = [v for v in tested if p["mt_lower"] <= v <= p["mt_depressed_upper"]]
        fall = [v for v in tested if p["mt_lower"] <= v <= p["mt_anchor_upper"]]
        if len(dep) >= 3:
            mode, est = "depressed_anchors", 1 - sum(dep) / len(dep)
        elif fall:
            mode, est = "fallback_anchors", 1 - sum(fall) / len(fall)
        else:
            mode, est = "no_anchor_observed", None

        source_count = sum(bool(lowkeys & ks) for ks in otherhigh.values())
        if not others:
            level = "singleton_species"
        elif not anchor_pool:
            level = "no_loo_anchor_pool"
        elif not tested:
            level = "no_anchor_observed_in_sample_A"
        elif len(dep) >= 3:
            level = "strong_multi_source_anchor_support" if source_count > 1 else "strong_single_source_anchor_support"
        elif len(tested) >= 2:
            level = "moderate_multi_source_anchor_support" if source_count > 1 else "moderate_single_source_anchor_support"
        else:
            level = "weak_anchor_support"

        npair, nmir, mfrac, norm = mirror_stats(own, p)
        overlap_candidate = (
            len(lowkeys) >= p["min_n_lowA"]
            and best_overlap >= p["min_overlap"]
            and frac is not None
            and frac >= p["min_frac_lowA_in_highB_candidate"]
        )
        overlap_high = overlap_candidate and frac >= p["min_frac_lowA_in_highB_highconf"]
        mtcand = est is not None and est >= p["contam_threshold_candidate"]
        mthigh = est is not None and est >= p["contam_threshold_highconf"]
        candidate = overlap_candidate and mtcand
        highconf = overlap_high and mthigh
        mcand = p95 is not None and norm >= p95
        mhigh = p99 is not None and norm >= p99

        if not others:
            status = "insufficient_singleton_species"
        elif not own:
            status = "insufficient_variant_data"
        elif highconf:
            status = "high_confidence_contaminated"
        elif candidate:
            status = "candidate_contaminated"
        elif overlap_candidate and not mtcand:
            status = "lowA_highB_overlap_only"
        elif mtcand and not overlap_candidate:
            status = "mt_high_hets_only"
        elif est is None:
            status = "insufficient_anchor_data"
        else:
            status = "no_strong_evidence"
        qc = "FAIL" if status == "high_confidence_contaminated" else "PASS" if status == "no_strong_evidence" else "WARN"

        base.update(
            n_lowA=len(lowkeys), best_source_sample=best_source,
            best_source_qc_status=qc_status.get(best_source, "MISSING") if best_source else "",
            best_overlap=best_overlap, best_frac_lowA_in_highB=frac,
            **dispersion,
            **sensitivity,
            n_anchor_pool_excluding_A=len(anchor_pool),
            n_anchor_tested_in_A=len(tested), n_depressed_anchor=len(dep),
            mt_high_hets_contamination=est, mt_high_hets_mode=mode,
            anchor_evidence_level=level, anchor_source_count=source_count,
            n_mirror_pairs=npair, n_low_variants_with_mirror=nmir,
            mirror_low_fraction=mfrac, normalized_mirror_support=norm,
            mirror_support_candidate=mcand, mirror_support_highconf=mhigh,
            contamination_status=status, contamination_flag_candidate=candidate,
            contamination_flag_highconf=highconf, qc_status=qc, qc_reason=status,
        )
        out.append(base)
    return out, eligibility


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path)
    ap.add_argument("--variant-table", type=Path)
    ap.add_argument("--outdir", type=Path)
    ap.add_argument("--negative-control-pairs", type=Path)
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()
    sec = {}
    if a.config:
        cfg = read_simple_yaml(a.config)
        if "intraspecies_contamination" not in cfg:
            raise ValueError("missing 'intraspecies_contamination' section in configuration")
        sec = cfg.get("intraspecies_contamination") or {}
        if not isinstance(sec, dict):
            raise ValueError("'intraspecies_contamination' must be a YAML mapping")
        if not truth(sec.get("enabled"), False):
            print("[intraspecies] enabled=false")
            return 0

    out = path(a.outdir or sec.get("outdir", "results/qc/intraspecies_contamination"))
    report = out / "reports/intraspecies_contamination_report.tsv"
    if report.exists():
        print(f"[intraspecies] replacing existing report: {report}", file=sys.stderr)
    for d in (out / "logs", out / "reports"):
        d.mkdir(parents=True, exist_ok=True)

    p = {
        **DEFAULTS,
        **{k: v for k, v in sec.items() if k in DEFAULTS},
        "use_snv_only": truth(sec.get("use_snv_only"), True),
        "pass_only": truth(sec.get("pass_only"), True),
    }

    source_list = path(sec.get(
        "source_sample_ref_file",
        "results/qc/sample_deduplication/reports/deduplicated_sample_ref_file.tsv",
    ))
    source_pairs = read_sample_pairs(source_list)
    allowed_samples = {s for _, s in source_pairs}

    sample_qc_path = path(sec.get(
        "sample_qc_report",
        "results/qc/sample_variant_filtering/reports/sample_qc.tsv",
    ))
    qc_status = read_sample_qc(sample_qc_path)

    source_qc_path = path(sec.get(
        "source_variant_qc_report",
        "results/qc/pre_liftover_variant_qc/reports/source_variant_qc.tsv",
    ))
    het_counts = count_strict_hets(source_qc_path, allowed_samples, p)

    cluster_report = path(sec.get(
        "local_cluster_variant_report",
        "results/qc/local_heteroplasmy_qc/reports/local_heteroplasmy_variant_detail.tsv",
    ))
    clustered_keys, cluster_sensitivity_available = load_clustered_variant_keys(cluster_report)

    table = a.variant_table or (path(sec["variant_table"]) if sec.get("variant_table") else None)
    if table is None and truth(sec.get("build_variant_table"), True):
        vcf = path(sec.get("vcf_dir", "results/qc/collected_variant_calling_results/collected_vcf"))
        meta = path(sec.get("sample_summary", "results/qc/collected_variant_calling_results/reports/variant_calling_collection_summary.tsv"))
        table = out / ".work/all_PASS_variants_core_table.tsv"
        table.parent.mkdir(parents=True, exist_ok=True)
        cmd = [
            sys.executable, str(ROOT / "qc_analysis/scripts/build_intraspecies_variant_table.py"),
            "--vcf-dir", str(vcf), "--metadata", str(meta), "--output", str(table),
            "--min-dp", str(sec.get("dp_min", 100)), "--pass-only", "--overwrite",
            "--log-file", str(out / "logs/variant_table_build.log"),
        ]
        if truth(sec.get("use_snv_only"), True):
            cmd.append("--snv-only")
        subprocess.run(cmd, check=True)
    if table is None:
        raise ValueError("build_variant_table=false requires variant_table")

    nc = a.negative_control_pairs or sec.get("negative_control_pairs")
    findings, eligibility = analyse(
        load_rows(path(table)), p, source_pairs, qc_status, het_counts,
        clustered_keys=clustered_keys,
        cluster_sensitivity_available=cluster_sensitivity_available,
        negative_control_pairs=nc,
    )
    write_rows(report, findings, REPORT_COLUMNS)
    write_rows(out / "reports/target_sample_eligibility.tsv", eligibility, ELIGIBILITY_COLUMNS)

    with (out / "run_parameters.tsv").open("w", newline="", encoding="utf-8") as h:
        w = csv.writer(h, delimiter="\t")
        w.writerow(("parameter", "value"))
        w.writerows(sorted(p.items()))
        w.writerow(("source_sample_ref_file", source_list))
        w.writerow(("sample_qc_report", sample_qc_path))
        w.writerow(("source_variant_qc_report", source_qc_path))
        w.writerow(("local_cluster_variant_report", cluster_report))
        w.writerow(("local_cluster_sensitivity_available", cluster_sensitivity_available))
        w.writerow(("variant_table", table))
        w.writerow(("timestamp", dt.datetime.now(dt.timezone.utc).isoformat()))

    n_eligible = sum(str(x.get("target_eligible", "")).lower() == "true" for x in eligibility)
    (out / "logs/intraspecies_contamination.log").write_text(
        f"deduplicated_samples={len(source_pairs)}\n"
        f"eligible_targets={n_eligible}\n"
        f"source_high_vaf_min={p['source_high_vaf_min']}\n"
        f"lowA={p['low_vaf_min']}-{p['low_vaf_max']}\n"
        f"min_alt_reads={p['min_alt_reads']}\n"
        f"dispersion_bin_bp={p['dispersion_bin_bp']}\n"
        f"dispersion_window_bp={p['dispersion_window_bp']}\n"
        f"cluster_sensitivity_available={cluster_sensitivity_available}\n"
        f"cluster_report={cluster_report}\n"
        f"report={report}\n",
        encoding="utf-8",
    )
    print(
        f"[intraspecies] report={report} deduplicated_samples={len(source_pairs)} "
        f"eligible_targets={n_eligible} source_AF>={p['source_high_vaf_min']} "
        f"lowA={p['low_vaf_min']}-{p['low_vaf_max']} AD_alt>={p['min_alt_reads']} "
        f"cluster_sensitivity={'yes' if cluster_sensitivity_available else 'no'}"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, OSError, subprocess.CalledProcessError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        raise SystemExit(2)
