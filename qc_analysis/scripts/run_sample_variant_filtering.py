#!/usr/bin/env python3
"""Evaluate biological sample QC criteria without removing samples.

The historical configuration uses Percent_100 / percent_100_min.  Sensitivity
experiments may select another precomputed Percent_X column via
sample_variant_filtering.percent_coverage_column and threshold it with
thresholds.percent_coverage_min.  The historical fields remain in the report
for backward compatibility.
"""
from __future__ import annotations
import argparse, csv, math, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from qc_analysis.lib.simple_yaml import read_simple_yaml

FIELDS = (
    "sample species mt_median_coverage Percent_100 percent_coverage_metric "
    "percent_coverage_value nuclear_median_coverage mtcn_median MAD "
    "pass_mt_coverage pass_percent_100 pass_percent_coverage "
    "pass_nuclear_coverage pass_mtcn pass_mad n_failed_criteria "
    "failed_criteria qc_status collection_status"
).split()

ALIASES = {
    "mt_median_coverage": ("mt_median_coverage", "mt_median_cov", "mt_coverage_median"),
    "Percent_100": ("Percent_100", "percent_100", "percent_mt_coverage_100"),
    "nuclear_median_coverage": ("nuclear_median_coverage", "nuclear_median_cov"),
    "mtcn_median": ("mtcn_median", "mtCN_median", "median_mtcn"),
    "MAD": ("MAD", "mad"),
}


def resolve(v):
    p = Path(str(v)).expanduser()
    return p if p.is_absolute() else ROOT / p


def value(row, names):
    low = {k.lower(): v for k, v in row.items()}
    for name in names:
        raw = low.get(name.lower(), "")
        try:
            x = float(raw)
            if math.isfinite(x):
                return x
        except (TypeError, ValueError):
            pass
    return None


def evaluate(row, thresholds, mad_enabled=True, percent_column="Percent_100"):
    vals = {k: value(row, v) for k, v in ALIASES.items()}
    selected_percent = value(row, (percent_column,))
    percent_min = float(thresholds.get("percent_coverage_min", thresholds["percent_100_min"]))

    selected_percent_pass = (
        selected_percent is not None and selected_percent >= percent_min
    )

    legacy_percent_pass = (
        vals["Percent_100"] is not None
        and vals["Percent_100"] >= float(thresholds["percent_100_min"])
    )

    checks = {
        "pass_mt_coverage": (
            vals["mt_median_coverage"] is not None
            and vals["mt_median_coverage"] >= float(thresholds["mt_median_coverage_min"])
        ),
        "pass_percent_coverage": selected_percent_pass,
        "pass_nuclear_coverage": (
            vals["nuclear_median_coverage"] is not None
            and vals["nuclear_median_coverage"] >= float(thresholds["nuclear_median_coverage_min"])
        ),
        "pass_mtcn": (
            vals["mtcn_median"] is not None
            and vals["mtcn_median"] >= float(thresholds["mtcn_min"])
        ),
        "pass_mad": (
            True
            if not mad_enabled
            else vals["MAD"] is not None and vals["MAD"] < float(thresholds["mad_max"])
        ),
    }

    percent_reason = (
        "low_percent_100"
        if percent_column.lower() == "percent_100"
        else f"low_{percent_column.lower()}"
    )
    reasons = {
        "pass_mt_coverage": "low_mt_coverage",
        "pass_percent_coverage": percent_reason,
        "pass_nuclear_coverage": "low_nuclear_coverage",
        "pass_mtcn": "low_mtcn",
        "pass_mad": "high_MAD",
    }
    failed = [reasons[k] for k, ok in checks.items() if not ok]

    return {
        "sample": row.get("sample", row.get("Sample", "")),
        "species": row.get("species", row.get("Species", "")),
        **{k: "NA" if v is None else v for k, v in vals.items()},
        "percent_coverage_metric": percent_column,
        "percent_coverage_value": "NA" if selected_percent is None else selected_percent,
        "pass_mt_coverage": checks["pass_mt_coverage"],
        "pass_percent_100": legacy_percent_pass,
        "pass_percent_coverage": selected_percent_pass,
        "pass_nuclear_coverage": checks["pass_nuclear_coverage"],
        "pass_mtcn": checks["pass_mtcn"],
        "pass_mad": checks["pass_mad"],
        "n_failed_criteria": len(failed),
        "failed_criteria": ";".join(failed),
        "qc_status": "PASS" if not failed else "FAIL",
        "collection_status": row.get("status", ""),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()

    sec = read_simple_yaml(a.config).get("sample_variant_filtering") or {}
    if sec.get("enabled", True) is False:
        print("[sample_variant_filtering] disabled; skipping.")
        return 0

    inp = resolve(
        sec.get(
            "input_summary",
            "results/qc/collected_variant_calling_results/reports/variant_calling_collection_summary.tsv",
        )
    )
    out = resolve(sec.get("output_dir", "results/qc/sample_variant_filtering"))
    report = out / "reports/sample_qc.tsv"
    if not inp.is_file():
        raise ValueError(f"missing collection summary: {inp}")
    if report.exists():
        print(
            f"[sample_variant_filtering] replacing existing report: {report}",
            file=sys.stderr,
        )

    t = {
        "mt_median_coverage_min": 100,
        "percent_100_min": 90,
        "nuclear_median_coverage_min": 20,
        "mtcn_min": 40,
        "mad_max": 0.5,
        **(sec.get("thresholds") or {}),
    }
    percent_column = str(sec.get("percent_coverage_column", "Percent_100")).strip() or "Percent_100"
    mad_enabled = sec.get("mad_enabled", True) is not False

    with inp.open(newline="", encoding="utf-8") as h:
        rows = [
            evaluate(
                r,
                t,
                mad_enabled=mad_enabled,
                percent_column=percent_column,
            )
            for r in csv.DictReader(h, delimiter="\t")
        ]

    out.joinpath("reports").mkdir(parents=True, exist_ok=True)
    out.joinpath("logs").mkdir(exist_ok=True)
    with report.open("w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    out.joinpath("logs/sample_variant_filtering.log").write_text(
        f"samples={len(rows)} "
        f"pass={sum(r['qc_status'] == 'PASS' for r in rows)} "
        f"mad_enabled={mad_enabled} "
        f"percent_coverage_column={percent_column} "
        f"percent_coverage_min={t.get('percent_coverage_min', t['percent_100_min'])}\n"
    )
    print(
        f"[sample_variant_filtering] report={report} samples={len(rows)} "
        f"percent_metric={percent_column}"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, KeyError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        raise SystemExit(2)
