#!/usr/bin/env python3
"""Run six mtDNA coverage / Percent_X QC sensitivity scenarios.

Fixed:
  nuclear median coverage >= 5x
  mtCN >= 40
  MAD disabled
  all downstream biological/QC rules unchanged

Scenarios:
  mtDNA median coverage >= 60x with Percent_60 >= 85% or 90%
  mtDNA median coverage >= 80x with Percent_80 >= 85% or 90%
  mtDNA median coverage >=100x with Percent_100 >=85% or 90%

Percent_X is recomputed from the collected per-base merged coverage files as
100 * (# bases with coverage >= X) / (# covered-table rows).

Each scenario reruns:
sample_variant_filtering -> local_heteroplasmy_qc ->
intraspecies_contamination -> interspecies_contamination -> final_filter.

Upstream variant calling, pre-liftover QC, liftover and functional annotation
outputs are reused read-only.
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WRAPPER = ROOT / "qc_analysis/scripts/run_qc_preprocessing.sh"

NUCLEAR_MIN = 5
MTCN_MIN = 40
MT_COV_VALUES = (60, 80, 100)
FRACTION_VALUES = (85, 90)

STEPS = (
    "sample_variant_filtering",
    "local_heteroplasmy_qc",
    "intraspecies_contamination",
    "interspecies_contamination",
    "final_filter",
)

YAML_KEY = re.compile(r"^(\s*)([^#\s][^:]*):(?:\s*(.*))?$")


def scenarios():
    out = []
    for mtcov in MT_COV_VALUES:
        for fraction in FRACTION_VALUES:
            out.append({
                "scenario": f"mtcov{mtcov}_percent{mtcov}_p{fraction}",
                "nuclear": NUCLEAR_MIN,
                "mtcn": MTCN_MIN,
                "mtcov": mtcov,
                "percent_depth": mtcov,
                "percent_column": f"Percent_{mtcov}",
                "fraction": fraction,
            })
    return out


def yaml_scalar(value) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if value is None:
        return "null"
    return str(value)


def line_paths(lines):
    stack = []
    for i, line in enumerate(lines):
        m = YAML_KEY.match(line)
        if not m:
            continue
        indent = len(m.group(1).replace("\t", "  "))
        key = m.group(2).strip()
        value = (m.group(3) or "").strip()
        while stack and stack[-1][0] >= indent:
            stack.pop()
        path = tuple(k for _, k in stack) + (key,)
        yield i, indent, path, value
        if value == "" or value.startswith("#"):
            stack.append((indent, key))


def set_yaml_path(text: str, path, value) -> str:
    lines = text.splitlines()
    target = tuple(path)
    for i, indent, current, _ in line_paths(lines):
        if current == target:
            lines[i] = " " * indent + f"{path[-1]}: {yaml_scalar(value)}"
            return "\n".join(lines) + "\n"
    raise KeyError(f"YAML path not found: {'.'.join(path)}")


def add_direct_child(text: str, section: str, key: str, value) -> str:
    direct = (section, key)
    lines = text.splitlines()
    for _, _, path, _ in line_paths(lines):
        if path == direct:
            return set_yaml_path(text, direct, value)

    section_idx = None
    section_indent = None
    for i, indent, path, _ in line_paths(lines):
        if path == (section,):
            section_idx = i
            section_indent = indent
            break
    if section_idx is None:
        raise KeyError(f"YAML section not found: {section}")

    lines.insert(
        section_idx + 1,
        " " * (section_indent + 2) + f"{key}: {yaml_scalar(value)}",
    )
    return "\n".join(lines) + "\n"


def read_tsv(path: Path):
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def as_float(value):
    try:
        x = float(value)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def quantile_type7(values, p):
    vals = sorted(v for v in values if v is not None and math.isfinite(v))
    n = len(vals)
    if n == 0:
        return None
    if n == 1:
        return vals[0]
    h = (n - 1) * p
    lo = int(math.floor(h))
    hi = int(math.ceil(h))
    if lo == hi:
        return vals[lo]
    return vals[lo] + (h - lo) * (vals[hi] - vals[lo])


def variant_key(row):
    human = tuple(str(row.get(k, "")) for k in ("human_chrom", "human_pos", "human_ref", "human_alt"))
    if all(human):
        return ("human",) + human
    source = tuple(str(row.get(k, "")) for k in ("source_chrom", "source_pos", "source_ref", "source_alt"))
    if all(source):
        return ("source",) + source
    original = tuple(str(row.get(k, "")) for k in ("original_chrom", "original_pos", "original_ref", "original_alt"))
    return ("original",) + original


def source_summary_path() -> Path:
    return ROOT / "results/qc/collected_variant_calling_results/reports/variant_calling_collection_summary.tsv"


def resolve_cov_path(row: dict) -> Path:
    raw = str(row.get("cov_file", "")).strip()
    candidates = []
    if raw and raw not in {"NA", "."}:
        p = Path(raw).expanduser()
        candidates.append(p if p.is_absolute() else ROOT / p)
    sample = str(row.get("sample", "")).strip()
    candidates.append(
        ROOT / "results/qc/collected_variant_calling_results/collected_cov"
        / f"{sample}.merged.max_depth.per_base_coverage.tsv"
    )
    for path in candidates:
        if path.is_file():
            return path
    raise FileNotFoundError(
        f"merged per-base coverage not found for sample={sample}; checked: "
        + ", ".join(str(p) for p in candidates)
    )


def percent_at_depth(cov_path: Path, depth: int) -> float:
    n = 0
    passed = 0
    with cov_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 4:
                raise ValueError(f"coverage row has fewer than 4 columns: {cov_path}: {row}")
            if row[0].lower() == "chrom" and row[1].lower() == "pos":
                continue
            cov = float(row[3])
            n += 1
            if cov >= depth:
                passed += 1
    if n == 0:
        raise ValueError(f"coverage table has zero data rows: {cov_path}")
    return 100.0 * passed / n


def build_augmented_summary(output_root: Path) -> Path:
    src = source_summary_path()
    if not src.is_file():
        raise FileNotFoundError(src)
    rows = read_tsv(src)
    out_rows = []
    for i, row in enumerate(rows, 1):
        sample = str(row.get("sample", "")).strip()
        if not sample:
            continue
        cov_path = resolve_cov_path(row)
        new = dict(row)
        for depth in MT_COV_VALUES:
            new[f"Percent_{depth}"] = f"{percent_at_depth(cov_path, depth):.6g}"
        out_rows.append(new)
        if i % 250 == 0:
            print(f"[mtcov_percent_sensitivity] percent metrics {i}/{len(rows)}", flush=True)

    out = output_root / "sensitivity_input_summary.tsv"
    write_tsv(out, out_rows)
    print(f"[mtcov_percent_sensitivity] augmented_summary={out}", flush=True)
    return out


def scenario_config(base_text: str, params: dict, output_root: Path, input_summary: Path) -> str:
    name = params["scenario"]
    scenario_root = output_root / name
    sample_dir = scenario_root / "sample_variant_filtering"
    local_dir = scenario_root / "local_heteroplasmy_qc"
    intra_dir = scenario_root / "intraspecies_contamination"
    inter_dir = scenario_root / "interspecies_contamination"
    final_dir = scenario_root / "final_filter"

    text = base_text
    text = add_direct_child(text, "sample_variant_filtering", "mad_enabled", False)
    text = add_direct_child(
        text,
        "sample_variant_filtering",
        "percent_coverage_column",
        params["percent_column"],
    )

    for path, value in (
        (("sample_variant_filtering", "input_summary"), input_summary),
        (("sample_variant_filtering", "output_dir"), sample_dir),
        (("sample_variant_filtering", "thresholds", "mt_median_coverage_min"), params["mtcov"]),
        (("sample_variant_filtering", "thresholds", "percent_100_min"), params["fraction"]),
        (("sample_variant_filtering", "thresholds", "nuclear_median_coverage_min"), params["nuclear"]),
        (("sample_variant_filtering", "thresholds", "mtcn_min"), params["mtcn"]),
    ):
        text = set_yaml_path(text, path, value)

    text = set_yaml_path(text, ("local_heteroplasmy_qc", "output_dir"), local_dir)
    text = add_direct_child(
        text,
        "local_heteroplasmy_qc",
        "sample_qc_report",
        sample_dir / "reports/sample_qc.tsv",
    )

    text = set_yaml_path(text, ("intraspecies_contamination", "outdir"), intra_dir)
    text = set_yaml_path(
        text,
        ("intraspecies_contamination", "artifact_removal_report"),
        local_dir / "reports/numt_variants_to_remove.tsv",
    )
    text = add_direct_child(
        text,
        "intraspecies_contamination",
        "sample_qc_report",
        sample_dir / "reports/sample_qc.tsv",
    )

    text = set_yaml_path(
        text,
        ("interspecies_contamination", "paths", "artifact_removal_report"),
        local_dir / "reports/numt_variants_to_remove.tsv",
    )
    text = set_yaml_path(
        text,
        ("interspecies_contamination", "paths", "output_dir"),
        inter_dir,
    )

    text = set_yaml_path(text, ("final_filter", "output_dir"), final_dir)
    text = set_yaml_path(
        text,
        ("final_filter", "sample_reports", "intraspecies", "path"),
        intra_dir / "reports/intraspecies_contamination_report.tsv",
    )
    text = set_yaml_path(
        text,
        ("final_filter", "sample_reports", "interspecies", "path"),
        inter_dir / "reports/interspecies_contamination_report.tsv",
    )
    text = set_yaml_path(
        text,
        ("final_filter", "sample_reports", "sample_qc", "path"),
        sample_dir / "reports/sample_qc.tsv",
    )
    return text


def prepare(base_config: Path, output_root: Path):
    base_text = base_config.read_text(encoding="utf-8")
    input_summary = build_augmented_summary(output_root)
    config_dir = output_root / "configs"
    config_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for index, params in enumerate(scenarios()):
        cfg = config_dir / f"{params['scenario']}.yaml"
        cfg.write_text(
            scenario_config(base_text, params, output_root, input_summary),
            encoding="utf-8",
        )
        rows.append({
            "index": index,
            **params,
            "input_summary": str(input_summary),
            "config": str(cfg),
        })
    write_tsv(output_root / "scenario_manifest.tsv", rows)
    return rows


def run_step(step: str, config: Path):
    cmd = ["bash", str(WRAPPER), step, str(config)]
    print("[mtcov_percent_sensitivity] RUN", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=ROOT, check=True)


def run_index(index: int, output_root: Path):
    manifest = read_tsv(output_root / "scenario_manifest.tsv")
    matches = [r for r in manifest if int(r["index"]) == index]
    if len(matches) != 1:
        raise ValueError(f"scenario index {index} not found uniquely in manifest")
    row = matches[0]
    scenario_root = output_root / row["scenario"]
    if scenario_root.exists():
        raise FileExistsError(
            f"{scenario_root} already exists; remove it before rerunning this scenario"
        )
    cfg = Path(row["config"])
    print(
        f"[mtcov_percent_sensitivity] scenario={row['scenario']} "
        f"nuclear>={row['nuclear']} mtCN>={row['mtcn']} "
        f"mtcov>={row['mtcov']} {row['percent_column']}>={row['fraction']}% MAD=OFF",
        flush=True,
    )
    for step in STEPS:
        run_step(step, cfg)


def build_fixed_heteroplasmy_reference(source_rows):
    by_species = {}
    for row in source_rows:
        species = str(row.get("species", "")).strip()
        nuc = as_float(row.get("nuclear_median_coverage"))
        het = as_float(row.get("n_hetero"))
        if not species or nuc is None or het is None or nuc < NUCLEAR_MIN:
            continue
        by_species.setdefault(species, []).append(het)
    ref = {}
    for species, vals in by_species.items():
        if len(vals) >= 5:
            ref[species] = {
                "n_reference": len(vals),
                "p95_heteroplasmy": quantile_type7(vals, 0.95),
            }
    return ref


def summarize_scenario(params, scenario_root: Path, source_by_sample, hetero_ref):
    sample_rows = read_tsv(scenario_root / "sample_variant_filtering/reports/sample_qc.tsv")
    local_rows = read_tsv(scenario_root / "local_heteroplasmy_qc/reports/numt_variants_to_remove.tsv")
    intra_rows = read_tsv(scenario_root / "intraspecies_contamination/reports/intraspecies_contamination_report.tsv")
    inter_rows = read_tsv(scenario_root / "interspecies_contamination/reports/interspecies_contamination_report.tsv")
    final_samples = read_tsv(scenario_root / "final_filter/reports/final_sample_qc.tsv")
    final_vars = read_tsv(scenario_root / "final_filter/reports/final_variant_qc.tsv")

    sample_qc_pass = [r for r in sample_rows if str(r.get("qc_status", "")).upper() == "PASS"]
    final_pass_samples = [r for r in final_samples if str(r.get("final_sample_status", "")).upper() == "PASS"]
    final_pass_ids = {str(r.get("sample", "")) for r in final_pass_samples}
    final_pass_vars = [r for r in final_vars if str(r.get("final_variant_status", "")).upper() == "PASS"]

    pass_species = {
        str(r.get("species", "")).strip()
        for r in final_pass_samples
        if str(r.get("species", "")).strip()
    }
    artifact_remove = [r for r in local_rows if str(r.get("filter_action", "")).upper() == "REMOVE"]
    artifact_samples = {
        str(r.get("sample", ""))
        for r in artifact_remove
        if str(r.get("sample", ""))
    }

    intra_fail = {
        str(r.get("sample", ""))
        for r in intra_rows
        if str(r.get("contamination_status", "")).lower() == "high_confidence_contaminated"
        or str(r.get("qc_status", "")).upper() == "FAIL"
    }
    inter_fail = {
        str(r.get("sample", ""))
        for r in inter_rows
        if str(r.get("qc_status", "")).upper() == "FAIL"
        or str(r.get("classification", "")).upper() == "INTERSPECIES_CONTAMINATION"
    }

    het_values = []
    homo_values = []
    outlier_assessable = 0
    outlier_n = 0
    for sample in final_pass_ids:
        src = source_by_sample.get(sample)
        if not src:
            continue
        het = as_float(src.get("n_hetero"))
        homo = as_float(src.get("n_homo"))
        if het is not None:
            het_values.append(het)
        if homo is not None:
            homo_values.append(homo)
        species = str(src.get("species", "")).strip()
        href = hetero_ref.get(species)
        if het is not None and href is not None:
            outlier_assessable += 1
            if het > href["p95_heteroplasmy"]:
                outlier_n += 1

    return {
        "scenario": params["scenario"],
        "nuclear_median_coverage_min": params["nuclear"],
        "mtcn_min": params["mtcn"],
        "mt_median_coverage_min": params["mtcov"],
        "percent_depth_threshold": params["percent_depth"],
        "percent_coverage_column": params["percent_column"],
        "percent_fraction_min": params["fraction"],
        "mad_enabled": "false",
        "sample_qc_pass": len(sample_qc_pass),
        "sample_qc_pass_species": len({
            str(r.get("species", "")).strip()
            for r in sample_qc_pass
            if str(r.get("species", "")).strip()
        }),
        "final_pass_samples": len(final_pass_ids),
        "final_pass_species": len(pass_species),
        "final_pass_variant_rows": len(final_pass_vars),
        "final_unique_variants": len({variant_key(r) for r in final_pass_vars}),
        "artifact_variants_removed": len(artifact_remove),
        "artifact_samples": len(artifact_samples),
        "intra_fail_samples": len({x for x in intra_fail if x}),
        "inter_fail_samples": len({x for x in inter_fail if x}),
        "median_heteroplasmy": quantile_type7(het_values, 0.50),
        "p95_heteroplasmy": quantile_type7(het_values, 0.95),
        "median_homoplasmy": quantile_type7(homo_values, 0.50),
        "heteroplasmy_outlier_assessable": outlier_assessable,
        "heteroplasmy_outlier_n": outlier_n,
        "heteroplasmy_outlier_pct": (
            100.0 * outlier_n / outlier_assessable if outlier_assessable else None
        ),
    }


def plot_results(summary_rows, output_root: Path):
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError(
            "matplotlib is required for automatic sensitivity plots. "
            "Install/load matplotlib and rerun --summarize."
        ) from exc

    rows = sorted(
        summary_rows,
        key=lambda r: (
            int(r["percent_fraction_min"]),
            int(r["mt_median_coverage_min"]),
        ),
    )
    plots_dir = output_root / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    fractions = sorted({int(r["percent_fraction_min"]) for r in rows})
    fig, axes = plt.subplots(1, len(fractions), figsize=(6 * len(fractions), 5), squeeze=False)
    for ax, frac in zip(axes[0], fractions):
        rr = [r for r in rows if int(r["percent_fraction_min"]) == frac]
        x = [int(r["mt_median_coverage_min"]) for r in rr]
        ax.plot(x, [int(r["final_pass_samples"]) for r in rr], marker="o", label="Final samples")
        ax.plot(x, [int(r["final_pass_species"]) for r in rr], marker="^", label="Final species")
        ax.set_title(f"Percent_X >= {frac}%")
        ax.set_xlabel("Matched mtDNA coverage / Percent_X depth")
        ax.set_ylabel("Count")
        ax.set_xticks(MT_COV_VALUES)
        ax.grid(True, alpha=0.25)
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False)
    fig.suptitle("mtDNA coverage / Percent_X sensitivity")
    fig.tight_layout(rect=(0, 0.08, 1, 0.93))
    fig.savefig(plots_dir / "mtcov_percent_retention.png", dpi=200)
    fig.savefig(plots_dir / "mtcov_percent_retention.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, len(fractions), figsize=(6 * len(fractions), 5), squeeze=False)
    for ax, frac in zip(axes[0], fractions):
        rr = [r for r in rows if int(r["percent_fraction_min"]) == frac]
        x = [int(r["mt_median_coverage_min"]) for r in rr]
        ax.plot(
            x,
            [as_float(r.get("heteroplasmy_outlier_pct")) for r in rr],
            marker="o",
            label="Species-specific P95 outliers",
        )
        ax.set_title(f"Percent_X >= {frac}%")
        ax.set_xlabel("Matched mtDNA coverage / Percent_X depth")
        ax.set_ylabel("Heteroplasmy outliers (%)")
        ax.set_xticks(MT_COV_VALUES)
        ax.grid(True, alpha=0.25)
    fig.suptitle("Fixed species-P95 heteroplasmy outliers")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(plots_dir / "mtcov_percent_heteroplasmy_outliers.png", dpi=200)
    fig.savefig(plots_dir / "mtcov_percent_heteroplasmy_outliers.pdf")
    plt.close(fig)


def summarize(output_root: Path):
    manifest = read_tsv(output_root / "scenario_manifest.tsv")
    source_rows = read_tsv(source_summary_path())
    source_by_sample = {
        str(r.get("sample", "")): r
        for r in source_rows
        if str(r.get("sample", ""))
    }
    hetero_ref = build_fixed_heteroplasmy_reference(source_rows)
    write_tsv(
        output_root / "fixed_species_heteroplasmy_p95.tsv",
        [
            {
                "species": species,
                "n_reference_samples": data["n_reference"],
                "p95_heteroplasmy": data["p95_heteroplasmy"],
            }
            for species, data in sorted(hetero_ref.items())
        ],
    )

    summary_rows = []
    for row in manifest:
        params = {
            "scenario": row["scenario"],
            "nuclear": int(float(row["nuclear"])),
            "mtcn": int(float(row["mtcn"])),
            "mtcov": int(float(row["mtcov"])),
            "percent_depth": int(float(row["percent_depth"])),
            "percent_column": row["percent_column"],
            "fraction": int(float(row["fraction"])),
        }
        summary_rows.append(
            summarize_scenario(
                params,
                output_root / params["scenario"],
                source_by_sample,
                hetero_ref,
            )
        )

    write_tsv(output_root / "scenario_summary.tsv", summary_rows)
    plot_results(summary_rows, output_root)
    print(f"[mtcov_percent_sensitivity] summary={output_root / 'scenario_summary.tsv'}")
    print(f"[mtcov_percent_sensitivity] plots={output_root / 'plots'}")


def submit(base_config: Path, output_root: Path, force: bool, concurrency: int):
    if shutil.which("sbatch") is None:
        raise RuntimeError("sbatch not found; run on a Slurm login node")

    if output_root.exists() and force:
        shutil.rmtree(output_root)
    elif output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(
            f"{output_root} already exists and is not empty; use --force to replace it"
        )

    output_root.mkdir(parents=True, exist_ok=True)
    manifest = prepare(base_config, output_root)
    logs = output_root / "logs"
    logs.mkdir(exist_ok=True)

    shell = ROOT / "qc_analysis/scripts/run_mtcov_percent_sensitivity.sh"
    time = os.environ.get("MTCOV_PERCENT_SENS_TIME", "24:00:00")
    mem = os.environ.get("MTCOV_PERCENT_SENS_MEM", "24G")
    cpus = os.environ.get("MTCOV_PERCENT_SENS_CPUS", "4")

    array_cmd = [
        "sbatch", "--parsable",
        "--job-name=qc_mtcov_pct_sens",
        f"--array=0-{len(manifest)-1}%{concurrency}",
        f"--output={logs}/scenario_%A_%a.out",
        f"--error={logs}/scenario_%A_%a.err",
        f"--time={time}",
        f"--mem={mem}",
        f"--cpus-per-task={cpus}",
        str(shell),
        "--array-task",
        "--base-config", str(base_config),
        "--output-root", str(output_root),
    ]
    result = subprocess.run(array_cmd, cwd=ROOT, check=True, text=True, capture_output=True)
    array_job = result.stdout.strip().split(";")[0]
    print(
        f"[mtcov_percent_sensitivity] submitted {len(manifest)}-scenario array job={array_job}"
    )

    summary_cmd = [
        "sbatch", "--parsable",
        "--job-name=qc_mtcov_pct_summary",
        f"--dependency=afterok:{array_job}",
        f"--output={logs}/summary_%j.out",
        f"--error={logs}/summary_%j.err",
        "--time=01:00:00",
        "--mem=8G",
        "--cpus-per-task=1",
        str(shell),
        "--summarize",
        "--base-config", str(base_config),
        "--output-root", str(output_root),
    ]
    result = subprocess.run(summary_cmd, cwd=ROOT, check=True, text=True, capture_output=True)
    summary_job = result.stdout.strip().split(";")[0]
    print(
        f"[mtcov_percent_sensitivity] submitted summary job={summary_job} afterok:{array_job}"
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base-config", type=Path, default=Path("config/qc_preprocessing.yaml"))
    ap.add_argument(
        "--output-root",
        type=Path,
        default=Path("results/qc_threshold_sensitivity/mtcov_percent"),
    )
    ap.add_argument("--submit", action="store_true")
    ap.add_argument("--array-task", action="store_true")
    ap.add_argument("--summarize", action="store_true")
    ap.add_argument("--run-index", type=int)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--concurrency", type=int, default=3)
    args = ap.parse_args()

    base_config = (
        args.base_config
        if args.base_config.is_absolute()
        else (ROOT / args.base_config).resolve()
    )
    output_root = (
        args.output_root
        if args.output_root.is_absolute()
        else (ROOT / args.output_root).resolve()
    )

    if args.submit:
        submit(base_config, output_root, args.force, args.concurrency)
        return 0

    if args.array_task:
        raw = os.environ.get("SLURM_ARRAY_TASK_ID")
        if raw is None:
            raise RuntimeError("--array-task requires SLURM_ARRAY_TASK_ID")
        run_index(int(raw), output_root)
        return 0

    if args.run_index is not None:
        run_index(args.run_index, output_root)
        return 0

    if args.summarize:
        summarize(output_root)
        return 0

    if output_root.exists() and args.force:
        shutil.rmtree(output_root)
    elif output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"{output_root} already exists and is not empty; use --force")

    output_root.mkdir(parents=True, exist_ok=True)
    prepare(base_config, output_root)
    for i in range(len(scenarios())):
        run_index(i, output_root)
    summarize(output_root)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, KeyError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
