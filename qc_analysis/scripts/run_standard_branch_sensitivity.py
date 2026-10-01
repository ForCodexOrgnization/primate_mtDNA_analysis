#!/usr/bin/env python3
"""Run Standard-branch QC threshold sensitivity at fixed nuclear coverage >=5x.

Experiment:
  nuclear median coverage >= 5x (fixed)
  MAD disabled
  mtCN minimum:               20, 40, 60
  mtDNA median coverage min:  75, 100, 150
  Percent_100 minimum:        85, 90, 95

This creates 27 isolated scenarios. Each scenario reruns the QC steps whose
outputs depend on sample membership or the scenario-specific artifact list:
sample_variant_filtering -> local_heteroplasmy_qc ->
intraspecies_contamination -> interspecies_contamination -> final_filter.

Upstream calling/liftover/codon/tRNA/rRNA outputs are reused read-only.

After all scenarios finish, the script writes scenario_summary.tsv and two
standard plots:
  1) samples/species/unique-variant retention
  2) sample retention vs fixed species-specific P95 heteroplasmy outliers

For the two plots, one threshold is swept at a time while the other two are
fixed at the historical baseline (mtCN=40, mtDNA coverage=100, Percent_100=90).
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
MTCN_VALUES = (20, 40, 60)
MT_COV_VALUES = (75, 100, 150)
PERCENT100_VALUES = (85, 90, 95)
BASELINE = {"mtcn": 40, "mtcov": 100, "pct100": 90}

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
    for mtcn in MTCN_VALUES:
        for mtcov in MT_COV_VALUES:
            for pct100 in PERCENT100_VALUES:
                name = f"mtcn{mtcn}_mtcov{mtcov}_p{pct100}"
                out.append({
                    "scenario": name,
                    "nuclear": NUCLEAR_MIN,
                    "mtcn": mtcn,
                    "mtcov": mtcov,
                    "pct100": pct100,
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


def scenario_config(base_text: str, params: dict, output_root: Path) -> str:
    name = params["scenario"]
    scenario_root = output_root / name
    sample_dir = scenario_root / "sample_variant_filtering"
    local_dir = scenario_root / "local_heteroplasmy_qc"
    intra_dir = scenario_root / "intraspecies_contamination"
    inter_dir = scenario_root / "interspecies_contamination"
    final_dir = scenario_root / "final_filter"

    text = base_text
    text = add_direct_child(text, "sample_variant_filtering", "mad_enabled", False)

    for path, value in (
        (("sample_variant_filtering", "output_dir"), sample_dir),
        (("sample_variant_filtering", "thresholds", "mt_median_coverage_min"), params["mtcov"]),
        (("sample_variant_filtering", "thresholds", "percent_100_min"), params["pct100"]),
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
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
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


def prepare(base_config: Path, output_root: Path):
    base_text = base_config.read_text(encoding="utf-8")
    config_dir = output_root / "configs"
    config_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for index, params in enumerate(scenarios()):
        cfg = config_dir / f"{params['scenario']}.yaml"
        cfg.write_text(scenario_config(base_text, params, output_root), encoding="utf-8")
        rows.append({
            "index": index,
            **params,
            "config": str(cfg),
        })
    write_tsv(output_root / "scenario_manifest.tsv", rows)
    return rows


def run_step(step: str, config: Path):
    cmd = ["bash", str(WRAPPER), step, str(config)]
    print("[standard_sensitivity] RUN", " ".join(cmd), flush=True)
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
            f"{scenario_root} already exists; remove it or rerun submission with --force"
        )
    cfg = Path(row["config"])
    print(
        f"[standard_sensitivity] scenario={row['scenario']} "
        f"nuclear>={row['nuclear']} mtCN>={row['mtcn']} "
        f"mtcov>={row['mtcov']} Percent_100>={row['pct100']} MAD=OFF",
        flush=True,
    )
    for step in STEPS:
        run_step(step, cfg)


def load_source_summary(base_config: Path):
    # The production sample filter currently defaults to this report. Keeping the
    # lookup explicit makes the heteroplasmy reference fixed across all scenarios.
    # A future config parser can replace this without changing the experiment.
    path = ROOT / "results/qc/collected_variant_calling_results/reports/variant_calling_collection_summary.tsv"
    if not path.is_file():
        raise FileNotFoundError(path)
    return read_tsv(path)


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

    pass_species = {str(r.get("species", "")).strip() for r in final_pass_samples if str(r.get("species", "")).strip()}
    artifact_remove = [r for r in local_rows if str(r.get("filter_action", "")).upper() == "REMOVE"]
    artifact_samples = {str(r.get("sample", "")) for r in artifact_remove if str(r.get("sample", ""))}

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
        "percent_100_min": params["pct100"],
        "mad_enabled": "false",
        "sample_qc_pass": len(sample_qc_pass),
        "sample_qc_pass_species": len({str(r.get("species", "")).strip() for r in sample_qc_pass if str(r.get("species", "")).strip()}),
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
        "heteroplasmy_outlier_pct": (100.0 * outlier_n / outlier_assessable) if outlier_assessable else None,
    }


def plot_results(summary_rows, output_root: Path):
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError(
            "matplotlib is required for automatic sensitivity plots. "
            "Install/load matplotlib and rerun --summarize."
        ) from exc

    plots_dir = output_root / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    def f(row, key):
        x = row.get(key)
        return None if x in (None, "") else float(x)

    panels = [
        {
            "label": "mtDNA copy number",
            "xkey": "mtcn_min",
            "filter": lambda r: int(float(r["mt_median_coverage_min"])) == BASELINE["mtcov"]
                               and int(float(r["percent_100_min"])) == BASELINE["pct100"],
        },
        {
            "label": "mtDNA median coverage",
            "xkey": "mt_median_coverage_min",
            "filter": lambda r: int(float(r["mtcn_min"])) == BASELINE["mtcn"]
                               and int(float(r["percent_100_min"])) == BASELINE["pct100"],
        },
        {
            "label": "Percent bases >=100x",
            "xkey": "percent_100_min",
            "filter": lambda r: int(float(r["mtcn_min"])) == BASELINE["mtcn"]
                               and int(float(r["mt_median_coverage_min"])) == BASELINE["mtcov"],
        },
    ]

    retention_data = []
    outlier_data = []
    for panel in panels:
        rows = sorted([r for r in summary_rows if panel["filter"](r)], key=lambda r: f(r, panel["xkey"]))
        if not rows:
            continue
        base_samples = f(rows[0], "final_pass_samples")
        base_species = f(rows[0], "final_pass_species")
        base_unique = f(rows[0], "final_unique_variants")
        for r in rows:
            retention_data.append({
                "metric": panel["label"],
                "threshold": f(r, panel["xkey"]),
                "samples_retained_pct": 100 * f(r, "final_pass_samples") / base_samples if base_samples else None,
                "species_retained_pct": 100 * f(r, "final_pass_species") / base_species if base_species else None,
                "unique_variants_retained_pct": 100 * f(r, "final_unique_variants") / base_unique if base_unique else None,
            })
            outlier_data.append({
                "metric": panel["label"],
                "threshold": f(r, panel["xkey"]),
                "samples_retained_pct": 100 * f(r, "final_pass_samples") / base_samples if base_samples else None,
                "heteroplasmy_outlier_pct": f(r, "heteroplasmy_outlier_pct"),
            })

    write_tsv(plots_dir / "retention_plot_data.tsv", retention_data)
    write_tsv(plots_dir / "heteroplasmy_outlier_plot_data.tsv", outlier_data)

    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    axes = axes.ravel()
    for ax, panel in zip(axes, panels):
        rows = [r for r in retention_data if r["metric"] == panel["label"]]
        x = [r["threshold"] for r in rows]
        ax.plot(x, [r["samples_retained_pct"] for r in rows], marker="o", linestyle="-", label="Samples retained")
        ax.plot(x, [r["species_retained_pct"] for r in rows], marker="^", linestyle="--", label="Species retained")
        ax.plot(x, [r["unique_variants_retained_pct"] for r in rows], marker="s", linestyle="--", label="Unique variants retained")
        ax.set_title(panel["label"])
        ax.set_xlabel("QC threshold")
        ax.set_ylabel("Retention (%)")
        ax.set_ylim(0, 105)
        ax.grid(True, alpha=0.25)
    axes[-1].axis("off")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False)
    fig.suptitle("Standard branch: retention across QC thresholds (nuclear coverage >=5x)")
    fig.tight_layout(rect=(0, 0.07, 1, 0.95))
    fig.savefig(plots_dir / "standard_branch_retention.png", dpi=200)
    fig.savefig(plots_dir / "standard_branch_retention.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    axes = axes.ravel()
    for ax, panel in zip(axes, panels):
        rows = [r for r in outlier_data if r["metric"] == panel["label"]]
        x = [r["threshold"] for r in rows]
        sample_ret = [r["samples_retained_pct"] for r in rows]
        out_pct = [r["heteroplasmy_outlier_pct"] for r in rows]

        ax.plot(x, sample_ret, marker="o", linestyle="-", label="Samples retained")
        ax.set_title(panel["label"])
        ax.set_xlabel("QC threshold")
        ax.set_ylabel("Samples retained (%)")
        ax.set_ylim(0, 105)
        ax.grid(True, alpha=0.25)

        ax2 = ax.twinx()
        ax2.plot(x, out_pct, marker="^", linestyle="--", label="Heteroplasmy outliers")
        ax2.set_ylabel("Species-specific heteroplasmy outliers (%)")
        finite = [v for v in out_pct if v is not None and math.isfinite(v)]
        if finite:
            ax2.set_ylim(0, max(finite) * 1.25 if max(finite) > 0 else 1)
    axes[-1].axis("off")
    # Custom combined legend from the first panel.
    h1, l1 = axes[0].get_legend_handles_labels()
    h2, l2 = axes[0].twinx().get_legend_handles_labels()
    # Recreate handles explicitly because the temporary twin axis is empty.
    from matplotlib.lines import Line2D
    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="-", label="Samples retained"),
        Line2D([0], [0], marker="^", linestyle="--", label="Heteroplasmy outliers"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=2, frameon=False)
    fig.suptitle("Standard branch: sample retention and heteroplasmy outliers (fixed species P95)")
    fig.tight_layout(rect=(0, 0.07, 1, 0.95))
    fig.savefig(plots_dir / "standard_branch_heteroplasmy_outliers.png", dpi=200)
    fig.savefig(plots_dir / "standard_branch_heteroplasmy_outliers.pdf")
    plt.close(fig)


def summarize(output_root: Path, base_config: Path):
    manifest = read_tsv(output_root / "scenario_manifest.tsv")
    source_rows = load_source_summary(base_config)
    source_by_sample = {str(r.get("sample", "")): r for r in source_rows if str(r.get("sample", ""))}
    hetero_ref = build_fixed_heteroplasmy_reference(source_rows)

    ref_rows = [
        {
            "species": species,
            "n_reference_samples": data["n_reference"],
            "p95_heteroplasmy": data["p95_heteroplasmy"],
        }
        for species, data in sorted(hetero_ref.items())
    ]
    write_tsv(output_root / "fixed_species_heteroplasmy_p95.tsv", ref_rows)

    summary_rows = []
    for row in manifest:
        params = {
            "scenario": row["scenario"],
            "nuclear": int(float(row["nuclear"])),
            "mtcn": int(float(row["mtcn"])),
            "mtcov": int(float(row["mtcov"])),
            "pct100": int(float(row["pct100"])),
        }
        scenario_root = output_root / params["scenario"]
        summary_rows.append(summarize_scenario(params, scenario_root, source_by_sample, hetero_ref))

    write_tsv(output_root / "scenario_summary.tsv", summary_rows)
    plot_results(summary_rows, output_root)
    print(f"[standard_sensitivity] summary={output_root / 'scenario_summary.tsv'}")
    print(f"[standard_sensitivity] plots={output_root / 'plots'}")


def submit(base_config: Path, output_root: Path, force: bool, concurrency: int):
    if shutil.which("sbatch") is None:
        raise RuntimeError("sbatch not found; run on a Slurm login node")

    if output_root.exists() and force:
        shutil.rmtree(output_root)
    elif output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"{output_root} already exists and is not empty; use --force to replace it")

    output_root.mkdir(parents=True, exist_ok=True)
    manifest = prepare(base_config, output_root)
    logs = output_root / "logs"
    logs.mkdir(exist_ok=True)

    shell = ROOT / "qc_analysis/scripts/run_standard_branch_sensitivity.sh"
    time = os.environ.get("STANDARD_SENS_TIME", "24:00:00")
    mem = os.environ.get("STANDARD_SENS_MEM", "24G")
    cpus = os.environ.get("STANDARD_SENS_CPUS", "4")

    array_cmd = [
        "sbatch", "--parsable",
        "--job-name=qc_standard_sens",
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
    print(f"[standard_sensitivity] submitted 27-scenario array job={array_job}")

    summary_cmd = [
        "sbatch", "--parsable",
        "--job-name=qc_standard_summary",
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
    print(f"[standard_sensitivity] submitted summary/plot job={summary_job} afterok:{array_job}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base-config", type=Path, default=Path("config/qc_preprocessing.yaml"))
    ap.add_argument("--output-root", type=Path, default=Path("results/qc_threshold_sensitivity/standard_branch"))
    ap.add_argument("--submit", action="store_true")
    ap.add_argument("--array-task", action="store_true")
    ap.add_argument("--summarize", action="store_true")
    ap.add_argument("--run-index", type=int)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--concurrency", type=int, default=4)
    args = ap.parse_args()

    base_config = args.base_config if args.base_config.is_absolute() else (ROOT / args.base_config).resolve()
    output_root = args.output_root if args.output_root.is_absolute() else (ROOT / args.output_root).resolve()

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
        summarize(output_root, base_config)
        return 0

    if output_root.exists() and args.force:
        shutil.rmtree(output_root)
    elif output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"{output_root} already exists and is not empty; use --force")

    output_root.mkdir(parents=True, exist_ok=True)
    prepare(base_config, output_root)
    for i in range(len(scenarios())):
        run_index(i, output_root)
    summarize(output_root, base_config)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, KeyError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
