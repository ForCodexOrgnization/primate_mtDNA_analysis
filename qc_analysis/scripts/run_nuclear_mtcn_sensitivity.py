#!/usr/bin/env python3
"""Run A/B/C/D sensitivity analysis for universal nuclear-coverage and mtCN QC.

The experiment fixes mtDNA median coverage >=100 and Percent_100 >=90, disables
MAD QC, and varies only nuclear coverage and mtCN:
  A_current                 nuclear >=20, mtCN >=40
  B_no_nuclear              nuclear OFF,  mtCN >=40
  C_no_mtcn                 nuclear >=20, mtCN OFF
  D_no_nuclear_no_mtcn      nuclear OFF,  mtCN OFF

Each scenario gets an isolated config/output tree. Upstream variant calling,
pre-liftover QC, liftover, codon/tRNA/rRNA annotations are reused. Steps whose
results depend on sample-QC membership or the scenario-specific artifact list are
rerun: sample_variant_filtering, local_heteroplasmy_qc,
intraspecies_contamination, interspecies_contamination, and final_filter.
"""
from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

SCENARIOS = {
    "A_current": {"nuclear": 20, "mtcn": 40},
    "B_no_nuclear": {"nuclear": 0, "mtcn": 40},
    "C_no_mtcn": {"nuclear": 20, "mtcn": 0},
    "D_no_nuclear_no_mtcn": {"nuclear": 0, "mtcn": 0},
}

STEPS = [
    "sample_variant_filtering",
    "local_heteroplasmy_qc",
    "intraspecies_contamination",
    "interspecies_contamination",
    "final_filter",
]

YAML_KEY = re.compile(r"^(\s*)([^#\s][^:]*):(?:\s*(.*))?$")


def yaml_scalar(value) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if value is None:
        return "null"
    return str(value)


def line_paths(lines: list[str]):
    stack: list[tuple[int, str]] = []
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


def set_yaml_path(text: str, path: tuple[str, ...], value) -> str:
    lines = text.splitlines()
    target = tuple(path)
    for i, indent, current, _ in line_paths(lines):
        if current == target:
            lines[i] = " " * indent + f"{path[-1]}: {yaml_scalar(value)}"
            return "\n".join(lines) + "\n"
    raise KeyError(f"YAML path not found: {'.'.join(path)}")


def add_direct_child(text: str, section: str, key: str, value) -> str:
    lines = text.splitlines()
    direct = (section, key)
    for _, _, path, _ in line_paths(lines):
        if path == direct:
            return set_yaml_path(text, direct, value)

    section_idx = None
    section_indent = None
    insert_idx = None
    for i, indent, path, _ in line_paths(lines):
        if path == (section,):
            section_idx = i
            section_indent = indent
            continue
        if section_idx is not None and i > section_idx and indent <= section_indent:
            insert_idx = i
            break
    if section_idx is None:
        raise KeyError(f"YAML section not found: {section}")
    if insert_idx is None:
        insert_idx = len(lines)

    # Insert immediately after the section header for visibility.
    lines.insert(section_idx + 1, " " * (section_indent + 2) + f"{key}: {yaml_scalar(value)}")
    return "\n".join(lines) + "\n"


def scenario_config(base_text: str, name: str, nuclear: int, mtcn: int, out_root: Path) -> str:
    scenario_root = out_root / name
    sample_dir = scenario_root / "sample_variant_filtering"
    local_dir = scenario_root / "local_heteroplasmy_qc"
    intra_dir = scenario_root / "intraspecies_contamination"
    inter_dir = scenario_root / "interspecies_contamination"
    final_dir = scenario_root / "final_filter"

    text = base_text

    # Experiment-controlled sample QC.
    text = add_direct_child(text, "sample_variant_filtering", "mad_enabled", False)
    for path, value in [
        (("sample_variant_filtering", "output_dir"), sample_dir),
        (("sample_variant_filtering", "thresholds", "mt_median_coverage_min"), 100),
        (("sample_variant_filtering", "thresholds", "percent_100_min"), 90),
        (("sample_variant_filtering", "thresholds", "nuclear_median_coverage_min"), nuclear),
        (("sample_variant_filtering", "thresholds", "mtcn_min"), mtcn),
    ]:
        text = set_yaml_path(text, path, value)

    # Scenario-specific local artifact analysis.
    text = set_yaml_path(text, ("local_heteroplasmy_qc", "output_dir"), local_dir)
    text = add_direct_child(
        text, "local_heteroplasmy_qc", "sample_qc_report",
        sample_dir / "reports/sample_qc.tsv",
    )

    # Scenario-specific intra-species contamination.
    text = set_yaml_path(text, ("intraspecies_contamination", "outdir"), intra_dir)
    text = set_yaml_path(
        text, ("intraspecies_contamination", "artifact_removal_report"),
        local_dir / "reports/numt_variants_to_remove.tsv",
    )
    text = add_direct_child(
        text, "intraspecies_contamination", "sample_qc_report",
        sample_dir / "reports/sample_qc.tsv",
    )

    # Inter-species analysis must also use the scenario-specific artifact list.
    text = set_yaml_path(
        text, ("interspecies_contamination", "paths", "artifact_removal_report"),
        local_dir / "reports/numt_variants_to_remove.tsv",
    )
    text = set_yaml_path(
        text, ("interspecies_contamination", "paths", "output_dir"),
        inter_dir,
    )

    # Final filter reads scenario-specific sample/intra/inter reports.
    text = set_yaml_path(text, ("final_filter", "output_dir"), final_dir)
    text = set_yaml_path(
        text, ("final_filter", "sample_reports", "intraspecies", "path"),
        intra_dir / "reports/intraspecies_contamination_report.tsv",
    )
    text = set_yaml_path(
        text, ("final_filter", "sample_reports", "interspecies", "path"),
        inter_dir / "reports/interspecies_contamination_report.tsv",
    )
    text = set_yaml_path(
        text, ("final_filter", "sample_reports", "sample_qc", "path"),
        sample_dir / "reports/sample_qc.tsv",
    )
    # run_final_filter_with_heteroplasmy.py infers these from local_heteroplasmy_qc.output_dir.

    return text


def read_tsv(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as h:
        return list(csv.DictReader(h, delimiter="\t"))


def truthy(value) -> bool:
    return str(value).strip().upper() in {"TRUE", "YES", "1", "PASS"}


def summarize(name: str, params: dict, scenario_root: Path) -> dict:
    sample_rows = read_tsv(scenario_root / "sample_variant_filtering/reports/sample_qc.tsv")
    local_rows = read_tsv(scenario_root / "local_heteroplasmy_qc/reports/numt_variants_to_remove.tsv")
    intra_rows = read_tsv(scenario_root / "intraspecies_contamination/reports/intraspecies_contamination_report.tsv")
    inter_rows = read_tsv(scenario_root / "interspecies_contamination/reports/interspecies_contamination_report.tsv")
    final_samples = read_tsv(scenario_root / "final_filter/reports/final_sample_qc.tsv")
    final_vars = read_tsv(scenario_root / "final_filter/reports/final_variant_qc.tsv")

    pass_samples = [r for r in sample_rows if str(r.get("qc_status", "")).upper() == "PASS"]
    final_pass_samples = [r for r in final_samples if str(r.get("final_sample_status", "")).upper() == "PASS"]
    final_pass_vars = [r for r in final_vars if str(r.get("final_variant_status", "")).upper() == "PASS"]

    def variant_key(r):
        if r.get("human_pos"):
            return (r.get("human_chrom", ""), r.get("human_pos", ""), r.get("human_ref", ""), r.get("human_alt", ""))
        return (r.get("source_chrom", ""), r.get("source_pos", ""), r.get("source_ref", ""), r.get("source_alt", ""))

    return {
        "scenario": name,
        "nuclear_median_coverage_min": params["nuclear"],
        "mtcn_min": params["mtcn"],
        "mt_median_coverage_min": 100,
        "percent_100_min": 90,
        "mad_enabled": "false",
        "sample_qc_pass": len(pass_samples),
        "sample_qc_pass_species": len({r.get("species", "") for r in pass_samples if r.get("species", "")}),
        "artifact_remove_rows": sum(str(r.get("filter_action", "")).upper() == "REMOVE" for r in local_rows),
        "intra_high_conf_contaminated": sum(str(r.get("contamination_status", "")).lower() == "high_confidence_contaminated" for r in intra_rows),
        "inter_fail": sum(str(r.get("qc_status", "")).upper() == "FAIL" for r in inter_rows),
        "final_pass_samples": len(final_pass_samples),
        "final_pass_species": len({r.get("species", "") for r in final_pass_samples if r.get("species", "")}),
        "final_pass_variant_rows": len(final_pass_vars),
        "final_unique_variants": len({variant_key(r) for r in final_pass_vars}),
    }


def write_summary(rows: list[dict], path: Path) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def run_step(step: str, config: Path) -> None:
    cmd = [
        "bash",
        str(ROOT / "qc_analysis/scripts/run_qc_preprocessing.sh"),
        step,
        str(config),
    ]
    print("[sensitivity] RUN", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=ROOT, check=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base-config", type=Path, default=Path("config/qc_preprocessing.yaml"))
    ap.add_argument("--output-root", type=Path, default=Path("results/qc_threshold_sensitivity/nuclear_mtcn"))
    ap.add_argument("--force", action="store_true", help="Delete and rerun existing scenario output directories.")
    ap.add_argument("--scenarios", default="A,B,C,D", help="Comma-separated subset of A,B,C,D.")
    args = ap.parse_args()

    base = args.base_config if args.base_config.is_absolute() else ROOT / args.base_config
    out_root = args.output_root if args.output_root.is_absolute() else ROOT / args.output_root
    if not base.is_file():
        raise FileNotFoundError(base)

    aliases = {
        "A": "A_current",
        "B": "B_no_nuclear",
        "C": "C_no_mtcn",
        "D": "D_no_nuclear_no_mtcn",
    }
    requested = []
    for token in (x.strip() for x in args.scenarios.split(",")):
        key = aliases.get(token.upper(), token)
        if key not in SCENARIOS:
            raise ValueError(f"unknown scenario: {token}")
        requested.append(key)

    out_root.mkdir(parents=True, exist_ok=True)
    config_dir = out_root / "configs"
    config_dir.mkdir(exist_ok=True)
    base_text = base.read_text(encoding="utf-8")

    summary_rows = []
    for name in requested:
        params = SCENARIOS[name]
        scenario_root = out_root / name
        if scenario_root.exists():
            if args.force:
                print(f"[sensitivity] removing existing {scenario_root}", flush=True)
                shutil.rmtree(scenario_root)
            else:
                raise FileExistsError(
                    f"{scenario_root} already exists. Use --force to rerun and replace this test output."
                )

        cfg_text = scenario_config(base_text, name, params["nuclear"], params["mtcn"], out_root)
        cfg_path = config_dir / f"{name}.yaml"
        cfg_path.write_text(cfg_text, encoding="utf-8")

        print(
            f"\n[sensitivity] ===== {name}: nuclear>={params['nuclear']} mtCN>={params['mtcn']} "
            "mt>=100 Percent_100>=90 MAD=OFF =====",
            flush=True,
        )
        for step in STEPS:
            run_step(step, cfg_path)

        summary_rows.append(summarize(name, params, scenario_root))
        write_summary(summary_rows, out_root / "scenario_summary.tsv")

    print(f"\n[sensitivity] complete: {out_root}", flush=True)
    print(f"[sensitivity] summary: {out_root / 'scenario_summary.tsv'}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
