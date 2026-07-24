#!/usr/bin/env python3
"""Run intra-species contamination QC from config or direct CLI arguments."""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Any, Mapping

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from qc_analysis.lib.simple_yaml import read_simple_yaml

R_PARAMETERS = {
    "dp_min", "use_snv_only", "low_vaf_min", "low_vaf_max", "high_vaf_min",
    "mt_lower", "mt_depressed_upper", "mt_anchor_upper", "min_n_lowA", "min_overlap",
    "min_frac_lowA_in_highB_candidate", "min_frac_lowA_in_highB_highconf",
    "contam_threshold_candidate", "contam_threshold_highconf", "mirror_low_vaf_min",
    "mirror_low_vaf_max", "mirror_high_vaf_min", "mirror_high_vaf_max", "mirror_tolerance",
    "min_mirror_pairs_for_raw_flag", "min_low_variants_with_mirror_for_flag",
    "target_negative_control_tier",
}


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", type=Path)
    p.add_argument("--variant-table", type=Path)
    p.add_argument("--outdir", type=Path)
    p.add_argument("--negative-control-pairs", type=Path)
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--rscript", default=None)
    return p


def boolean(section: Mapping[str, Any], name: str, default: bool) -> bool:
    value = section.get(name, default)
    if isinstance(value, bool):
        return value
    raise ValueError(f"intraspecies_contamination.{name} must be true or false")


def text_value(section: Mapping[str, Any], name: str, default: Any = None) -> str | None:
    value = section.get(name, default)
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (str, int, float)):
        raise ValueError(f"intraspecies_contamination.{name} must be a scalar or null")
    return str(value)


def resolve_path(value: str | None) -> Path | None:
    if value is None or value == "":
        return None
    path = Path(value).expanduser()
    return path if path.is_absolute() else REPO_ROOT / path


def config_values(config_path: Path) -> dict[str, Any]:
    try:
        data = read_simple_yaml(config_path)
    except (OSError, ValueError) as error:
        raise ValueError(f"failed to parse configuration: {error}") from error
    if not isinstance(data, dict):
        raise ValueError("YAML root must be a mapping")
    section = data.get("intraspecies_contamination")
    if section is None:
        raise ValueError("missing 'intraspecies_contamination' section in configuration")
    if not isinstance(section, dict):
        raise ValueError("'intraspecies_contamination' must be a YAML mapping")
    values: dict[str, Any] = dict(section)
    values.update({
        "enabled": boolean(section, "enabled", False),
        "build_variant_table": boolean(section, "build_variant_table", True),
        "overwrite": boolean(section, "overwrite", False),
        "use_snv_only": boolean(section, "use_snv_only", True),
        "vcf_dir": text_value(section, "vcf_dir"),
        "metadata": text_value(section, "metadata"),
        "variant_table": text_value(section, "variant_table"),
        "negative_control_pairs": text_value(section, "negative_control_pairs"),
        "outdir": text_value(section, "outdir"),
        "dp_min": text_value(section, "dp_min", 100),
    })
    if not values["outdir"]:
        raise ValueError("intraspecies_contamination.outdir must be nonempty")
    try:
        if int(values["dp_min"]) < 0:
            raise ValueError
    except (TypeError, ValueError):
        raise ValueError("intraspecies_contamination.dp_min must be a non-negative integer") from None
    return values


def print_values(config: Path | None, values: Mapping[str, Any]) -> None:
    if config:
        print(f"[intraspecies] config={config}")
    for name in ("enabled", "build_variant_table", "vcf_dir", "metadata", "variant_table", "negative_control_pairs", "outdir", "dp_min", "use_snv_only"):
        value = values.get(name)
        if isinstance(value, bool): value = str(value).lower()
        print(f"[intraspecies] {name}={value if value not in (None, '') else '<not set>'}")


def run(args: argparse.Namespace) -> int:
    if args.config:
        values = config_values(args.config.resolve())
        print_values(args.config.resolve(), values)
        if not values["enabled"]:
            print("[intraspecies] disabled; skipping.")
            return 0
        build = values["build_variant_table"]
        outdir = resolve_path(values["outdir"])
        variant_table = resolve_path(values["variant_table"])
        negative_controls = resolve_path(values["negative_control_pairs"])
        overwrite = values["overwrite"]
    else:
        if not args.variant_table or not args.outdir:
            raise ValueError("direct CLI mode requires --variant-table and --outdir")
        values = {}
        build = False
        outdir, variant_table = args.outdir.resolve(), args.variant_table.resolve()
        negative_controls, overwrite = args.negative_control_pairs, args.overwrite
        print_values(None, {"enabled": True, "build_variant_table": False, "variant_table": variant_table, "outdir": outdir})
    assert outdir is not None
    (outdir / "input").mkdir(parents=True, exist_ok=True)
    (outdir / "logs").mkdir(parents=True, exist_ok=True)
    if build:
        vcf_dir, metadata = resolve_path(values["vcf_dir"]), resolve_path(values["metadata"])
        if not vcf_dir or not metadata:
            raise ValueError("build mode requires vcf_dir and metadata")
        variant_table = outdir / "input" / "all_PASS_variants_core_table.tsv"
        command = [sys.executable, str(REPO_ROOT / "qc_analysis/scripts/build_intraspecies_variant_table.py"), "--vcf-dir", str(vcf_dir), "--metadata", str(metadata), "--output", str(variant_table), "--min-dp", str(values["dp_min"]), "--pass-only", "--log-file", str(outdir / "input/variant_table_build_warnings.log")]
        if values["use_snv_only"]: command.append("--snv-only")
        if overwrite: command.append("--overwrite")
        subprocess.run(command, check=True)
    elif variant_table is None:
        raise ValueError("pre-built mode requires variant_table")
    command = [args.rscript or __import__("os").environ.get("RSCRIPT", "Rscript"), str(REPO_ROOT / "qc_analysis/scripts/run_intraspecies_contamination.R"), "--variant-table", str(variant_table), "--outdir", str(outdir)]
    if negative_controls: command += ["--negative-control-pairs", str(negative_controls)]
    if overwrite: command.append("--overwrite")
    for name in sorted(R_PARAMETERS & values.keys()):
        value = values[name]
        if value is not None:
            command += ["--parameter", f"{name}={str(value).lower() if isinstance(value, bool) else value}"]
    subprocess.run(command, check=True)
    return 0


def main() -> int:
    try:
        return run(parser().parse_args())
    except (ValueError, subprocess.CalledProcessError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2

if __name__ == "__main__":
    raise SystemExit(main())
