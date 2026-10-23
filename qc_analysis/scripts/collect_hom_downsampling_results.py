#!/usr/bin/env python3
"""Collect independently named HOM pilot outputs for downstream QC analyses.

On McCleary:
  python3 qc_analysis/scripts/collect_hom_downsampling_results.py
  python3 qc_analysis/scripts/collect_hom_downsampling_results.py --dry-run
  python3 qc_analysis/scripts/collect_hom_downsampling_results.py --mode symlink

Read-only with respect to source results. Files are copied atomically (or linked)
to unique names including sample, raw depth, seed and consensus DP.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import sys
import tempfile
from collections import Counter
from pathlib import Path

DEFAULT_ROOT = Path(
    "/vast/palmer/scratch/lake_nicole/lt692/"
    "primate_mtDNA_analysis_coverage_only/results/qc_coverage_only/downsampling_pilot"
)
RUN_RE = re.compile(r"^ds(?P<depth>\d+)x_seed(?P<seed>\d+)(?:_cdp(?P<cdp>\d+))?$")

# Exact, published end products, not any WDL/Cromwell intermediate with the
# same basename. Optional files are useful for diagnosing liftback and NUMTs.
FILES = [
    ("final_vcf", "round_2_variant_calling_original_coords", "round2.original_coords.clean.final.split.vcf.gz", True),
    ("final_vcf_index", "round_2_variant_calling_original_coords", "round2.original_coords.clean.final.split.vcf.gz.tbi", True),
    ("per_base_coverage", "round_2_variant_calling_original_coords", "round2.original_coords.per_base_coverage.tsv", True),
    ("mtcn", "round_2/mtcn", "round2.mtcn.tsv", True),
    ("liftover_summary", "round_2_variant_calling_original_coords", "round2.original_coords.liftback.summary.tsv", False),
    ("coverage_liftover_summary", "round_2_variant_calling_original_coords", "round2.original_coords.coverage_liftback.summary.tsv", False),
    ("merge_policy_summary", "round_2_variant_calling_original_coords", "round2.original_coords.merge_policy.summary.tsv", False),
    ("unresolved_liftover", "round_2_variant_calling_original_coords", "round2.original_coords.unresolved.tsv", False),
    ("skipped_consensus_overlap", "round_2_variant_calling_original_coords", "round2.original_coords.skipped_overlap_consensus_sites.tsv", False),
    ("round2_dropped_by_round1", "round_2_variant_calling_original_coords", "round2.original_coords.round2_dropped_by_round1.tsv", False),
    ("stacked_indel_filter_summary", "round_2_variant_calling_original_coords", "round2.original_coords.stacked_indel_filter.summary.tsv", False),
    ("nuclear_coverage_qc", "round_2/mtcn", "round2.nuclear_coverage_qc.tsv", False),
    ("nuclear_reference_validation", "round_2/mtcn", "round2.nuclear_reference_validation.tsv", False),
    ("nuclear_cov_method", "round_2/mtcn", "round2.nuclear_cov_method.tsv", False),
]

def read_manifest(path: Path):
    if not path.is_file():
        raise FileNotFoundError(f"Missing manifest: {path}")
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if not rows:
        raise ValueError(f"Empty manifest: {path}")
    required = {"sample", "species", "target_depth", "seed"}
    if not required.issubset(rows[0]):
        raise ValueError(f"Missing manifest columns: {required - set(rows[0])}")
    return rows

def publish(src: Path, dst: Path, mode: str, dry_run: bool):
    if dry_run:
        return "would_collect"
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        if mode == "symlink":
            if dst.is_symlink() and dst.resolve() == src.resolve():
                return "already_collected"
        else:
            if dst.is_file() and dst.stat().st_size == src.stat().st_size:
                # Do not overwrite anything of the same name, even if sizes match.
                return "already_collected"
        raise RuntimeError(f"Destination conflict (refusing overwrite): {dst}")
    if mode == "symlink":
        os.symlink(str(src.resolve()), str(dst))
        return "collected"
    fd, tmp = tempfile.mkstemp(prefix=f".{dst.name}.", dir=str(dst.parent))
    try:
        with os.fdopen(fd, "wb") as out, src.open("rb") as inp:
            shutil.copyfileobj(inp, out, length=1024 * 1024)
        os.replace(tmp, dst)
    finally:
        if os.path.lexists(tmp):
            os.unlink(tmp)
    return "collected"

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base", type=Path, default=DEFAULT_ROOT)
    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--mode", choices=["copy", "symlink"], default="copy")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--only-cdp", type=int, default=None,
                    help="Restrict collection to one consensus DP, e.g. 100")
    args = ap.parse_args()
    base = args.base.resolve()
    runs_root = base / "variant_calling_pilot"
    out = (args.output or base / "collected_hom_results").resolve()
    manifest = read_manifest(base / "candidate_downsampling_jobs.tsv")
    jobs = {(x["sample"], int(x["target_depth"]), int(x["seed"])): x
            for x in manifest}
    if len(jobs) != len(manifest):
        raise ValueError("Duplicate combinations in candidate manifest")
    if not runs_root.is_dir():
        raise FileNotFoundError(f"Missing calling results: {runs_root}")
    versions = []
    for run in sorted(runs_root.iterdir()):
        if not run.is_dir():
            continue
        match = RUN_RE.fullmatch(run.name)
        if not match:
            continue
        depth = int(match.group("depth"))
        seed = int(match.group("seed"))
        # The original untagged experiment ran at DP100.
        cdp = int(match.group("cdp") or 100)
        if args.only_cdp is not None and cdp != args.only_cdp:
            continue
        versions.append((run, depth, seed, cdp))
    if not versions:
        raise RuntimeError(f"No matching ds* calling directories in {runs_root}")

    inventory = []
    overview = []
    seen = set()
    for run, depth, seed, cdp in versions:
        root = run / "round2"
        if not root.exists():
            continue
        for sample_dir in sorted(root.iterdir()):
            if not sample_dir.is_dir():
                continue
            sample = sample_dir.name
            key = (sample, depth, seed)
            if key not in jobs:
                continue
            trial = (sample, depth, seed, cdp)
            if trial in seen:
                raise RuntimeError(f"Two directories map to same experimental condition: {trial}")
            seen.add(trial)
            species = jobs[key]["species"]
            prefix = f"{sample}__ds{depth}x__seed{seed}__cdp{cdp}"
            missing = []
            exported = 0
            for kind, subdir, suffix, required in FILES:
                src = sample_dir / subdir / (sample + "." + suffix)
                dst = out / kind / (prefix + "." + suffix)
                if not src.is_file() or (required and src.stat().st_size == 0):
                    state = "MISSING" if required else "not_available"
                    if required:
                        missing.append(kind)
                else:
                    state = publish(src, dst, args.mode, args.dry_run)
                    exported += 1
                inventory.append({
                    "sample": sample, "species": species,
                    "target_depth": depth, "seed": seed, "consensus_min_dp": cdp,
                    "file_type": kind, "required": int(required),
                    "status": state, "source": str(src),
                    "collected_path": str(dst) if state not in ("MISSING", "not_available") else "",
                    "size_bytes": src.stat().st_size if src.is_file() else "",
                })
            # Retain sampling metrics and explicit consensus-parameter receipts.
            supplementary = [
                ("downsampling_metrics", base / "candidate_downsampled" / sample /
                 f"{sample}.candidate.ds{depth}x.seed{seed}.metrics.tsv"),
                ("round1_parameters", run / "round1" / "pilot_consensus_parameters.tsv"),
                ("round2_parameters", run / "round2" / "pilot_consensus_parameters.tsv"),
            ]
            for kind, src in supplementary:
                dst = out / kind / (prefix + "." + kind + ".tsv")
                if not src.is_file():
                    state = "not_available"
                else:
                    state = publish(src, dst, args.mode, args.dry_run)
                    exported += 1
                inventory.append({
                    "sample": sample, "species": species,
                    "target_depth": depth, "seed": seed, "consensus_min_dp": cdp,
                    "file_type": kind, "required": 0, "status": state,
                    "source": str(src),
                    "collected_path": str(dst) if state != "not_available" else "",
                    "size_bytes": src.stat().st_size if src.is_file() else "",
                })
            overview.append({
                "sample": sample, "species": species,
                "target_depth": depth, "seed": seed, "consensus_min_dp": cdp,
                "required_files_complete": "YES" if not missing else "NO",
                "missing_required_types": ",".join(missing),
                "files_available": exported, "round2_directory": str(sample_dir)
            })
    if not overview:
        raise RuntimeError("No completed or partially completed Round2 sample folders found")

    # The inventory includes missing/optional files; it can be re-run safely
    # after jobs finish. Atomic TSV replacements prevent truncated manifests.
    def write_tsv(name, rows):
        if args.dry_run:
            return
        out.mkdir(parents=True, exist_ok=True)
        target = out / name
        fd, temp = tempfile.mkstemp(prefix=f".{name}.", dir=out)
        try:
            with os.fdopen(fd, "w", newline="") as stream:
                wr = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
                wr.writeheader()
                wr.writerows(rows)
            os.replace(temp, target)
        finally:
            if os.path.lexists(temp):
                os.unlink(temp)

    write_tsv("files_manifest.tsv", inventory)
    write_tsv("experiment_overview.tsv", overview)
    counts = Counter(x["required_files_complete"] for x in overview)
    by_cdp = Counter(x["consensus_min_dp"] for x in overview)
    print(f"Found {len(overview)} sample/depth/seed/consensus-DP runs; DP breakdown: {dict(sorted(by_cdp.items()))}")
    print(f"Required outputs complete: {counts['YES']}; incomplete: {counts['NO']}")
    print(f"Files inventoried: {len(inventory)}; collection mode: {args.mode}; dry_run={args.dry_run}")
    print(f"Output: {out}")
    if counts["NO"]:
        for x in overview:
            if x["required_files_complete"] == "NO":
                print(f"INCOMPLETE: {x['sample']} {x['target_depth']}x seed{x['seed']} cdp{x['consensus_min_dp']}: {x['missing_required_types']}")
    return 0

if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
