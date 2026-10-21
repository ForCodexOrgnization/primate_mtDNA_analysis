#!/usr/bin/env python3
"""Downsample the existing Round-1 chrM+NUMT candidate BAM by read template.

Two modes:
  prepare: measure raw chrM median depth once per candidate BAM, make 60-job TSV.
  run:     downsample one manifest row, measure actual depth, and write metrics.

Requires samtools >=1.15 on PATH. Does not change original WGS CRAM or candidate BAM.
This is stage-1 data production, NOT Round-1/Round-2 variant calling.
"""
import argparse
import csv
import os
from pathlib import Path
import statistics
import subprocess
import sys

ROOT = Path(os.getenv("PILOT_BASE", "/vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis_coverage_only/results/qc_coverage_only/downsampling_pilot"))
DEFAULT_CANDIDATE_ROOT = Path(os.getenv("CANDIDATE_ROOT", "/vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis_coverage_only/results"))
FIELDS = ["sample", "species", "input_bam", "mt_contig", "original_raw_median", "target_depth", "seed", "fraction"]


def command(*args):
    print("+", " ".join(map(str, args)), flush=True)
    subprocess.run([str(x) for x in args], check=True)


def locate(sample, candidate_root):
    # Canonical publishDir of the upstream EXTRACT_CANDIDATE_READS process.
    possibilities = [
        candidate_root / sample / "round_1" / "candidate_reads" / f"{sample}.with_mates.bam",
        candidate_root / sample / f"{sample}.with_mates.bam",
        candidate_root / f"{sample}.with_mates.bam",
    ]
    found = [p for p in possibilities if p.is_file()]
    if len(found) != 1:
        raise RuntimeError(
            f"{sample}: expected exactly one candidate BAM in known layout; found {len(found)}. "
            "Supply --candidate-root pointing to the pipeline outdir or build a symlink. "
            "Checked: " + ", ".join(map(str, possibilities))
        )
    return found[0].resolve()


def depth_median(bam, contig):
    # Include all reference bases (including zeros); use same tool/settings for
    # source and destination. Do not confuse this with post-calling DP.
    proc = subprocess.Popen(
        ["samtools", "depth", "-aa", "-r", contig, str(bam)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    vals = []
    for row in proc.stdout:
        columns = row.rstrip("\n").split("\t")
        if len(columns) != 3 or columns[0] != contig:
            proc.kill()
            raise RuntimeError(f"Unexpected samtools depth output: {row[:160]}")
        vals.append(int(columns[2]))
    err = proc.stderr.read()
    if proc.wait() != 0:
        raise RuntimeError(f"samtools depth failed on {bam}: {err}")
    if not vals:
        raise RuntimeError(f"No depth values for contig {contig} in {bam}")
    return statistics.median(vals), vals


def idxstats(bam, contig):
    proc = subprocess.run(["samtools", "idxstats", str(bam)], capture_output=True, text=True, check=True)
    m = n = 0
    contigs = set()
    for row in proc.stdout.splitlines():
        c, length, mapped, unmapped = row.split("\t")
        contigs.add(c)
        if c == contig:
            m += int(mapped)
        elif c != "*":
            n += int(mapped)
    if contig not in contigs:
        raise RuntimeError(f"{contig} not in BAM reference dictionary: {bam}")
    return m, n


def prepare(args):
    candidate_root = args.candidate_root
    source = ROOT / "pilot_samples.tsv"
    if not source.is_file():
        raise RuntimeError(f"Missing pilot sample table: {source}")
    rows = []
    with source.open(newline="") as fp:
        samples = list(csv.DictReader(fp, delimiter="\t"))
    if not samples:
        raise RuntimeError("No samples")
    for rec in samples:
        sample = rec["sample"]
        bam = locate(sample, candidate_root)
        if not Path(str(bam) + ".bai").is_file():
            command("samtools", "index", "-@", str(args.threads), str(bam))
        mt = args.mt_contig  # pilot source BAMs should have species-verified contig.
        before, _ = depth_median(bam, mt)
        if before < max(args.targets):
            raise RuntimeError(f"{sample}: original raw mt depth {before} below requested maximum {max(args.targets)}")
        print(f"CALIBRATED sample={sample} chrM={mt} raw_median={before}", flush=True)
        for target in args.targets:
            for seed in args.seeds:
                fraction = min(1.0, target / before)
                rows.append({
                    "sample": sample, "species": rec["species"], "input_bam": str(bam),
                    "mt_contig": mt, "original_raw_median": before,
                    "target_depth": target, "seed": seed, "fraction": f"{fraction:.12f}",
                })
    ROOT.mkdir(parents=True, exist_ok=True)
    # Write the whole manifest only after all samples calibrate successfully.
    manifest = ROOT / "candidate_downsampling_jobs.tsv"
    tmp = manifest.with_suffix(".tsv.tmp")
    with tmp.open("w", newline="") as fp:
        writer = csv.DictWriter(fp, delimiter="\t", fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, manifest)
    print(f"PREPARED {len(rows)} jobs: {manifest}", flush=True)


def run(args):
    manifest = ROOT / "candidate_downsampling_jobs.tsv"
    with manifest.open(newline="") as fp:
        jobs = list(csv.DictReader(fp, delimiter="\t"))
    if args.task < 1 or args.task > len(jobs):
        raise RuntimeError(f"Invalid task {args.task}; valid 1..{len(jobs)}")
    job = jobs[args.task - 1]
    sample, target, seed = job["sample"], job["target_depth"], job["seed"]
    fraction = float(job["fraction"])
    if not (0 < fraction <= 1):
        raise RuntimeError(f"Invalid fraction {fraction}")
    source = Path(job["input_bam"])
    if not source.is_file():
        raise RuntimeError(f"Missing candidate BAM {source}")
    mt = job["mt_contig"]
    dest = ROOT / "candidate_downsampled" / sample
    dest.mkdir(parents=True, exist_ok=True)
    prefix = f"{sample}.candidate.ds{target}x.seed{seed}"
    out = dest / f"{prefix}.bam"
    metrics = dest / f"{prefix}.metrics.tsv"
    if not out.is_file() or not Path(str(out) + ".bai").is_file():
        temp = dest / f"{prefix}.partial.{os.getpid()}.bam"
        try:
            command("samtools", "view", "-@", str(args.threads), "-b",
                    "--subsample-seed", seed, "--subsample", job["fraction"],
                    "-o", str(temp), str(source))
            command("samtools", "quickcheck", "-v", str(temp))
            command("samtools", "index", "-@", str(args.threads), str(temp))
            os.replace(temp, out)
            os.replace(Path(str(temp) + ".bai"), Path(str(out) + ".bai"))
        finally:
            temp.unlink(missing_ok=True)
            Path(str(temp) + ".bai").unlink(missing_ok=True)
    command("samtools", "quickcheck", "-v", str(out))
    measured, vals = depth_median(out, mt)
    mt_before, nonmt_before = idxstats(source, mt)
    mt_after, nonmt_after = idxstats(out, mt)
    def fraction_num(n):
        return n / len(vals)
    row = {
        **job, "actual_raw_median": measured,
        "actual_raw_mean": f"{statistics.mean(vals):.4f}",
        "Percent_20": f"{100*fraction_num(sum(x>=20 for x in vals)):.5f}",
        "Percent_30": f"{100*fraction_num(sum(x>=30 for x in vals)):.5f}",
        "Percent_50": f"{100*fraction_num(sum(x>=50 for x in vals)):.5f}",
        "Percent_75": f"{100*fraction_num(sum(x>=75 for x in vals)):.5f}",
        "Percent_100": f"{100*fraction_num(sum(x>=100 for x in vals)):.5f}",
        "mt_mapped_before": mt_before, "nonmt_mapped_before": nonmt_before,
        "mt_mapped_after": mt_after, "nonmt_mapped_after": nonmt_after,
        "bam_path": str(out),
    }
    cols = FIELDS + [c for c in row if c not in FIELDS]
    tmp = metrics.with_suffix(".tsv.tmp")
    with tmp.open("w", newline="") as fp:
        writer = csv.DictWriter(fp, delimiter="\t", fieldnames=cols)
        writer.writeheader()
        writer.writerow(row)
    os.replace(tmp, metrics)
    print(f"DONE {sample} target={target}x seed={seed} raw_median={measured} output={out}", flush=True)


def main():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="mode", required=True)
    q = sub.add_parser("prepare")
    q.add_argument("--candidate-root", type=Path, default=DEFAULT_CANDIDATE_ROOT)
    q.add_argument("--mt-contig", default="chrM")
    q.add_argument("--targets", type=int, nargs="+", default=[30, 50, 100])
    q.add_argument("--seeds", type=int, nargs="+", default=[101, 202])
    q.add_argument("--threads", type=int, default=4)
    q = sub.add_parser("run")
    q.add_argument("--task", type=int, required=True)
    q.add_argument("--threads", type=int, default=4)
    args = p.parse_args()
    try:
        if args.mode == "prepare":
            prepare(args)
        else:
            run(args)
    except (RuntimeError, subprocess.CalledProcessError, OSError, KeyError, ValueError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
