#!/usr/bin/env python3
"""Merge one-row per-sample rRNA match summaries into a cohort summary."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from qc_analysis.scripts.merge_sample_summaries import merge_sample_summaries,run_merge_cli


def merge_summaries(reports_dir,output):
    return merge_sample_summaries(
        reports_dir,output,'*.rrna_match_summary.tsv',
        'all_samples.rrna_match_summary.tsv',label='rRNA match'
    )


def main():
    return run_merge_cli(
        'rrna_match','*.rrna_match_summary.tsv',
        'all_samples.rrna_match_summary.tsv','rRNA match'
    )


if __name__ == '__main__':
    raise SystemExit(main())
