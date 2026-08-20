#!/usr/bin/env python3
"""Merge one-row per-sample tRNA match summaries into a cohort summary."""
import argparse
import csv
import os
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from qc_analysis.lib.simple_yaml import read_simple_yaml


def merge_summaries(reports_dir, output):
    reports_dir=Path(reports_dir)
    output=Path(output)
    candidates=sorted(
        path for path in reports_dir.glob('*.trna_match_summary.tsv')
        if path.name != 'all_samples.trna_match_summary.tsv'
        and path.resolve() != output.resolve()
    )
    if not candidates:
        raise ValueError(f'no per-sample tRNA match summaries found in {reports_dir}')

    expected_header=None
    sample_index=None
    rows=[]
    seen_samples={}
    for path in candidates:
        with path.open(newline='') as handle:
            records=[record for record in csv.reader(handle,delimiter='\t') if record]
        if not records:
            raise ValueError(f'per-sample summary has no header or data row: {path}')
        header=records[0]
        data=records[1:]
        if expected_header is None:
            expected_header=header
            if 'sample' not in header:
                raise ValueError(f'per-sample summary header lacks sample column: {path}')
            sample_index=header.index('sample')
        elif header != expected_header:
            raise ValueError(f'inconsistent per-sample summary header: {path}')
        if len(data) != 1:
            raise ValueError(f'per-sample summary must contain exactly one data row; found {len(data)}: {path}')
        row=data[0]
        if len(row) != len(expected_header):
            raise ValueError(
                f'per-sample summary row has {len(row)} columns; expected {len(expected_header)}: {path}'
            )
        sample=row[sample_index].strip()
        if not sample:
            raise ValueError(f'per-sample summary has a blank sample value: {path}')
        if sample in seen_samples:
            raise ValueError(
                f'duplicate sample {sample!r} in {seen_samples[sample]} and {path}'
            )
        seen_samples[sample]=path
        rows.append(row)

    rows.sort(key=lambda row: row[sample_index])
    output.parent.mkdir(parents=True,exist_ok=True)
    temporary=None
    try:
        with tempfile.NamedTemporaryFile(
            mode='w',newline='',dir=output.parent,
            prefix=output.name+'.tmp.',delete=False
        ) as handle:
            temporary=Path(handle.name)
            writer=csv.writer(handle,delimiter='\t',lineterminator='\n')
            writer.writerow(expected_header)
            writer.writerows(rows)
        os.replace(temporary,output)
    finally:
        if temporary and temporary.exists():
            temporary.unlink()
    return len(rows)


def resolve_paths(args):
    reports_dir=Path(args.reports_dir) if args.reports_dir else None
    if args.config:
        config=read_simple_yaml(Path(args.config))
        configured=(config.get('trna_match',{}).get('paths',{}) or {}).get('reports_dir')
        if not reports_dir and configured:
            reports_dir=Path(configured)
    if not reports_dir:
        raise ValueError('provide --reports-dir or a config with trna_match.paths.reports_dir')
    output=Path(args.output) if args.output else reports_dir/'all_samples.trna_match_summary.tsv'
    return reports_dir,output


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config')
    parser.add_argument('--reports-dir')
    parser.add_argument('--output')
    args=parser.parse_args()
    try:
        reports_dir,output=resolve_paths(args)
        count=merge_summaries(reports_dir,output)
    except (OSError,ValueError) as error:
        print(f'ERROR: {error}',file=sys.stderr)
        return 1
    print(f'Merged {count} tRNA match summaries into {output}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
