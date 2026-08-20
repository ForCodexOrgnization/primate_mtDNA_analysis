import csv
import subprocess
import sys
from pathlib import Path

import pytest

from qc_analysis.scripts.merge_trna_match_summaries import merge_summaries


ROOT=Path(__file__).resolve().parents[2]
SCRIPT=ROOT/'qc_analysis/scripts/merge_trna_match_summaries.py'
HEADER='sample\ttotal_records\tn_trna_id_match\n'


def write_summary(directory,name,row=None,header=HEADER,extra_rows=()):
    path=directory/f'{name}.trna_match_summary.tsv'
    text=header
    if row is not None:
        text+='\t'.join(row)+'\n'
    text+=''.join('\t'.join(item)+'\n' for item in extra_rows)
    path.write_text(text)
    return path


def read_rows(path):
    with path.open(newline='') as handle:
        return list(csv.reader(handle,delimiter='\t'))


def test_merges_two_summaries_with_one_header_in_sample_order(tmp_path):
    write_summary(tmp_path,'zeta',('zeta','2','2'))
    write_summary(tmp_path,'alpha',('alpha','1','1'))
    output=tmp_path/'all_samples.trna_match_summary.tsv'
    assert merge_summaries(tmp_path,output) == 2
    assert read_rows(output) == [
        ['sample','total_records','n_trna_id_match'],
        ['alpha','1','1'],
        ['zeta','2','2'],
    ]


def test_existing_all_samples_summary_is_ignored(tmp_path):
    write_summary(tmp_path,'sample',('S1','1','1'))
    existing=tmp_path/'all_samples.trna_match_summary.tsv'
    existing.write_text('stale\ncontent\nwith\nmultiple\nrows\n')
    assert merge_summaries(tmp_path,existing) == 1
    assert read_rows(existing)[1] == ['S1','1','1']


def test_inconsistent_headers_fail_clearly(tmp_path):
    write_summary(tmp_path,'one',('S1','1','1'))
    write_summary(tmp_path,'two',('S2','1'),header='sample\ttotal_records\n')
    with pytest.raises(ValueError,match='inconsistent.*header'):
        merge_summaries(tmp_path,tmp_path/'all_samples.trna_match_summary.tsv')


@pytest.mark.parametrize('row,extra,expected',[
    (None,(),'exactly one data row'),
    (('S1','1','1'),(('S1','2','2'),),'exactly one data row'),
])
def test_zero_or_multiple_data_rows_fail_clearly(tmp_path,row,extra,expected):
    write_summary(tmp_path,'bad',row,extra_rows=extra)
    with pytest.raises(ValueError,match=expected):
        merge_summaries(tmp_path,tmp_path/'all_samples.trna_match_summary.tsv')


def test_no_per_sample_summaries_fails_clearly(tmp_path):
    with pytest.raises(ValueError,match='no per-sample'):
        merge_summaries(tmp_path,tmp_path/'all_samples.trna_match_summary.tsv')


def test_cli_resolves_reports_directory_from_config(tmp_path):
    reports=tmp_path/'reports';reports.mkdir()
    write_summary(reports,'sample',('S1','1','1'))
    config=tmp_path/'config.yaml'
    config.write_text(f'trna_match:\n  paths:\n    reports_dir: {reports}\n')
    result=subprocess.run(
        [sys.executable,str(SCRIPT),'--config',str(config)],
        cwd=ROOT,text=True,capture_output=True,
    )
    assert result.returncode == 0,result.stderr
    assert (reports/'all_samples.trna_match_summary.tsv').is_file()
