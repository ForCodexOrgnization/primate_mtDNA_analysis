import os, subprocess
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[2]
WRAPPER=ROOT/'qc_analysis/scripts/run_qc_preprocessing.sh'

def run(*args, env=None):
    return subprocess.run(['bash',str(WRAPPER),*args],cwd=ROOT,text=True,capture_output=True,env=os.environ|dict(env or {}))

def config(tmp_path):
    p=tmp_path/'config with spaces.yaml';p.write_text('genbank_mitos2_comparison:\n  enabled: true\n')
    return p

def test_singleton_and_concurrency_validation(tmp_path):
    c=config(tmp_path);x=run('--dry-run-submit','compare_genbank_mitos2',str(c))
    assert x.returncode==0 and '--array=1-1' in x.stdout
    for value in ('0','-1','abc'):
        x=run('--dry-run-submit','compare_genbank_mitos2',str(c),env={'SLURM_ARRAY_CONCURRENCY':value})
        assert x.returncode==2 and 'positive integer' in x.stderr

def test_single_sample_is_singleton_and_paths_are_quoted(tmp_path):
    c=config(tmp_path);c.write_text(f'coordinate_liftover:\n  paths:\n    output_dir: {tmp_path}/out\n')
    x=run('--dry-run-submit','--sample','S1','coordinate_liftover',str(c))
    assert x.returncode==0 and '--array=1-1' in x.stdout
    assert 'config\\ with\\ spaces.yaml' in x.stdout

def test_array_expression_default_and_override(tmp_path):
    samples=tmp_path/'samples.tsv';samples.write_text('sample\tspecies\n'+''.join(f'S{i}\tsp\n' for i in range(100)))
    c=tmp_path/'array.yaml';c.write_text(f'coordinate_liftover:\n  paths:\n    sample_ref_file: {samples}\n    output_dir: {tmp_path}/out\n')
    for limit, expected in [('20','--array=1-100%20'),('7','--array=1-100%7')]:
        x=run('--dry-run-submit','coordinate_liftover',str(c),env={'SLURM_ARRAY_CONCURRENCY':limit})
        assert x.returncode==0 and expected in x.stdout


@pytest.mark.parametrize(
    "internal_options",
    [
        ("--array-task", "--task-file"),
        ("--task-file", "--array-task"),
    ],
)
def test_array_worker_parses_internal_options_in_either_order(tmp_path, internal_options):
    tasks = tmp_path / "tasks.txt"
    tasks.write_text("first-sample\nselected-sample\n")
    c = config(tmp_path)
    python = tmp_path / "python"
    call_log = tmp_path / "calls.log"
    python.write_text("#!/usr/bin/env bash\nprintf '%s\\n' \"$*\" > \"$CALL_LOG\"\n")
    python.chmod(0o755)

    if internal_options[0] == "--array-task":
        options = ["--array-task", "--task-file", str(tasks)]
    else:
        options = ["--task-file", str(tasks), "--array-task"]
    x = run(
        *options,
        "codon_match",
        str(c),
        env={
            "SLURM_ARRAY_TASK_ID": "2",
            "SLURM_ARRAY_JOB_ID": "1234",
            "PYTHON": str(python),
            "CALL_LOG": str(call_log),
        },
    )

    assert x.returncode == 0, x.stderr
    assert "Usage:" not in x.stderr
    assert "step=codon_match" in x.stderr
    assert f"config={c}" in x.stderr
    assert "selected_item=selected-sample" in x.stderr
    assert "Running codon_match" in x.stderr
    assert call_log.read_text().strip() == (
        f"qc_analysis/scripts/run_codon_match.py --config {c} --sample selected-sample"
    )


@pytest.mark.parametrize("option", ["--sample", "--task-file", "--array-concurrency"])
def test_options_requiring_values_report_clear_error(option):
    x = run(option)
    assert x.returncode == 2
    assert x.stderr.startswith(f"ERROR: {option} requires a value")


def test_array_task_requires_task_file(tmp_path):
    x = run("--array-task", "codon_match", str(config(tmp_path)))
    assert x.returncode == 2
    assert "ERROR: --array-task requires --task-file" in x.stderr


def test_unknown_option_reports_clear_error():
    x = run("--not-an-option")
    assert x.returncode == 2
    assert x.stderr.startswith("ERROR: unknown option: --not-an-option\nUsage:")
