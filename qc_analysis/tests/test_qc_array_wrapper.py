import os, subprocess
from pathlib import Path

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
