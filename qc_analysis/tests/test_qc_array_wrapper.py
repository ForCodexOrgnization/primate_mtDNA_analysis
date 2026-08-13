import os, subprocess
from pathlib import Path

import pytest
from qc_analysis.scripts.qc_array_manifest import resolve_runtime_paths

ROOT=Path(__file__).resolve().parents[2]
WRAPPER=ROOT/'qc_analysis/scripts/run_qc_preprocessing.sh'

@pytest.mark.parametrize('step,expected', [
    ('collect_variant_calling_results','results/qc/variant_calling_collection'),
    ('discover_global_anchor','results/qc/coordinate_liftover/global_anchor'),
    ('coordinate_liftover','results/qc/coordinate_liftover'),
    ('mitos2_prepare_tasks','results/qc/mitos2_annotation'),
    ('mitos2_annotation','results/qc/mitos2_annotation'),
    ('mitos2_merge','results/qc/mitos2_annotation'),
    ('build_primate_codon_table','results/qc/codon_table_build'),
    ('compare_genbank_mitos2','results/qc/genbank_mitos2_comparison'),
    ('codon_match_validate','results/qc/codon_match'), ('codon_match','results/qc/codon_match'),
    ('codon_match_merge','results/qc/codon_match'), ('trna_match','results/qc/trna_match'),
    ('rrna_match','results/qc/rrna_match'), ('intraspecies_contamination','results/qc/intraspecies_contamination'),
    ('sample_variant_filtering','results/qc/sample_variant_filtering'),
])
def test_documented_step_metadata_defaults(step, expected):
    resolved=resolve_runtime_paths(step,{})
    assert resolved['output_dir']==expected
    assert resolved['job_array_dir']==f'{expected}/job_arrays'
    assert resolved['log_dir']==f'{expected}/logs/job_arrays'

def run(*args, env=None):
    return subprocess.run(['bash',str(WRAPPER),*args],cwd=ROOT,text=True,capture_output=True,env=os.environ|dict(env or {}))

def config(tmp_path):
    p=tmp_path/'config with spaces.yaml';p.write_text(
        f'genbank_mitos2_comparison:\n  paths:\n    output_dir: {tmp_path}/comparison output\n  enabled: true\n')
    return p

def test_singleton_and_concurrency_validation(tmp_path):
    c=config(tmp_path);x=run('--dry-run-submit','compare_genbank_mitos2',str(c))
    assert x.returncode==0 and '--array=1-1' in x.stdout
    for value in ('0','-1','abc'):
        x=run('--dry-run-submit','compare_genbank_mitos2',str(c),env={'SLURM_ARRAY_CONCURRENCY':value})
        assert x.returncode==2 and 'positive integer' in x.stderr

def test_sample_variant_filtering_singleton_manifest_and_dry_run(tmp_path):
    c=config(tmp_path)
    result=run('--dry-run-submit','sample_variant_filtering',str(c))
    assert result.returncode==0, result.stderr
    assert '--array=1-1' in result.stdout
    manifests=list((ROOT/'results/qc/sample_variant_filtering/job_arrays').glob('sample_variant_filtering.*.manifest.tsv'))
    assert manifests
    assert '\tsingleton\tsample_variant_filtering\t' in manifests[-1].read_text()

def test_submit_all_preserves_initial_dependency_chain(tmp_path):
    result=run('--dry-run-submit','all',str(config(tmp_path)))
    # This deliberately minimal config may stop at a later sample-array node,
    # but the requested global chain must be constructed without an unsupported
    # sample_variant_filtering error.
    assert 'Unsupported array step' not in result.stderr
    lines=[line for line in result.stdout.splitlines() if line.startswith('DRY RUN:')]
    wanted=['collect_variant_calling_results','intraspecies_contamination','sample_variant_filtering','discover_global_anchor']
    positions=[]
    for step in wanted:
        positions.append(next(i for i,line in enumerate(lines) if f'qc_preprocessing_{step}' in line))
    assert positions==sorted(positions)
    for previous,current in zip(wanted,wanted[1:]):
        line=next(line for line in lines if f'qc_preprocessing_{current}' in line)
        assert f'afterok:dry_{previous}' in line

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


@pytest.mark.parametrize(
    "step,item,expected_sample",
    [
        ("build_primate_codon_table", "build_primate_codon_table", None),
        ("codon_match", "SAMPLE1", "SAMPLE1"),
        ("coordinate_liftover", "SAMPLE1", "SAMPLE1"),
        ("human_contamination", "human_contamination", None),
        ("codon_match_validate", "codon_match_validate", None),
        ("build_primate_homo_background", "build_primate_homo_background", None),
    ],
)
def test_array_item_is_a_sample_only_for_sample_classified_steps(
    tmp_path, step, item, expected_sample
):
    tasks = tmp_path / "tasks.txt"
    tasks.write_text(f"{item}\n")
    c = config(tmp_path)
    python = tmp_path / "python"
    call_log = tmp_path / "calls.log"
    python.write_text(
        "#!/usr/bin/env bash\n"
        "[[ $# -eq 1 && $1 == - ]] && cat >/dev/null\n"
        "printf '%s\\n' \"$*\" >> \"$CALL_LOG\"\n"
    )
    python.chmod(0o755)

    result = run(
        "--array-task", "--task-file", str(tasks), step, str(c),
        env={
            "SLURM_ARRAY_TASK_ID": "1",
            "PYTHON": str(python),
            "CALL_LOG": str(call_log),
            "BIOPYTHON_USE_MODULE": "0",
        },
    )

    assert result.returncode == 0, result.stderr
    calls = call_log.read_text().splitlines()
    if expected_sample:
        assert any(f"--sample {expected_sample}" in call for call in calls)
    else:
        assert all("--sample" not in call for call in calls)


def test_mitos2_reference_array_item_is_not_forwarded_as_a_sample(tmp_path):
    tasks = tmp_path / "tasks.txt"
    tasks.write_text("reference:NC_012920.1\n")
    conda_base = tmp_path / "conda"
    (conda_base / "etc/profile.d").mkdir(parents=True)
    prefix = tmp_path / "mitos2"
    (prefix / "bin").mkdir(parents=True)
    call_log = tmp_path / "calls.log"
    python = prefix / "bin/python"
    python.write_text(
        "#!/usr/bin/env bash\n"
        "if [[ $1 == --version ]]; then echo 'Python 3'; exit; fi\n"
        "if [[ $1 == -c ]]; then echo '1.83'; exit; fi\n"
        "printf '%s\\n' \"$*\" >> \"$CALL_LOG\"\n"
    )
    python.chmod(0o755)
    bindir = tmp_path / "bin"
    bindir.mkdir()
    module = bindir / "module"
    module.write_text("#!/usr/bin/env bash\nexit 0\n")
    module.chmod(0o755)
    conda = bindir / "conda"
    conda.write_text(
        f"#!/usr/bin/env bash\n[[ $1 == info ]] && printf '%s\\n' '{conda_base}'\n"
    )
    conda.chmod(0o755)
    (conda_base / "etc/profile.d/conda.sh").write_text(
        f"conda() {{ [[ $1 == activate ]] && export CONDA_PREFIX='{prefix}' CONDA_DEFAULT_ENV=mitos2; }}\n"
    )
    c = tmp_path / "mitos.yaml"
    c.write_text(
        "mitos2_annotation:\n  settings:\n    conda_module: test\n    conda_env: mitos2\n"
    )

    result = run(
        "--array-task", "--task-file", str(tasks), "mitos2_annotation", str(c),
        env={
            "SLURM_ARRAY_TASK_ID": "1",
            "PATH": f"{bindir}:{os.environ['PATH']}",
            "CALL_LOG": str(call_log),
        },
    )

    assert result.returncode == 0, result.stderr
    call = call_log.read_text().strip()
    assert "--reference NC_012920.1" in call
    assert "--sample" not in call


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


@pytest.mark.parametrize('step,section,expected', [
    ('coordinate_liftover', 'coordinate_liftover', 'coordinate_liftover'),
    ('codon_match', 'codon_match', 'codon_match'),
    ('trna_match', 'trna_match', 'trna_match'),
    ('rrna_match', 'rrna_match', 'rrna_match'),
])
def test_step_specific_metadata_and_log_directories(tmp_path, step, section, expected):
    samples=tmp_path/'samples.tsv'; samples.write_text('sample\nS1\n')
    output=tmp_path/expected; inputs=tmp_path/'inputs'; inputs.mkdir(); (inputs/'S1.vcf').write_text('input\n')
    paths=f'    output_dir: {output}\n'
    settings=''
    if step=='coordinate_liftover': paths+=f'    sample_ref_file: {samples}\n'
    elif step=='codon_match':
        paths+=f'    sample_reference_map: {samples}\n    input_vcf_dir: {inputs}\n'
        settings='  settings:\n    input_vcf_pattern: "{sample}.vcf"\n    output_suffix: .out.vcf\n'
    elif step=='trna_match':
        paths+=f'    input_vcf_dir: {inputs}\n    fallback_input_vcf_dir: {inputs}\n'
        settings='  settings:\n    input_vcf_pattern: "{sample}.vcf"\n    fallback_input_vcf_pattern: "{sample}.vcf"\n    output_suffix: .out.vcf\n'
    else:
        paths+=f'    input_vcf_dir: {inputs}\n    fallback_codon_vcf_dir: {inputs}\n    fallback_raw_vcf_dir: {inputs}\n'
        settings='  settings:\n    input_vcf_pattern: "{sample}.vcf"\n    fallback_codon_vcf_pattern: "{sample}.vcf"\n    fallback_raw_vcf_pattern: "{sample}.vcf"\n    output_suffix: .out.vcf\n'
    cfg=tmp_path/f'{step}.yaml'; cfg.write_text(f'{section}:\n  paths:\n{paths}{settings}')
    legacy=ROOT/'results/qc/job_arrays'; before=set(legacy.iterdir()) if legacy.exists() else set()
    result=run('--dry-run-submit','--sample','S1',step,str(cfg),env={'AUTO_SUBMIT_MERGE':'false'})
    assert result.returncode==0, result.stderr
    assert f'task_file={output}/job_arrays/' in result.stderr
    assert f'logs={output}/logs/job_arrays/%A_%a.{{out,err}}' in result.stderr
    assert (set(legacy.iterdir()) if legacy.exists() else set()) == before


def test_explicit_metadata_and_global_log_overrides_and_immutable_worker_path(tmp_path):
    output=tmp_path/'output with spaces'; metadata=tmp_path/'metadata with spaces'; logs=tmp_path/'logs with spaces'
    inputs=tmp_path/'inputs'; inputs.mkdir(); (inputs/'S1.vcf').write_text('input\n')
    cfg=tmp_path/'override.yaml'; cfg.write_text(
        f'codon_match:\n  paths:\n    output_dir: {output}\n    job_array_dir: {metadata}\n    input_vcf_dir: {inputs}\n'
        '  settings:\n    input_vcf_pattern: "{sample}.vcf"\n    output_suffix: .out.vcf\n')
    result=run('--dry-run-submit','--sample','S1','codon_match',str(cfg),
               env={'SLURM_LOG_DIR':str(logs),'AUTO_SUBMIT_MERGE':'false'})
    assert result.returncode==0, result.stderr
    task=next(metadata.glob('codon_match.*.tasks.txt'))
    assert f'task_file={task}' in result.stderr
    assert f'logs={logs}/%A_%a.{{out,err}}' in result.stderr
    # The shell-escaped dry-run command embeds the timestamped file, not current.tsv.
    assert str(task).replace(' ', '\\ ') in result.stdout
    assert 'current.tsv' not in result.stdout


def test_retry_manifest_stays_with_step_metadata(tmp_path):
    output=tmp_path/'codon'; cfg=tmp_path/'retry.yaml'; inputs=tmp_path/'inputs'; inputs.mkdir(); (inputs/'S1.vcf').write_text('input\n')
    cfg.write_text(f'codon_match:\n  paths:\n    output_dir: {output}\n    input_vcf_dir: {inputs}\n'
                   '  settings:\n    input_vcf_pattern: "{sample}.vcf"\n    output_suffix: .out.vcf\n')
    result=run('--dry-run-submit','--prepare-retry','--sample','S1','codon_match',str(cfg),
               env={'AUTO_SUBMIT_MERGE':'false'})
    assert result.returncode==0, result.stderr
    assert list((output/'job_arrays').glob('codon_match.*.retry.tasks.txt'))
    assert list((output/'job_arrays').glob('codon_match.*.retry.manifest.tsv'))


def test_direct_execution_creates_no_array_metadata(tmp_path):
    cfg=tmp_path/'direct.yaml'; output=tmp_path/'direct-output'; python=tmp_path/'python'; calls=tmp_path/'calls'
    cfg.write_text(f'codon_match:\n  paths:\n    output_dir: {output}\n')
    python.write_text('#!/usr/bin/env bash\nprintf "%s\\n" "$*" > "$CALL_LOG"\n'); python.chmod(0o755)
    result=run('--sample','S1','codon_match',str(cfg),env={'PYTHON':str(python),'CALL_LOG':str(calls)})
    assert result.returncode==0, result.stderr
    assert not (output/'job_arrays').exists()
