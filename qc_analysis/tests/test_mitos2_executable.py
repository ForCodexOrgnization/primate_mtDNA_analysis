import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'qc_analysis/scripts/run_mitos2_annotation.py'


def load_module():
    spec = importlib.util.spec_from_file_location('run_mitos2_annotation', SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_mitos2_uses_runmitos_after_conda_activation(monkeypatch):
    module = load_module()
    captured = {}

    class Result:
        returncode = 0
        stdout = 'CONDA_PREFIX=/tmp/mitos2\nMITOS executable=/tmp/mitos2/bin/runmitos\nUsing MITOS2 executable: /tmp/mitos2/bin/runmitos\n'
        stderr = ''

    def fake_run(args, **kwargs):
        captured['command'] = args[2]
        return Result()

    monkeypatch.setattr(module.subprocess, 'run', fake_run)
    executable, output = module.command({'conda_module': 'miniconda', 'conda_env': 'mitos2'})

    assert executable == 'runmitos'
    assert output == Result.stdout
    assert 'module load miniconda' in captured['command']
    assert 'source "$(conda info --base)/etc/profile.d/conda.sh"' in captured['command']
    assert 'conda activate mitos2' in captured['command']
    assert 'command -v runmitos' in captured['command']
    assert 'command -v ' + 'MITOS' + '2' not in captured['command']
    assert 'command -v ' + 'runmitos' + '.py' not in captured['command']


def test_mitos2_command_templates_invoke_runmitos():
    module = load_module()
    commands = module.templates('runmitos', 'input.fa', 'output', {'genetic_code': 2, 'threads': 4})

    assert all(command.startswith('runmitos ') for command in commands)
