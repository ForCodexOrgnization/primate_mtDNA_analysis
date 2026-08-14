"""Regression coverage for explicit MITOS2 interpreter selection in the wrapper."""

from pathlib import Path
import os
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[2]
WRAPPER = ROOT / "qc_analysis/scripts/run_qc_preprocessing.sh"


def write_executable(path: Path, contents: str) -> None:
    path.write_text(contents)
    path.chmod(0o755)


def test_mitos2_wrapper_uses_conda_prefix_python_when_shell_python_lacks_biopython(tmp_path):
    """The activated environment's Python must win even with a stale shell PATH."""
    shell_bin = tmp_path / "shell-bin"
    shell_bin.mkdir()
    conda_base = tmp_path / "conda-base"
    profile_dir = conda_base / "etc/profile.d"
    profile_dir.mkdir(parents=True)
    mitos_prefix = tmp_path / "mitos2"
    mitos_bin = mitos_prefix / "bin"
    mitos_bin.mkdir(parents=True)

    # This is deliberately the default python and cannot import Bio.
    write_executable(
        shell_bin / "python3",
        "#!/usr/bin/env bash\necho 'ModuleNotFoundError: No module named Bio' >&2\nexit 1\n",
    )
    write_executable(shell_bin / "python", "#!/usr/bin/env bash\nexec \"$(dirname \"$0\")/python3\" \"$@\"\n")
    write_executable(shell_bin / "module", "#!/usr/bin/env bash\nexit 0\n")
    write_executable(
        shell_bin / "conda",
        f"#!/usr/bin/env bash\n[[ \"$1 $2\" == 'info --base' ]] && printf '%s\\n' '{conda_base}'\n",
    )
    (profile_dir / "conda.sh").write_text(
        f"conda() {{\n  if [[ \"$1\" == activate ]]; then\n"
        f"    export CONDA_PREFIX='{mitos_prefix}' CONDA_DEFAULT_ENV=mitos2\n  fi\n}}\n"
    )
    write_executable(
        mitos_bin / "python",
        "#!/usr/bin/env bash\n"
        "if [[ \"$1\" == --version ]]; then echo 'Python 3.11.9'; exit 0; fi\n"
        "if [[ \"$1\" == -c ]]; then\n"
        "  if [[ \"$2\" == *'import sys'* ]]; then echo \"$0\"; fi\n"
        "  echo '1.81'; exit 0\n"
        "fi\n"
        "echo \"$@\" >> \"$MITOS2_CALL_LOG\"\n",
    )
    write_executable(mitos_bin / "runmitos", "#!/usr/bin/env bash\nexit 0\n")
    config = tmp_path / "qc.yaml"
    config.write_text(
        "mitos2_annotation:\n"
        "  settings:\n"
        "    conda_module: miniconda/test\n"
        "    conda_env: mitos2\n"
    )
    call_log = tmp_path / "mitos2-calls.log"
    env = os.environ | {
        "PATH": f"{shell_bin}:{os.environ['PATH']}",
        "PYTHON": str(shell_bin / "python3"),
        "MITOS2_CALL_LOG": str(call_log),
        # Exercise the case in which the wrapper starts in that environment.
        "CONDA_PREFIX": str(mitos_prefix),
        "CONDA_DEFAULT_ENV": "mitos2",
    }

    completed = subprocess.run(
        ["bash", str(WRAPPER), "mitos2_prepare_tasks", str(config)],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    assert f"MITOS2_PYTHON={mitos_bin / 'python'}" in completed.stderr
    assert f"command -v python={shell_bin / 'python'}" in completed.stderr
    assert "Biopython version=1.81" in completed.stderr
    assert call_log.read_text().startswith("qc_analysis/scripts/run_mitos2_annotation.py")


@pytest.mark.parametrize(
    ("step", "expected_script"),
    [
        ("build_trna_indexes", "build_all_trna_indexes.py"),
        ("trna_match", "run_trna_match.py"),
    ],
)
def test_trna_steps_activate_configured_environment_and_binary(tmp_path, step, expected_script):
    shell_bin = tmp_path / "shell-bin"
    shell_bin.mkdir()
    conda_base = tmp_path / "conda-base"
    (conda_base / "etc/profile.d").mkdir(parents=True)
    trna_prefix = tmp_path / "custom-trna-env"
    trna_bin = trna_prefix / "bin"
    trna_bin.mkdir(parents=True)
    call_log = tmp_path / "calls.log"
    module_log = tmp_path / "module.log"

    write_executable(shell_bin / "module", "#!/usr/bin/env bash\nprintf '%s\\n' \"$*\" >> \"$MODULE_LOG\"\n")
    write_executable(
        shell_bin / "conda",
        f"#!/usr/bin/env bash\n[[ \"$1 $2\" == 'info --base' ]] && printf '%s\\n' '{conda_base}'\n",
    )
    write_executable(shell_bin / "python3", "#!/usr/bin/env bash\nprintf '%s\\n' \"$*\" >> \"$CALL_LOG\"\n")
    (conda_base / "etc/profile.d/conda.sh").write_text(
        "conda() {\n"
        "  if [[ \"$1\" == activate ]]; then\n"
        f"    export CONDA_PREFIX='{trna_prefix}' CONDA_DEFAULT_ENV=bespoke_trna\n"
        f"    export PATH='{trna_bin}':\"$PATH\"\n"
        "  fi\n"
        "}\n"
    )
    write_executable(
        trna_bin / "custom-trnascan",
        "#!/usr/bin/env bash\n[[ \"$1\" == --version ]] && echo 'custom tRNAscan-SE 2.0'\n",
    )
    config = tmp_path / "qc.yaml"
    config.write_text(
        "trna_match:\n"
        "  settings:\n"
        "    conda_module: miniconda/trna-test\n"
        "    conda_env: bespoke_trna\n"
        "    trnascan_bin: custom-trnascan\n"
        "mitos2_annotation:\n"
        "  settings:\n"
        "    conda_module: miniconda/mitos-untouched\n"
        "    conda_env: mitos2\n"
    )
    env = os.environ | {
        "PATH": f"{shell_bin}:{os.environ['PATH']}",
        "PYTHON": str(shell_bin / "python3"),
        "CALL_LOG": str(call_log),
        "MODULE_LOG": str(module_log),
    }

    completed = subprocess.run(
        ["bash", str(WRAPPER), step, str(config)], cwd=ROOT, env=env,
        text=True, capture_output=True, check=False,
    )

    assert completed.returncode == 0, completed.stderr
    assert module_log.read_text() == "load miniconda/trna-test\n"
    assert "CONDA_DEFAULT_ENV=bespoke_trna" in completed.stderr
    assert f"tRNAscan-SE executable={trna_bin / 'custom-trnascan'}" in completed.stderr
    assert "tRNAscan-SE version=custom tRNAscan-SE 2.0" in completed.stderr
    assert expected_script in call_log.read_text()


def test_missing_configured_trnascan_binary_fails_before_python(tmp_path):
    shell_bin = tmp_path / "shell-bin"
    shell_bin.mkdir()
    conda_base = tmp_path / "conda-base"
    (conda_base / "etc/profile.d").mkdir(parents=True)
    prefix = tmp_path / "empty-env"
    prefix.mkdir()
    call_log = tmp_path / "calls.log"
    write_executable(shell_bin / "module", "#!/usr/bin/env bash\nexit 0\n")
    write_executable(
        shell_bin / "conda",
        f"#!/usr/bin/env bash\n[[ \"$1 $2\" == 'info --base' ]] && echo '{conda_base}'\n",
    )
    write_executable(shell_bin / "python3", "#!/usr/bin/env bash\necho called >> \"$CALL_LOG\"\n")
    (conda_base / "etc/profile.d/conda.sh").write_text(
        f"conda() {{ [[ \"$1\" == activate ]] && export CONDA_PREFIX='{prefix}' CONDA_DEFAULT_ENV=trnascan_env; }}\n"
    )
    config = tmp_path / "qc.yaml"
    config.write_text(
        "trna_match:\n  settings:\n    conda_module: miniconda/test\n"
        "    conda_env: trnascan_env\n    trnascan_bin: definitely-missing-trnascan\n"
    )
    completed = subprocess.run(
        ["bash", str(WRAPPER), "build_trna_indexes", str(config)], cwd=ROOT,
        env=os.environ | {"PATH": f"{shell_bin}:{os.environ['PATH']}", "PYTHON": str(shell_bin / "python3"), "CALL_LOG": str(call_log)},
        text=True, capture_output=True, check=False,
    )

    assert completed.returncode == 127
    assert "tRNAscan-SE was not found after activating conda environment: trnascan_env" in completed.stderr
    assert "Configured executable: definitely-missing-trnascan" in completed.stderr
    assert not call_log.exists()
