"""Regression coverage for explicit MITOS2 interpreter selection in the wrapper."""

from pathlib import Path
import os
import subprocess


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
