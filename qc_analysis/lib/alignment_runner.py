"""Environment-aware execution helpers for external sequence aligners."""
from __future__ import annotations

import os
import shlex
import shutil
import subprocess
from pathlib import Path
from typing import Sequence


def with_threads(options: Sequence[str], threads: int) -> list[str]:
    """Add MAFFT's thread option unless callers already supplied one."""
    result = list(options)
    if threads > 0 and not any(option == "--thread" or option.startswith("--thread=") for option in result):
        result.extend(["--thread", str(threads)])
    return result


def effective_thread_count(options: Sequence[str], threads: int) -> str:
    for index, option in enumerate(options):
        if option == "--thread" and index + 1 < len(options):
            return options[index + 1]
        if option.startswith("--thread="):
            return option.partition("=")[2]
    return str(threads)


def _diagnostics(aligner: str, module_load: str, conda_env: str, command_result: str = "") -> str:
    return (
        f"Unable to resolve aligner {aligner!r}. PATH={os.environ.get('PATH', '')!r}; "
        f"module={module_load or '(none)'!r}; conda_env={conda_env or '(none)'!r}; "
        f"command -v {aligner}={command_result or '(not found)'!r}. "
        "Test the requested environment with --check-environment."
    )


def _shell_environment_command(aligner: str, module_load: str, conda_env: str, command: str) -> str:
    lines = ["set -euo pipefail", "source /etc/profile >/dev/null 2>&1 || true"]
    if module_load:
        lines.append(f"module load {shlex.quote(module_load)}")
    lines.extend(["CONDA_BASE=$(conda info --base)", 'source "$CONDA_BASE/etc/profile.d/conda.sh"'])
    if conda_env:
        lines.append(f"conda activate {shlex.quote(conda_env)}")
    lines.append(f"command -v {shlex.quote(aligner)}")
    lines.append(command)
    return "\n".join(lines)


def resolve_aligner(aligner: str, use_conda_env: bool, module_load: str, conda_env: str) -> tuple[str, str]:
    """Resolve an aligner in direct, conda, then module+conda environments."""
    path = Path(aligner).expanduser()
    if (path.is_absolute() or path.parent != Path(".")) and path.is_file() and os.access(path, os.X_OK):
        return str(path), "direct"
    found = shutil.which(aligner)
    if found:
        return found, "PATH"
    if not use_conda_env:
        raise RuntimeError(_diagnostics(aligner, module_load, conda_env))
    conda = shutil.which("conda")
    if conda and conda_env:
        probe = subprocess.run([conda, "run", "--no-capture-output", "-n", conda_env, aligner, "--version"], text=True, capture_output=True)
        if probe.returncode == 0:
            return aligner, f"conda:{conda_env}"
    probe_command = _shell_environment_command(aligner, module_load, conda_env, "true")
    probe = subprocess.run(["/bin/bash", "-c", probe_command], text=True, capture_output=True)
    if probe.returncode == 0 and probe.stdout.strip():
        return probe.stdout.strip().splitlines()[-1], f"module+conda:{conda_env}"
    result = (probe.stdout + probe.stderr).strip()
    raise RuntimeError(_diagnostics(aligner, module_load, conda_env, result))


def run_aligner(aligner: str, options: list[str], input_fasta: Path, output_fasta: Path, threads: int, use_conda_env: bool, module_load: str, conda_env: str) -> tuple[str, str, list[str]]:
    """Run an aligner and redirect only its stdout to ``output_fasta``."""
    options = with_threads(options, threads)
    executable, environment = resolve_aligner(aligner, use_conda_env, module_load, conda_env)
    command = [executable, *options, str(input_fasta)]
    if environment.startswith("conda:"):
        conda = shutil.which("conda")
        command = [conda or "conda", "run", "--no-capture-output", "-n", conda_env, *command]
    with output_fasta.open("w") as output:
        try:
            if environment.startswith("module+conda:"):
                shell_command = " ".join(shlex.quote(part) for part in [aligner, *options, str(input_fasta)])
                subprocess.run(["/bin/bash", "-c", _shell_environment_command(aligner, module_load, conda_env, shell_command)], check=True, stdout=output, stderr=subprocess.PIPE, text=True)
            else:
                subprocess.run(command, check=True, stdout=output, stderr=subprocess.PIPE, text=True)
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(f"Aligner {aligner!r} failed with exit status {exc.returncode}: {exc.stderr.strip()}") from exc
    return executable, environment, command


def check_aligner_environment(aligner: str, options: list[str], threads: int, use_conda_env: bool, module_load: str, conda_env: str) -> dict[str, str]:
    executable, environment = resolve_aligner(aligner, use_conda_env, module_load, conda_env)
    command = [executable, "--version"]
    if environment.startswith("conda:"):
        command = [shutil.which("conda") or "conda", "run", "--no-capture-output", "-n", conda_env, *command]
    if environment.startswith("module+conda:"):
        version = subprocess.run(["/bin/bash", "-c", _shell_environment_command(aligner, module_load, conda_env, f"{shlex.quote(aligner)} --version")], text=True, capture_output=True)
    else:
        version = subprocess.run(command, text=True, capture_output=True)
    if version.returncode != 0:
        command[-1] = "-V"
        if environment.startswith("module+conda:"):
            version = subprocess.run(["/bin/bash", "-c", _shell_environment_command(aligner, module_load, conda_env, f"{shlex.quote(aligner)} -V")], text=True, capture_output=True)
        else:
            version = subprocess.run(command, text=True, capture_output=True)
    if version.returncode != 0:
        raise RuntimeError(f"Unable to determine version for {executable!r}: {version.stderr.strip()}")
    return {"aligner": aligner, "resolved_executable": executable, "version": (version.stdout or version.stderr).strip().splitlines()[0], "threads": effective_thread_count(options, threads), "environment": environment, "status": "PASS"}
