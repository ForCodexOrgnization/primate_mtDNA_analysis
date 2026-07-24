import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).parents[2]
SCRIPT = ROOT / "qc_analysis/scripts/run_intraspecies_contamination.sh"
YAML_PYTHON = "/usr/bin/python3"


def run_wrapper(config, *, cwd=ROOT, env=None):
    environment = os.environ.copy()
    environment["PYTHON"] = YAML_PYTHON
    if env:
        environment.update(env)
    return subprocess.run(
        ["bash", str(SCRIPT), str(config)], cwd=cwd, env=environment,
        text=True, capture_output=True,
    )


def write_config(path, section):
    path.write_text("intraspecies_contamination:\n" + section, encoding="utf-8")


def test_execution_does_not_require_pyyaml(tmp_path):
    config = tmp_path / "config.yaml"
    write_config(config, "  enabled: false\n  outdir: output\n")
    result = run_wrapper(config)
    assert result.returncode == 0, result.stderr
    assert "PyYAML" not in result.stderr


def test_malformed_yaml_exits_nonzero(tmp_path):
    config = tmp_path / "bad.yaml"
    config.write_text("intraspecies_contamination: [\n", encoding="utf-8")
    result = run_wrapper(config)
    assert result.returncode != 0
    assert "YAML mapping" in result.stderr
    assert "unbound variable" not in result.stderr


def test_missing_section_exits_nonzero(tmp_path):
    config = tmp_path / "missing.yaml"
    config.write_text("other: true\n", encoding="utf-8")
    result = run_wrapper(config)
    assert result.returncode != 0
    assert "missing 'intraspecies_contamination'" in result.stderr


def test_disabled_step_normalizes_values_and_nulls(tmp_path):
    config = tmp_path / "disabled config.yaml"
    write_config(config, """  enabled: false
  build_variant_table: true
  overwrite: false
  use_snv_only: true
  vcf_dir:
  metadata:
  variant_table:
  negative_control_pairs:
  outdir: "output directory"
""")
    result = run_wrapper(config)
    assert result.returncode == 0, result.stderr
    assert "[intraspecies] disabled; skipping." in result.stdout
    assert "enabled=false" in result.stdout
    assert "build_variant_table=true" in result.stdout
    assert "vcf_dir=<not set>" in result.stdout
    assert "outdir=output directory" in result.stdout


def test_paths_with_spaces_are_shell_safe(tmp_path):
    config = tmp_path / "paths.yaml"
    table = tmp_path / "table with spaces.tsv"
    outdir = tmp_path / "output with spaces"
    write_config(config, f"""  enabled: true
  build_variant_table: false
  variant_table: "{table}"
  outdir: "{outdir}"
""")
    rscript = tmp_path / "fake Rscript"
    rscript.write_text("#!/usr/bin/env bash\nprintf '%s\\n' \"$@\" > \"$CAPTURE\"\n", encoding="utf-8")
    rscript.chmod(0o755)
    capture = tmp_path / "arguments.txt"
    result = run_wrapper(config, env={"RSCRIPT": str(rscript), "CAPTURE": str(capture)})
    assert result.returncode == 0, result.stderr
    arguments = capture.read_text(encoding="utf-8").splitlines()
    assert str(table) in arguments
    assert str(outdir) in arguments


def test_default_config_works_outside_repository_root(tmp_path):
    result = subprocess.run(
        ["bash", str(SCRIPT)], cwd=tmp_path, text=True, capture_output=True,
        env={**os.environ, "PYTHON": YAML_PYTHON},
    )
    assert result.returncode == 0, result.stderr
    assert f"config={ROOT / 'config/qc_preprocessing.yaml'}" in result.stdout
    assert "[intraspecies] disabled; skipping." in result.stdout
