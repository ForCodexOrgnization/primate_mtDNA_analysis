import os

import pytest


@pytest.fixture(autouse=True)
def isolate_exported_lmod_functions(monkeypatch, tmp_path):
    """Prevent login-node Lmod shell functions from leaking into subprocess tests.

    Several wrapper tests intentionally place a fake `module` executable on PATH.
    Yale login shells can both export the real Lmod `module` function and recreate
    it during non-interactive Bash startup. Remove inherited exported functions,
    then force Bash to unset any recreated `module` function via BASH_ENV. This is
    test-only isolation; production wrapper behavior is unchanged.
    """
    for key in list(os.environ):
        if key.startswith("BASH_FUNC_module"):
            monkeypatch.delenv(key, raising=False)

    bash_env = tmp_path / "pytest_bash_env.sh"
    bash_env.write_text("unset -f module 2>/dev/null || true\n")
    monkeypatch.setenv("BASH_ENV", str(bash_env))
