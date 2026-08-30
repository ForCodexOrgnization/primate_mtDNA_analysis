import os

import pytest


@pytest.fixture(autouse=True)
def isolate_exported_lmod_functions(monkeypatch):
    """Prevent login-node Lmod shell functions from leaking into subprocess tests.

    Several wrapper tests intentionally place a fake `module` executable on PATH.
    Yale login shells export the real Lmod `module` function, and Bash gives an
    imported function precedence over PATH, causing those tests to hit the live
    module system. Remove only that inherited function from the pytest process;
    production wrapper behavior is unchanged.
    """
    for key in list(os.environ):
        if key.startswith("BASH_FUNC_module"):
            monkeypatch.delenv(key, raising=False)
