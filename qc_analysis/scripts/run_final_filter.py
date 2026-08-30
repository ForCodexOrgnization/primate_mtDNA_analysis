#!/usr/bin/env python3
"""Compatibility entrypoint for the streaming final-filter implementation."""
from __future__ import annotations

import sqlite3
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from qc_analysis.scripts.run_final_filter_streaming import *  # noqa: F401,F403,E402


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, KeyError, OSError, RuntimeError, sqlite3.Error, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
