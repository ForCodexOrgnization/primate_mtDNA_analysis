#!/usr/bin/env bash
set -euo pipefail

if ! command -v quarto >/dev/null 2>&1; then
  echo "ERROR: quarto is not available on PATH. Install/load Quarto before rendering reports." >&2
  exit 1
fi

shopt -s nullglob
reports=(preprocessing/analysis/*.qmd)
if [[ ${#reports[@]} -eq 0 ]]; then
  echo "ERROR: no preprocessing Quarto reports found under preprocessing/analysis/*.qmd" >&2
  exit 1
fi

for report in "${reports[@]}"; do
  echo "[preprocessing] Rendering ${report}" >&2
  quarto render "$report"
done
