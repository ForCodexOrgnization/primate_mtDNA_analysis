#!/usr/bin/env bash
# Stage an isolated incremental QC workspace from historical qc_back outputs.
# Usage: bash qc_analysis/scripts/prepare_incremental_strict.sh /path/to/updated_collection_summary.tsv
set -Eeuo pipefail

if [[ $# -ne 1 || ! -s "$1" ]]; then
  echo "Usage: bash qc_analysis/scripts/prepare_incremental_strict.sh /path/to/updated_collection_summary.tsv" >&2
  exit 2
fi
cd "$(dirname "$0")/../.."
summary_src="$(readlink -f "$1")"
backup="results/qc_back"
target="results/qc_incremental"
if [[ ! -d "$backup" ]]; then
  echo "ERROR: missing backup $backup" >&2
  exit 1
fi
if [[ -e "$target" ]]; then
  echo "ERROR: destination already exists: $target" >&2
  echo "For safety, do not reuse or overwrite a previous incremental workspace." >&2
  exit 1
fi
if ! command -v rsync >/dev/null 2>&1; then
  echo "ERROR: rsync is required" >&2
  exit 127
fi
mkdir -p "$target"
echo "Copying backup into isolated incremental workspace (may take a while)..."
rsync -a --exclude='job_arrays/' --exclude='logs/job_arrays/' "$backup/" "$target/"
cp "$summary_src" "$target/updated_collection_summary.tsv"
python3 - <<'PY'
import csv
from pathlib import Path
p=Path('results/qc_incremental/updated_collection_summary.tsv')
with p.open(newline='',encoding='utf-8') as h:
    reader=csv.DictReader(h,delimiter='\t')
    required={'sample','species','status','Percent_100','nuclear_median_coverage','mtcn_median'}
    if not reader.fieldnames or not required.issubset(reader.fieldnames):
        raise SystemExit(f'ERROR: updated summary missing columns: {sorted(required-set(reader.fieldnames or []))}')
    rows=list(reader)
names=[row['sample'] for row in rows]
if not rows or len(names)!=len(set(names)):
    raise SystemExit('ERROR: empty summary or duplicated sample identifiers')
print(f'[incremental] updated summary samples={len(rows)}')
PY
echo "[incremental] staged: $target"
echo "[incremental] backup remains untouched: $backup"
echo "Next: inspect and submit the incremental workflow only after confirming compatibility."
