#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
CONFIG="${1:-$REPO_ROOT/config/qc_preprocessing.yaml}"
PYTHON="${PYTHON:-python3}"
RSCRIPT="${RSCRIPT:-Rscript}"

# Defaults keep diagnostics safe under `set -u`; parser failures still exit below.
ENABLED=""
BUILD_VARIANT_TABLE=""
OVERWRITE="false"
VCF_DIR=""
METADATA=""
VARIANT_TABLE=""
NEGATIVE_CONTROL_PAIRS=""
OUTDIR=""
DP_MIN="100"
USE_SNV_ONLY="true"

[[ -s "$CONFIG" ]] || { echo "ERROR: missing config: $CONFIG" >&2; exit 2; }

if ! "$PYTHON" -c 'import yaml' >/dev/null 2>&1; then
    echo "ERROR: PyYAML is required by run_intraspecies_contamination.sh." >&2
    echo "Install it with one of:" >&2
    echo "  python3 -m pip install --user PyYAML" >&2
    echo "  conda install -c conda-forge pyyaml" >&2
    exit 2
fi

if ! CONFIG_ASSIGNMENTS="$("$PYTHON" - "$CONFIG" <<'PY'
import shlex
import sys

import yaml


def emit(name, value):
    if value is None:
        value = ""
    elif isinstance(value, bool):
        value = "true" if value else "false"
    else:
        value = str(value)
    print(f"{name}={shlex.quote(value)}")


def boolean(section, name, default):
    value = section.get(name, default)
    if isinstance(value, bool):
        return value
    if isinstance(value, str) and value.lower() in {"true", "false"}:
        return value.lower() == "true"
    raise SystemExit(f"ERROR: intraspecies_contamination.{name} must be true or false")


def text_value(section, name, default=None):
    value = section.get(name, default)
    if value is None:
        return None
    if not isinstance(value, (str, int, float)) or isinstance(value, bool):
        raise SystemExit(f"ERROR: intraspecies_contamination.{name} must be a scalar or null")
    return str(value)


config_path = sys.argv[1]
try:
    with open(config_path, encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
except (OSError, yaml.YAMLError) as error:
    raise SystemExit(f"ERROR: failed to parse configuration: {error}")

if not isinstance(config, dict):
    raise SystemExit("ERROR: YAML root must be a mapping")

section = config.get("intraspecies_contamination")
if section is None:
    raise SystemExit("ERROR: missing 'intraspecies_contamination' section in configuration")
if not isinstance(section, dict):
    raise SystemExit("ERROR: 'intraspecies_contamination' must be a YAML mapping")

values = {
    "ENABLED": boolean(section, "enabled", False),
    "BUILD_VARIANT_TABLE": boolean(section, "build_variant_table", True),
    "OVERWRITE": boolean(section, "overwrite", False),
    "VCF_DIR": text_value(section, "vcf_dir"),
    "METADATA": text_value(section, "metadata"),
    "VARIANT_TABLE": text_value(section, "variant_table"),
    "NEGATIVE_CONTROL_PAIRS": text_value(section, "negative_control_pairs"),
    "OUTDIR": text_value(section, "outdir"),
    "DP_MIN": text_value(section, "dp_min", 100),
    "USE_SNV_ONLY": boolean(section, "use_snv_only", True),
}

if not values["OUTDIR"]:
    raise SystemExit("ERROR: intraspecies_contamination.outdir must be nonempty")
try:
    if int(values["DP_MIN"]) < 0:
        raise ValueError
except (TypeError, ValueError):
    raise SystemExit("ERROR: intraspecies_contamination.dp_min must be a non-negative integer")

for name, value in values.items():
    emit(name, value)
PY
)"; then
    echo "ERROR: failed to parse configuration: $CONFIG" >&2
    exit 2
fi

eval "$CONFIG_ASSIGNMENTS"

echo "[intraspecies] config=$CONFIG"
echo "[intraspecies] enabled=$ENABLED"
echo "[intraspecies] build_variant_table=$BUILD_VARIANT_TABLE"
echo "[intraspecies] vcf_dir=${VCF_DIR:-<not set>}"
echo "[intraspecies] metadata=${METADATA:-<not set>}"
echo "[intraspecies] variant_table=${VARIANT_TABLE:-<not set>}"
echo "[intraspecies] negative_control_pairs=${NEGATIVE_CONTROL_PAIRS:-<not set>}"
echo "[intraspecies] outdir=$OUTDIR"

if [[ "$ENABLED" != "true" ]]; then
    echo "[intraspecies] disabled; skipping."
    exit 0
fi

mkdir -p "$OUTDIR/input" "$OUTDIR/logs"
if [[ "$BUILD_VARIANT_TABLE" == "true" ]]; then
    [[ -n "$VCF_DIR" && -n "$METADATA" ]] || { echo "ERROR: build mode requires vcf_dir and metadata" >&2; exit 2; }
    VARIANT_TABLE="$OUTDIR/input/all_PASS_variants_core_table.tsv"
    cmd=("$PYTHON" "$REPO_ROOT/qc_analysis/scripts/build_intraspecies_variant_table.py" --vcf-dir "$VCF_DIR" --metadata "$METADATA" --output "$VARIANT_TABLE" --min-dp "$DP_MIN" --pass-only --log-file "$OUTDIR/input/variant_table_build_warnings.log")
    [[ "$USE_SNV_ONLY" == "true" ]] && cmd+=(--snv-only)
    [[ "$OVERWRITE" == "true" ]] && cmd+=(--overwrite)
    "${cmd[@]}"
else
    [[ -n "$VARIANT_TABLE" ]] || { echo "ERROR: pre-built mode requires variant_table" >&2; exit 2; }
fi
cmd=("$RSCRIPT" "$REPO_ROOT/qc_analysis/scripts/run_intraspecies_contamination.R" --variant-table "$VARIANT_TABLE" --outdir "$OUTDIR" --config "$CONFIG")
[[ -n "$NEGATIVE_CONTROL_PAIRS" ]] && cmd+=(--negative-control-pairs "$NEGATIVE_CONTROL_PAIRS")
[[ "$OVERWRITE" == "true" ]] && cmd+=(--overwrite)
"${cmd[@]}"
