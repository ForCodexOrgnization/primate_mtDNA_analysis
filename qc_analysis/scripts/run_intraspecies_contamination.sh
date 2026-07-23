#!/usr/bin/env bash
set -euo pipefail
CONFIG="${1:-config/qc_preprocessing.yaml}"; PYTHON="${PYTHON:-python3}"; RSCRIPT="${RSCRIPT:-Rscript}"
[[ -s "$CONFIG" ]] || { echo "ERROR: missing config: $CONFIG" >&2; exit 2; }
eval "$("$PYTHON" - "$CONFIG" <<'PY'
import sys
try: import yaml
except ImportError: raise SystemExit('PyYAML is required by this wrapper')
x=yaml.safe_load(open(sys.argv[1]))['intraspecies_contamination']
for k in ('enabled','build_variant_table','overwrite','vcf_dir','metadata','variant_table','negative_control_pairs','outdir','dp_min','use_snv_only'):
 v=x.get(k); print('%s=%r' % (k.upper(), v))
PY
)"
[[ "$ENABLED" == True ]] || { echo "ERROR: intraspecies_contamination.enabled must be true" >&2; exit 2; }
mkdir -p "$OUTDIR/input" "$OUTDIR/logs"
if [[ "$BUILD_VARIANT_TABLE" == True ]]; then
 [[ -n "$VCF_DIR" && -n "$METADATA" ]] || { echo "ERROR: build mode requires vcf_dir and metadata" >&2; exit 2; }
 VARIANT_TABLE="$OUTDIR/input/all_PASS_variants_core_table.tsv"
 cmd=("$PYTHON" qc_analysis/scripts/build_intraspecies_variant_table.py --vcf-dir "$VCF_DIR" --metadata "$METADATA" --output "$VARIANT_TABLE" --min-dp "$DP_MIN" --pass-only --log-file "$OUTDIR/input/variant_table_build_warnings.log")
 [[ "$USE_SNV_ONLY" == True ]] && cmd+=(--snv-only); [[ "$OVERWRITE" == True ]] && cmd+=(--overwrite)
 echo "[intraspecies] vcf_dir=$VCF_DIR metadata=$METADATA variant_table=$VARIANT_TABLE outdir=$OUTDIR"; "${cmd[@]}"
else
 [[ -n "$VARIANT_TABLE" ]] || { echo "ERROR: pre-built mode requires variant_table" >&2; exit 2; }
fi
cmd=("$RSCRIPT" qc_analysis/scripts/run_intraspecies_contamination.R --variant-table "$VARIANT_TABLE" --outdir "$OUTDIR" --config "$CONFIG"); [[ -n "$NEGATIVE_CONTROL_PAIRS" ]] && cmd+=(--negative-control-pairs "$NEGATIVE_CONTROL_PAIRS"); [[ "$OVERWRITE" == True ]] && cmd+=(--overwrite); "${cmd[@]}"
