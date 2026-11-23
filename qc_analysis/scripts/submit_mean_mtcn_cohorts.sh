#!/usr/bin/env bash
# Submit either independent cohort without rerunning WGS alignment/liftover/annotation.
set -euo pipefail
cd "$(dirname "$0")/../.."

usage() {
  echo "Usage: bash qc_analysis/scripts/submit_mean_mtcn_cohorts.sh hom10|homhet40" >&2
  exit 2
}
[[ $# -eq 1 ]] || usage
cohort="$1"
case "$cohort" in
  hom10)
    config="config/qc_cohort_hom10_mean_mtcn.yaml"
    root="results/qc_cohort_hom10_mean_mtcn"
    ;;
  homhet40)
    config="config/qc_cohort_homhet40_mean_mtcn.yaml"
    root="results/qc_cohort_homhet40_mean_mtcn"
    ;;
  *) usage ;;
esac

python3 qc_analysis/scripts/prepare_mean_mtcn_cohorts.py --cohort "$cohort"
mkdir -p "$root/slurm_logs"
logs="$root/slurm_logs"
sample_script="qc_analysis/scripts/run_sample_variant_filtering.py"
source_script="qc_analysis/scripts/run_pre_liftover_variant_qc.py"
final_script="qc_analysis/scripts/run_final_filter_with_heteroplasmy.py"
hom_script="qc_analysis/scripts/run_homoplasmy_qc.py"
artifact_script="qc_analysis/scripts/run_local_heteroplasmy_qc_with_expansion.sh"

sample_job=$(sbatch --parsable \
  --job-name="${cohort}_sample" --mem=16G --time=24:00:00 \
  --output="$logs/%j.out" --error="$logs/%j.err" \
  --wrap="set -euo pipefail; cd '$PWD'; python3 '$sample_script' --config '$config'")
echo "Submitted sample QC: $sample_job"

upstream="$sample_job"
if [[ "$cohort" == "homhet40" ]]; then
  manifest="$root/source_samples.txt"
  python3 - "$manifest" <<'PY'
import sys
from pathlib import Path
directory = Path("results/qc_coverage_only/collected_variant_calling_results/collected_vcf")
samples = sorted({p.name.split(".round2", 1)[0] for p in directory.glob("*.round2.original_coords.clean.final.split.vcf*")})
if not samples:
    raise SystemExit("No source calling VCFs found")
Path(sys.argv[1]).write_text("\n".join(samples) + "\n")
print("Original samples for DP100 source QC:", len(samples))
PY
  n=$(wc -l < "$manifest")
  array_job=$(sbatch --parsable --dependency="afterok:$sample_job" \
    --job-name="${cohort}_source" --array="1-${n}%20" --mem=6G --time=24:00:00 \
    --output="$logs/%A_%a.out" --error="$logs/%A_%a.err" \
    --wrap="set -euo pipefail; cd '$PWD'; sample=\$(sed -n \"\${SLURM_ARRAY_TASK_ID}p\" '$manifest'); test -n \"\$sample\"; python3 '$source_script' --config '$config' --sample \"\$sample\"")
  echo "Submitted DP100 source QC array: $array_job"
  merge_job=$(sbatch --parsable --dependency="afterok:$array_job" \
    --job-name="${cohort}_merge" --mem=12G --time=12:00:00 \
    --output="$logs/%j.out" --error="$logs/%j.err" \
    --wrap="set -euo pipefail; cd '$PWD'; python3 '$source_script' --config '$config' --merge-reports")
  echo "Submitted source report merge: $merge_job"
  upstream="$merge_job"
fi

artifact_job=$(sbatch --parsable --dependency="afterok:$upstream" \
  --job-name="${cohort}_artifact" --mem=64G --time=24:00:00 \
  --output="$logs/%j.out" --error="$logs/%j.err" \
  --wrap="set -euo pipefail; cd '$PWD'; bash '$artifact_script' '$config'")
echo "Submitted local artifact QC: $artifact_job"

final_job=$(sbatch --parsable --dependency="afterok:$artifact_job" \
  --job-name="${cohort}_final" --mem=64G --time=24:00:00 \
  --output="$logs/%j.out" --error="$logs/%j.err" \
  --wrap="set -euo pipefail; cd '$PWD'; python3 '$final_script' --config '$config'")
echo "Submitted Final QC: $final_job"

if [[ "$cohort" == "hom10" ]]; then
  hom_job=$(sbatch --parsable --dependency="afterok:$final_job" \
    --job-name="${cohort}_hom" --mem=16G --time=12:00:00 \
    --output="$logs/%j.out" --error="$logs/%j.err" \
    --wrap="set -euo pipefail; cd '$PWD'; python3 '$hom_script' --config '$config'")
  echo "Submitted HOM adjusted report: $hom_job"
fi

echo "Submitted cohort=$cohort. Inspect: squeue -u \$USER"
echo "QC files in: $root"
