#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
config=config/qc_preprocessing_coverage_hom_only_dp50_selective.yaml
runner=qc_analysis/scripts/run_hom_only_dp50_selective.py
pre=qc_analysis/scripts/run_pre_liftover_variant_qc.py
log=results/qc_coverage_only_dp50/slurm_logs
python3 "$runner" --phase prepare
mkdir -p "$log"
manifest=results/qc_coverage_only_dp50/pre_liftover_samples.txt
python3 - "$manifest" <<'PY'
import sys
from pathlib import Path
root=Path("results/qc_coverage_only/collected_variant_calling_results/collected_vcf")
samples=sorted(set(p.name.split(".round2",1)[0] for p in root.glob("*.round2.original_coords.clean.final.split.vcf*")))
if not samples: raise SystemExit("ERROR: no collected VCF samples")
Path(sys.argv[1]).write_text("\n".join(samples)+"\n")
print("Samples:", len(samples))
PY
n=$(wc -l < "$manifest")
array=$(sbatch --parsable --job-name=dp50_source --array="1-$n%20" --mem=8G --time=24:00:00 --output="$log/%A_%a.out" --error="$log/%A_%a.err" --wrap="set -euo pipefail; cd '$PWD'; sample=\$(sed -n \"\${SLURM_ARRAY_TASK_ID}p\" '$manifest'); python3 '$pre' --config '$config' --sample \"\$sample\"")
merge=$(sbatch --parsable --dependency="afterok:$array" --job-name=dp50_merge --mem=8G --time=24:00:00 --output="$log/%j.out" --error="$log/%j.err" --wrap="set -euo pipefail; cd '$PWD'; python3 '$pre' --config '$config' --merge-reports")
artifact=$(sbatch --parsable --dependency="afterok:$merge" --job-name=dp50_artifact --mem=64G --time=24:00:00 --output="$log/%j.out" --error="$log/%j.err" --wrap="set -euo pipefail; cd '$PWD'; bash qc_analysis/scripts/run_local_heteroplasmy_qc_with_expansion.sh '$config'")
refresh=$(sbatch --parsable --dependency="afterok:$artifact" --job-name=dp50_refresh --mem=96G --time=24:00:00 --output="$log/%j.out" --error="$log/%j.err" --wrap="set -euo pipefail; cd '$PWD'; python3 '$runner' --phase refresh")
final=$(sbatch --parsable --dependency="afterok:$refresh" --job-name=dp50_final --mem=48G --time=24:00:00 --output="$log/%j.out" --error="$log/%j.err" --wrap="set -euo pipefail; cd '$PWD'; python3 '$runner' --phase final")
hom=$(sbatch --parsable --dependency="afterok:$final" --job-name=dp50_hom --mem=16G --time=24:00:00 --output="$log/%j.out" --error="$log/%j.err" --wrap="set -euo pipefail; cd '$PWD'; python3 '$runner' --phase hom")
printf 'Submitted: source=%s merge=%s artifact=%s refresh=%s final=%s hom=%s\n' "$array" "$merge" "$artifact" "$refresh" "$final" "$hom"
