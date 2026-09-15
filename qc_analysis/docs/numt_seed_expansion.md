# Limited NUMT seed expansion

This step follows `run_local_heteroplasmy_qc.py` and performs a conservative expansion around heteroplasmy clusters that already meet a REMOVE-level NUMT rule.

## Strong seeds allowed to expand

A strict cluster can drive expansion only when it is already marked `REMOVE` because of either:

1. direct sample-level NUMT support; or
2. species-level NUMT support plus same-species recurrence.

Species-level NUMT support alone (`FLAG`) never drives removal expansion.

## Expansion rule

Strict seed discovery remains unchanged:

- >=3 HET variants
- circular 250-bp local window
- strict seed/cluster AF span <=0.06

For a residual HET to be added to a strong NUMT-supported seed, all of the following must hold:

- same sample;
- assigned to the closest appropriate strong seed;
- inside the same NUMT -> chrM interval as that seed;
- distance from the strict seed cluster <=250 bp;
- `|residual AF - original seed median AF| <=0.07`;
- after adding variants greedily from smallest AF distance outward, total expanded AF span remains <=0.10.

The original strict-seed median is never updated during expansion, preventing chain drift.

## Why these defaults

Cohort sensitivity analysis showed that increasing the AF expansion threshold from 0.06 to 0.07 recovered additional NUMT-supported residual HETs without increasing the nearby-outside-NUMT control count, while thresholds above 0.07 showed diminishing returns. The 0.10 total-span cap prevents broad AF bands from being absorbed into a NUMT seed.

The current implementation therefore uses defaults:

- strict cluster AF span: `0.06`
- expansion distance: `250 bp`
- maximum AF distance from original seed median: `0.07`
- maximum total expanded AF span: `0.10`

These are configurable through optional keys under `local_heteroplasmy_qc:`:

```yaml
numt_seed_expansion_enabled: true
numt_seed_expansion_distance_bp: 250
numt_seed_expansion_max_delta_af: 0.07
numt_seed_expansion_max_total_af_span: 0.10
```

If these keys are absent, the values above are used automatically.

## Run

Run strict clustering plus limited expansion with:

```bash
bash qc_analysis/scripts/run_local_heteroplasmy_qc_with_expansion.sh \
  config/qc_preprocessing.yaml
```

Or run the expansion alone after a fresh strict-cluster run:

```bash
python3 qc_analysis/scripts/expand_numt_seed_variants.py \
  --config config/qc_preprocessing.yaml
```

## Added outputs / annotations

The expansion step updates the standard heteroplasmy reports and `numt_variants_to_remove.tsv`, and also writes:

- `numt_seed_expansion_variants.tsv`
- `numt_seed_expansion_summary.tsv`

Expanded variants are annotated with:

- `numt_expanded=YES`
- `numt_expansion_seed_id`
- `numt_expansion_delta_af`
- `numt_expansion_distance_bp`
- `numt_expansion_final_af_span`
- `filter_reason=NUMT_SEED_EXPANSION`

The step is rerunnable: previously expanded variants are restored to their pre-expansion state before the expansion is recalculated.
