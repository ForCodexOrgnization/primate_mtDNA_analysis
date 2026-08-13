# QC analysis pipeline

本目录包含变异检测结果收集、全局锚点发现、坐标转换以及密码子、tRNA、rRNA
注释等 QC 步骤。推荐统一通过
[`scripts/run_qc_preprocessing.sh`](scripts/run_qc_preprocessing.sh) 运行，而不是直接调用各个
Python/R 脚本。默认配置文件是
[`../config/qc_preprocessing.yaml`](../config/qc_preprocessing.yaml)。运行前请先检查其中的输入、
输出、参考序列和软件环境路径。

## 快速开始（推荐：Slurm）

在仓库根目录、集群登录节点运行：

```bash
bash qc_analysis/scripts/run_qc_preprocessing.sh --submit all config/qc_preprocessing.yaml
```

`--submit all` 会把各步骤作为独立 Slurm job/array 提交，并用 `afterok` 依赖保证下游只在
上游成功后启动。这是完整数据集的推荐运行方式；不要在计算节点或已有 Slurm job 中再次
使用 `--submit`。

如需在当前 shell 中按顺序串行执行（适合调试或小数据集）：

```bash
bash qc_analysis/scripts/run_qc_preprocessing.sh all config/qc_preprocessing.yaml
```

直接模式不会提交 job array，所有样本会由当前进程依次处理。MITOS2 等耗时步骤在生产数据
上通常应使用 `--submit`。

## `all` 的运行顺序

完整主流程按以下依赖顺序运行：

| 顺序 | wrapper step | 作用 | 主要下游依赖 |
|---:|---|---|---|
| 1 | `collect_variant_calling_results` | 收集并标准化 variant-calling 的 VCF、coverage 和 mtCN 结果 | 坐标转换输入 |
| 2 | `intraspecies_contamination` | 在原始物种坐标上生成样本级种内污染报告（不修改 VCF） | 最终样本判定 |
| 3 | `discover_global_anchor` | 对参考线粒体序列做全局比对并生成经过验证的 anchor 表 | 坐标转换 anchor |
| 4 | `coordinate_liftover` | 将每个样本的原始物种坐标转换到人 chrM 坐标 | 后续 VCF 注释 |
| 5 | `mitos2_prepare_tasks` | 生成每个目标参考序列一条记录的 MITOS2 任务表 | 仅 `--submit all` 中的显式节点 |
| 6 | `mitos2_annotation` | 按参考序列运行 MITOS2 | MITOS2 合并结果 |
| 7 | `mitos2_merge` | 合并并严格质控 MITOS2 结果，生成生产密码子表和样本映射 | `codon_match_validate` |
| 8 | `compare_genbank_mitos2` | 可选：按相同序列 SHA256 比较独立 GenBank 与 MITOS2 注释 | 验证证据（不阻塞生产流程） |
| 9 | `codon_match_validate` | 在读取 VCF 前校验并建立密码子输入索引 | 仅 `--submit all` 中的显式校验节点 |
| 10 | `codon_match` | 给 lifted VCF 添加密码子匹配注释 | tRNA 注释输入 |
| 11 | `codon_match_merge` | 原子合并每个样本的密码子汇总 | cohort 汇总 |

MITOS2 annotations of the exact variant-calling FASTAs are the only production
CDS/codon source. GenBank remains an independent sequence-hash-matched benchmark;
it never overrides or substitutes for MITOS2 in production.
MITOS2 input and coordinate FASTA sequence identity must be established; accession or
biological similarity alone is insufficient.

Codon matching retains overlapping CDS genes and evaluates all compatible source ×
human gene/phase candidate pairs instead of collapsing a position to one gene.
| 12 | `trna_match` | 给 VCF 添加 tRNA 匹配和结构相关注释 | rRNA 注释输入 |
| 13 | `rrna_match` | 给 VCF 添加 rRNA 区域/可选结构注释 | 最终注释报告 |
| 14 | `final_filter` | 汇总所有样本级和变异级报告并一次性生成最终文件 | `final_vcf/final_cov/final_mtcn` |

上表是 `--submit all` 创建的完整依赖图。直接运行 `all` 时，wrapper 会在单一进程中完成相同
的主要生物学步骤，但任务准备、输入校验和部分合并操作会由相应步骤内部处理，而不是作为
独立 Slurm 节点显示。

`collect_variant_calling_results` 的输入是 run-manager 聚合根目录，其中包含
`vcf/`、`mtcn/`、`round2_coverage/` 和 `numt_decoy_coverage/`（`receipts/` 会被忽略）。
样本由 metadata 表提供；两个 coverage 文件以 `(chrom, pos, target)` 为键逐位取最大深度，
并写入稳定的 `collected_cov/{sample}.merged.max_depth.per_base_coverage.tsv` 下游接口。

`intraspecies_contamination` 现在紧跟收集步骤并属于 `all`；它只写样本级报告。
`final_filter` 是唯一执行最终排除的终端步骤。详见
[`docs/intraspecies_contamination.md`](docs/intraspecies_contamination.md)。

## 分步骤运行

如果需要检查中间结果或重跑失败阶段，应保持上表顺序。例如：

```bash
CONFIG=config/qc_preprocessing.yaml

bash qc_analysis/scripts/run_qc_preprocessing.sh --submit collect_variant_calling_results "$CONFIG"
bash qc_analysis/scripts/run_qc_preprocessing.sh --submit intraspecies_contamination "$CONFIG"
bash qc_analysis/scripts/run_qc_preprocessing.sh --submit discover_global_anchor "$CONFIG"
bash qc_analysis/scripts/run_qc_preprocessing.sh --submit coordinate_liftover "$CONFIG"
bash qc_analysis/scripts/run_qc_preprocessing.sh --submit mitos2_annotation "$CONFIG"
bash qc_analysis/scripts/run_qc_preprocessing.sh --submit build_primate_codon_table "$CONFIG"
bash qc_analysis/scripts/run_qc_preprocessing.sh --submit codon_match_validate "$CONFIG"
bash qc_analysis/scripts/run_qc_preprocessing.sh --submit codon_match "$CONFIG"
bash qc_analysis/scripts/run_qc_preprocessing.sh --submit trna_match "$CONFIG"
bash qc_analysis/scripts/run_qc_preprocessing.sh --submit rrna_match "$CONFIG"
bash qc_analysis/scripts/run_qc_preprocessing.sh --submit final_filter "$CONFIG"
```

注意：单独以 `--submit` 运行 `mitos2_annotation` 或 `codon_match` 时，wrapper 默认分别自动提交
`mitos2_merge` 或 `codon_match_merge`，并设置 `afterok` 依赖。通常不需要再次手动提交 merge。
如需自行控制依赖，可在提交 producer 时设置 `AUTO_SUBMIT_MERGE=false`。

运行独立的种内污染检查：

```bash
bash qc_analysis/scripts/run_qc_preprocessing.sh --submit intraspecies_contamination config/qc_preprocessing.yaml
```

最终目录为 `results/qc/final_filter/`，其中 `reports/final_sample_qc.tsv`
汇总四类样本状态和最终原因。`reports/final_variant_qc.tsv` 使用
`sample,human_chrom,human_pos,human_ref,human_alt` 作为 post-liftover canonical join key。
下游 VCF 坐标只写入 `human_*`；无法可靠恢复时 `source_*` 以及兼容保留但已弃用的
`original_*` 写为 `NOT_AVAILABLE`，不会冒充物种原始坐标。变异报告中的 generic
`CHROM/POS/REF/ALT` 必须显式配置 `coordinate_system: human`，
`reports/final_filter_summary.tsv` 提供计数；仅样本 PASS 且变异 PASS 的记录进入
`final_vcf/`，其 coverage 与 mtCN 文件分别复制到 `final_cov/` 和 `final_mtcn/`。

## 单样本运行与 array 并发

`coordinate_liftover`、`codon_match`、`trna_match` 和 `rrna_match` 是按样本拆分的 array；
`mitos2_annotation` 按参考序列拆分。默认最多同时运行 20 个 array task：

```bash
# 修改整个提交的 array 并发上限
bash qc_analysis/scripts/run_qc_preprocessing.sh \
  --submit --array-concurrency 40 codon_match config/qc_preprocessing.yaml

# 只运行一个样本
bash qc_analysis/scripts/run_qc_preprocessing.sh \
  --submit --sample SAMPLE_NAME coordinate_liftover config/qc_preprocessing.yaml

# SAMPLE 环境变量写法等价
SAMPLE=SAMPLE_NAME \
  bash qc_analysis/scripts/run_qc_preprocessing.sh \
  --submit codon_match config/qc_preprocessing.yaml
```

已完成且通过完整性检查的 array 项默认会被跳过。需要强制纳入任务表时使用：

```bash
FORCE_RERUN=true \
  bash qc_analysis/scripts/run_qc_preprocessing.sh \
  --submit codon_match config/qc_preprocessing.yaml
```

只生成任务清单、不提交，或预览最终 `sbatch` 命令：

```bash
bash qc_analysis/scripts/run_qc_preprocessing.sh \
  --submit --prepare-only codon_match config/qc_preprocessing.yaml

bash qc_analysis/scripts/run_qc_preprocessing.sh \
  --dry-run-submit codon_match config/qc_preprocessing.yaml
```

任务清单、submission metadata 和日志默认写在各步骤 output 目录下的 `job_arrays/` 和
`logs/job_arrays/`。重试、资源覆盖和目录解析规则详见
[`docs/slurm_job_arrays.md`](docs/slurm_job_arrays.md)。

## 常用资源和环境覆盖

```bash
# 通用 Slurm 资源
SLURM_TIME=48:00:00 SLURM_MEM=32G SLURM_CPUS=8 \
  bash qc_analysis/scripts/run_qc_preprocessing.sh \
  --submit build_primate_codon_table config/qc_preprocessing.yaml

# 步骤专用资源（优先于通用值）
LIFTOVER_SLURM_TIME=12:00:00 LIFTOVER_SLURM_CPUS=8 \
  bash qc_analysis/scripts/run_qc_preprocessing.sh \
  --submit coordinate_liftover config/qc_preprocessing.yaml

# codon table 内部 worker 数应与申请的 CPU 数相符
SLURM_CPUS=8 CODON_TABLE_WORKERS=8 \
  bash qc_analysis/scripts/run_qc_preprocessing.sh \
  --submit build_primate_codon_table config/qc_preprocessing.yaml
```

可用的步骤专用前缀包括 `LIFTOVER_`、`CODON_MATCH_`、`TRNA_MATCH_`、`RRNA_MATCH_`
和 `MITOS2_`，后接 `SLURM_TIME`、`SLURM_MEM` 或 `SLURM_CPUS`。还可使用
`SLURM_PARTITION` 指定分区。Biopython 与 MITOS2 的 module/conda 设置来自配置文件，必要时
可用 wrapper `--help` 中列出的环境变量覆盖。

## 查看帮助与详细文档

```bash
bash qc_analysis/scripts/run_qc_preprocessing.sh --help
```

各步骤的输入、输出和算法说明位于 [`docs/`](docs/)；遇到失败时，先查看相应步骤的
`logs/job_arrays/*.err`、任务 manifest 以及 YAML 中的路径和 `enabled` 设置。

## Terminal sample QC and filtering

`intraspecies_contamination` runs in each species' **original coordinates**. The
validation truth is `validation/contamination_reference.R`; the production job
is the dependency-light Python implementation. Mirror evidence is the complete
within-sample low-AF × high-AF cross product and Tier-2 negative controls
calibrate its normalized p95/p99 evidence.

`sample_variant_filtering` is report-only. A sample passes when mt median depth
is at least 100, at least 90% of mtDNA bases have 100× depth, nuclear median
depth is at least 20, median mtCN is at least 40, and MAD is strictly below
0.5. Neither this step nor contamination analysis removes data.

`final_filter` is the only irreversible step. Intraspecies and sample-QC reports
are required by default; human and interspecies reports are optional. Passing
samples require sample QC PASS and must not be `high_confidence_contaminated`.
Candidate/insufficient contamination results remain warnings. VCFs are selected
in deterministic priority order: rRNA, tRNA, codon, then raw coordinate
liftover. There is intentionally no original-coordinate collection fallback.
Final VCFs are BGZF compressed and tabix indexed (pysam or bgzip+tabix required).
# Human mtDNA contamination QC

`run_human_contamination.py` runs immediately after coordinate liftover. It
matches only the canonical post-liftover Human rCRS `POS + ALT` allele; source
coordinates and source alleles are never used for PhyloTree matching. This is
important when liftover has performed an `ALT_REF_FLIP` and transformed the
allele-specific FORMAT values.

The Human-contamination cohort is the set of successfully emitted lifted VCFs
in `input_vcf_dir`. `sample_ref_file.tsv` supplies species labels only: rows
without a lifted VCF are reported by input validation as out-of-scope, not as
missing inputs, and a lifted sample absent from metadata is retained with an
empty species label.

The transparent baseline screen asks whether sufficiently deep, low-VAF
primate variants are enriched for Human PhyloTree SNVs. The default low-VAF
range is 0.01--0.50, with at least six denominator variants, six distinct
marker hits, and a marker fraction of at least 0.60. A strict `FAIL` additionally
requires at least 70% of marker AFs to lie within 0.03 of their median and at
least three eligible non-control-region hits. A baseline-positive sample that
lacks strict support is a `CANDIDATE`; fewer than six low-VAF variants is
`INSUFFICIENT_DATA`, not `PASS`.

HaploGrep 3 is optional supplementary phylogenetic characterization. It is
given a synthetic profile containing only selected, deep, low-VAF Human-marker
alleles—never the complete primate VCF. An explicitly configured wrapper
executable takes precedence over `jar + java`. Missing tools are reported as
`TOOL_UNAVAILABLE` unless strict tool availability is requested. HaploGrep
quality is not a probability and is not required for `FAIL`: low quality can
reflect a sparse contaminant profile with too few downstream markers for a
fine-scale assignment.

Production screening requires the complete PhyloTree/HaploGrep rCRS v17.1
marker export, not the small illustrative table formerly shipped here. Prepare
the normalized SNV table reproducibly (complex mutations are counted and
excluded, POS+ALT duplicates are merged) with:

```bash
python3 qc_analysis/scripts/prepare_human_phylotree_markers.py \
  --input data/reference_tables/haplogrep-rcrs-v17.1_uniq_SNP.txt \
  --output data/reference_tables/human_phylotree_rcrs_v17.1_snv.tsv \
  --summary data/reference_tables/human_phylotree_rcrs_v17.1_snv.qc.json
```

The minimum marker count guard prevents an accidental placeholder from being
used; `allow_small_marker_reference` is intended only for synthetic tests.
HaploGrep receives one headerless HSD row whose selected mutations are separate
tab fields. Such a contaminant-only profile may be sparse. Its quality is
supplementary and is not a probability. The recommended default
`quality_required_for_fail: false` leaves marker enrichment, VAF coherence, and
non-control-region support solely responsible for `FAIL`. When enabled, absent,
failed, low-quality, or insufficient phylogenetic evidence changes an otherwise
strict `FAIL` to `CANDIDATE`, rather than silently passing it. HaploGrep remains
optional unless `require_tool_when_enabled: true`. Keep `marker_version` and the
configured HaploGrep tree synchronized (currently rCRS v17.1).

Validate configuration and inputs without analysis:

```bash
python3 qc_analysis/scripts/run_human_contamination.py \
  --config config/qc_preprocessing.yaml --validate-inputs
```

Run directly or submit the singleton step:

```bash
bash qc_analysis/scripts/run_qc_preprocessing.sh human_contamination config/qc_preprocessing.yaml
bash qc_analysis/scripts/run_qc_preprocessing.sh --submit human_contamination config/qc_preprocessing.yaml
```

## Non-destructive orthology and primate background workflow

The post-liftover dependency order is:

```
coordinate_liftover
  -> MITOS2/reference preparation
  -> codon_match -> trna_match -> rrna_match
  -> build_primate_homo_background
  -> human_contamination
  -> interspecies contamination report (when configured)
  -> final_filter
```

Codon/tRNA/rRNA matching generates QC annotations and reports; variants are not
removed until `final_filter`. The fully annotated, unfiltered input is
`results/qc/rrna_match/vcf_rrna/{sample}.lifted.codon.trna.rrna.vcf`.
`build_primate_homo_background` consolidates those existing annotations into
`orthology_match_report.tsv`; it does not replace or rewrite the VCF.

The primate homoplasmic background is an intermediate calibration dataset used
to distinguish recurrent/evolutionary primate alleles from more informative
Human PhyloTree marker hits. By default it contains exact post-liftover
`human_pos + human_ref + human_alt` SNVs with VCF `FILTER=PASS`, DP >= 100,
AF >= 0.95, and orthology status `PASS`. It does not read Human or interspecies
contamination status and is never a final-output VCF source. Human QC retains
its historical all-marker screen and reports background and informative subsets;
the new corrected metrics do not affect classification by default. HaploGrep
remains supplementary (`quality_required_for_fail: false`).
