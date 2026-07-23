# Build primate codon table

`build_primate_codon_table.py` downloads (or reuses) GenBank records and writes a
reference-level coding-position table plus a sample-to-reference map used by
`codon_match`:

```bash
python qc_analysis/scripts/build_primate_codon_table.py --config config/qc_preprocessing.yaml
bash qc_analysis/scripts/run_qc_preprocessing.sh build_primate_codon_table config/qc_preprocessing.yaml
```

`codon_table_level: reference` is the default. In this mode the builder never
loads `mitos2_cds_table`, which is a potentially multi-gigabyte sample-level
artifact. MITOS2 fallback annotations come only from
`mitos2_reference_cds_table`, are grouped once by coordinate reference, and
are written once per selected reference. `sample_reference_map.tsv` contains
the per-sample association. Set `codon_table_level: sample` only for the
explicit legacy workflow.

On the HPC, the preprocessing wrapper loads the Biopython module only for the
`build_primate_codon_table` step. The default configuration uses:

```bash
module load Biopython/1.83-foss-2022b
```

Override the module without editing the wrapper with `BIOPYTHON_MODULE`, or skip
module loading with `BIOPYTHON_USE_MODULE=0` when the selected `PYTHON` already
has Biopython installed:

```bash
BIOPYTHON_MODULE=Biopython/1.83-foss-2022b \
  bash qc_analysis/scripts/run_qc_preprocessing.sh build_primate_codon_table config/qc_preprocessing.yaml

BIOPYTHON_USE_MODULE=0 \
  bash qc_analysis/scripts/run_qc_preprocessing.sh build_primate_codon_table config/qc_preprocessing.yaml
```

The wrapper performs a `from Bio import Entrez, SeqIO` preflight after any module
load and before it downloads or parses GenBank records. `PYTHON` remains
overrideable (for example, `PYTHON=/path/to/python`).

`sample_ref_file` requires `sample`; a headerless two-column `sample, species` file
is also supported for existing workflows. Useful optional fields are `species` and
`family`. Several samples can share one coordinate reference; they receive one map
row each, while codons are emitted once for that reference.

## Inferring accession from reference manifests

An accession in the first populated direct sample column is always preferred;
by default those columns are `accession`, `accession_version`, `reference_id`, and
`seq_name`. For a headerless two-column `sample_ref_file`, the script instead
infers the mitochondrial accession from
`references/manifests/reference_materialization_manifest.resolved.tsv` by matching
`sample_ref` `species` to manifest `target_species`. Matching is case-insensitive
and treats spaces and underscores equivalently (for example, `Cheracebus lugens`
and `Cheracebus_lugens`).

The preferred chrM annotation accession order is:

1. `final_chrM_genbank_accn`
2. `final_chrM_refseq_accn`
3. `final_chrM_accession`
4. `chrM_source_accession`

If the resolved manifest has no usable match, the configured materialization and
in-house manifests are tried next. The configured species FASTA directory is a
last fallback: its filename or first header can supply an accession, but GenBank is
still required to obtain CDS annotation. The summary table records the accession
source, manifest match, and FASTA path for auditing.

Validate mapping without downloading GenBank records:

```bash
module load Biopython/1.83-foss-2022b
python qc_analysis/scripts/build_primate_codon_table.py \
  --config config/qc_preprocessing.yaml \
  --dry-run
```

Run the full preprocessing step with:

```bash
module load Biopython/1.83-foss-2022b
bash qc_analysis/scripts/run_qc_preprocessing.sh \
  build_primate_codon_table \
  config/qc_preprocessing.yaml
```

Downloaded records are cached in `data/reference_tables/primate_genbank/<accession>.gb`.
The primary outputs are `data/reference_tables/reference_codon_table.tsv` and
`data/reference_tables/sample_reference_map.tsv`. The codon table includes
`reference_key`, original 1-based genomic `pos`, normalized/original gene names, genomic
reference base, coding-oriented codon and phase, all three genomic codon positions,
strand, and CDS qualifiers. Joined and minus-strand CDS features are emitted in
coding orientation while retaining their original genomic coordinates.

`reference_key` is constructed from the canonical coordinate-reference FASTA path,
coordinate-reference accession, and SHA256 of the normalized FASTA sequence. It is
therefore not species-based. `codon_match` first resolves `sample -> reference_key`,
then looks up `reference_key + pos`. Set
`write_legacy_sample_level_export: true` and configure `legacy_sample_codon_table`
only when an older consumer still needs the duplicated sample-level table; it is
disabled by default.

Final sample failures are recorded without stopping other samples in
`results/qc/codon_table_build/failed_genbank_downloads.tsv`; per-sample counts and
sanity-check warnings are recorded in the matching summary table. A GenBank
resolution or download problem is removed from this final-failure table if MITOS2
subsequently provides a usable fallback annotation. Set the optional
`build_primate_codon_table.settings.email` for NCBI Entrez requests. Use `--dry-run`
to resolve and audit every accession without downloading or parsing GenBank records,
or `--force-download` to refresh cached records.

## Parallel and HPC operation

The summary `status` is a final sample outcome:

* `completed` means GenBank CDS annotation successfully generated the codon table.
* `completed_mitos2_fallback` means GenBank annotation was unavailable, but MITOS2
  fallback successfully generated the sample codon table.
* `failed` means neither GenBank nor MITOS2 produced a usable codon table.

The final completed count includes both `completed` and
`completed_mitos2_fallback`. Before reporting it, the script also compares the
successful summary sample set with the output-table sample set and emits a warning
if they differ.

Accession and coordinate resolution remain serial. The script downloads each
unique missing GenBank accession at most once using `download_workers`, then parses
only eligible downloaded/cached GenBank files using `workers`. MITOS2 fallback
post-processing currently occurs in the main process. Thus, when every sample
requires MITOS2 fallback, the GenBank parsing thread pool has no useful work (and
is not created), although mixed runs still support parallel GenBank parsing.
Final TSV construction, comparisons, validation, and writing are deliberately
single-process, so the shared output files are written exactly once in
sample-reference order. Do **not** run several `--sample`
invocations concurrently against the same configuration: each invocation still
owns the same global output paths.

Use `--workers` for eligible cached/downloaded GenBank parsing. Its default is
CLI value, then YAML `settings.workers`, then `SLURM_CPUS_PER_TASK`, then 1.
`--download-workers` controls only concurrent Entrez downloads; without an API
key it defaults to 1. The process-wide rate limiter spaces all Entrez request
starts using `requests_per_second_without_api_key` (default 3) or
`requests_per_second_with_api_key` (default 10). `sleep_seconds` remains for
configuration compatibility, but the rate limiter is authoritative for parallel
downloads. Set `email` and, where permitted, `entrez_api_key` in private
configuration rather than committing credentials.

Recommended initial HPC settings are 4--8 workers, 16G memory, and one download
worker without an API key. With an API key, use at most 2--4 download workers
while retaining the global limiter. For example:

```bash
SLURM_CPUS=8 CODON_TABLE_WORKERS=8 \
  bash qc_analysis/scripts/run_qc_preprocessing.sh \
  --submit build_primate_codon_table config/qc_preprocessing.yaml
```

Cached `.gb` records substantially reduce runtime. More workers than unique
accessions or selected samples do not improve throughput.
