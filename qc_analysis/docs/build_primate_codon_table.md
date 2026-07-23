# Build primate codon table

`build_primate_codon_table.py` downloads (or reuses) GenBank records and writes the
sample-level coding-position table used by `codon_match`:

```bash
python qc_analysis/scripts/build_primate_codon_table.py --config config/qc_preprocessing.yaml
bash qc_analysis/scripts/run_qc_preprocessing.sh build_primate_codon_table config/qc_preprocessing.yaml
```

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
`family`. Several samples may use one accession and each receives its own rows so
`run_codon_match.py` can look up `sample + pos`.

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
The output `data/reference_tables/all_primate_position_codon_table.tsv` includes the
sample, original 1-based genomic `pos`, normalized/original gene names, genomic
reference base, coding-oriented codon and phase, all three genomic codon positions,
strand, and CDS qualifiers. Joined and minus-strand CDS features are emitted in
coding orientation while retaining their original genomic coordinates.

Failures are recorded without stopping other samples in
`results/qc/codon_table_build/failed_genbank_downloads.tsv`; per-sample counts and
sanity-check warnings are recorded in the matching summary table. Set the optional
`build_primate_codon_table.settings.email` for NCBI Entrez requests. Use `--dry-run`
to resolve and audit every accession without downloading or parsing GenBank records,
or `--force-download` to refresh cached records.

## Parallel and HPC operation

The script resolves accessions once, downloads each unique missing accession at
most once, and then prepares/parses samples internally with worker threads.
Final TSV construction, MITOS2 fallback selection, comparisons, validation, and
writing are deliberately single-process, so the shared output files are written
exactly once in sample-reference order. Do **not** run several `--sample`
invocations concurrently against the same configuration: each invocation still
owns the same global output paths.

Use `--workers` for sample preparation and cached GenBank parsing. Its default is
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
