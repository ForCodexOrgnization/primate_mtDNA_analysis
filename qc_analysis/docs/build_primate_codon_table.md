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
is also supported for existing workflows. To build annotations, each sample needs an
accession in the first populated configured column, by default `accession`,
`accession_version`, `reference_id`, then `seq_name`. Useful optional fields are
`species` and `family`. Several samples may use one accession and each receives its
own rows so `run_codon_match.py` can look up `sample + pos`.

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
to report planned downloads without creating the output table, or `--force-download`
to refresh cached records.
