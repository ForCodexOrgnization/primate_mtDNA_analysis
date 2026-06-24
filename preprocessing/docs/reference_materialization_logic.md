# Reference materialization logic

Materialization converts reviewed reference choices into local FASTA files. NCBI WG rows use `final_wg_ftp_path` to download `{basename}_genomic.fna.gz` and `{basename}_assembly_report.txt`, standardized to `references/wg/{assembly_accession}/{assembly_accession}.genome.fa.gz` and `.assembly_report.txt`.

Embedded chrM references are extracted from the local WG FASTA and indexed. If the selected NCBI GCA/GCF assembly partner does not contain the expected chrM record, materialization looks up the paired assembly from the NCBI assembly summaries, removes the non-chrM-bearing local WG/chrM files, downloads the paired assembly, and retries extraction from that partner. Independent chrM references are extracted from a local mitochondrion FASTA when possible or downloaded from NCBI nucleotide. DNA Zoo automatic download is not implemented and such rows are flagged for manual review.
