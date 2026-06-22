# Reference classification logic

`chrM_reference_context` is `embedded_in_wg_ref` when `final_chrM_source` contains `whole_genome_assembly`, `final_chrM_assembly_accession` equals `final_wg_assembly_accession`, and `final_chrM_length` is 14,000-25,000 bp. It is `independent_chrM_ref` when a complete chrM accession exists but does not come from the final WG assembly. It is `missing_chrM_ref` when the accession is blank or length is outside 14,000-25,000 bp.

`reference_pairing_status` compares normalized target, WG, and chrM species names to distinguish same-species and cross-species pairings. Rows are manually reviewed when references are missing, chrM length is invalid, only WG or chrM exists, no reference exists, DNA Zoo download is required, or materialization fails.

Embedded chrM must be extracted from the WG FASTA to keep the NUMT score reference pair internally consistent. Short mitochondrial-like fragments such as 204 bp records are never complete chrM references.
