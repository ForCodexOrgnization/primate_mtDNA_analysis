validate_chrM_context_paths <- function(df) {
  bad <- with(df, chrM_reference_context == "embedded_in_wg_ref" & !grepl("references/chrM/embedded_from_wg/", chrM_fasta_path)) |
    with(df, chrM_reference_context == "independent_chrM_ref" & !grepl("references/chrM/independent/", chrM_fasta_path))
  df$manual_review_required[bad] <- "yes"
  df$manual_review_reason[bad] <- ifelse(df$manual_review_reason[bad] == "", "chrM_context_path_mismatch", paste(df$manual_review_reason[bad], "chrM_context_path_mismatch", sep = ";"))
  df
}
