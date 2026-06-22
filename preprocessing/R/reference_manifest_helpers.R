`%||%` <- function(x, y) if (is.null(x)) y else x
blank_to_na <- function(x) { x <- as.character(x); x[is.na(x) | trimws(x) == ""] <- NA_character_; x }
coalesce_chr <- function(...) { args <- lapply(list(...), blank_to_na); out <- args[[1]]; for (a in args[-1]) out[is.na(out)] <- a[is.na(out)]; out[is.na(out)] <- ""; out }
normalize_species <- function(x) tolower(gsub("[^a-z0-9]+", "_", trimws(blank_to_na(x))))
safe_token <- function(x) { x <- normalize_species(x); x[is.na(x) | x == ""] <- "unknown"; x }
mito_len_ok <- function(x, min_len = 14000, max_len = 25000) suppressWarnings(!is.na(as.numeric(x)) & as.numeric(x) >= min_len & as.numeric(x) <= max_len)
classify_chrM_context <- function(source, chrM_acc, chrM_asm, wg_asm, chrM_len, min_len = 14000, max_len = 25000) {
  src <- tolower(coalesce_chr(source)); acc <- coalesce_chr(chrM_acc); len_ok <- mito_len_ok(chrM_len, min_len, max_len)
  same_asm <- coalesce_chr(chrM_asm) != "" & coalesce_chr(wg_asm) != "" & coalesce_chr(chrM_asm) == coalesce_chr(wg_asm)
  ifelse(acc == "" | !len_ok, "missing_chrM_ref", ifelse(grepl("whole_genome_assembly", src) & same_asm, "embedded_in_wg_ref", "independent_chrM_ref"))
}
reference_pairing_status <- function(target, wg_species, chrM_species) {
  t <- normalize_species(target); wg <- normalize_species(wg_species); mt <- normalize_species(chrM_species)
  has_wg <- !is.na(wg) & wg != ""; has_mt <- !is.na(mt) & mt != ""
  ifelse(has_wg & has_mt & wg == t & mt == t, "same_species_wg_same_species_chrM",
  ifelse(has_wg & has_mt & wg == t & mt != t, "same_species_wg_cross_species_chrM",
  ifelse(has_wg & has_mt & wg != t & mt == t, "cross_species_wg_same_species_chrM",
  ifelse(has_wg & has_mt & wg != t & mt != t, "cross_species_wg_cross_species_chrM",
  ifelse(has_wg & !has_mt, "wg_only_no_chrM", ifelse(!has_wg & has_mt, "chrM_only_no_wg", "no_reference_found"))))))
}
manual_review_reasons <- function(df, min_len = 14000, max_len = 25000) {
  apply(df, 1, function(r) {
    reasons <- character(); get <- function(n) if (n %in% names(r)) as.character(r[[n]]) else ""
    if (is.na(get("final_wg_assembly_accession")) || trimws(get("final_wg_assembly_accession")) == "") reasons <- c(reasons, "missing_wg_assembly_accession")
    if (is.na(get("final_chrM_accession")) || trimws(get("final_chrM_accession")) == "") reasons <- c(reasons, "missing_chrM_accession")
    len <- suppressWarnings(as.numeric(get("final_chrM_length"))); if (is.na(len) || len < min_len || len > max_len) reasons <- c(reasons, "chrM_length_outside_14000_25000")
    if (get("chrM_reference_context") == "missing_chrM_ref") reasons <- c(reasons, "missing_chrM_ref")
    if (get("reference_pairing_status") %in% c("wg_only_no_chrM", "chrM_only_no_wg", "no_reference_found")) reasons <- c(reasons, get("reference_pairing_status"))
    if (grepl("dnazoo|dna zoo", get("final_wg_ref_source"), ignore.case = TRUE)) reasons <- c(reasons, "dnazoo_download_not_implemented")
    paste(unique(reasons), collapse = ";")
  })
}
