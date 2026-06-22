read_tsv_flexible <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) stop("Missing or empty TSV: ", path)
  read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, quote = "")
}
write_tsv <- function(x, path) { dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE); write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "") }
ensure_columns <- function(df, cols) { for (nm in cols) if (!nm %in% names(df)) df[[nm]] <- ""; df }
