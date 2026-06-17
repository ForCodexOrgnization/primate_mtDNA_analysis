# Reusable helpers for human contamination QC reports.

read_optional_table <- function(path) {
  if (!file.exists(path)) {
    return(data.frame())
  }

  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

normalize_low_vaf <- function(vaf, low_vaf_min = 0.01, low_vaf_max = 0.05) {
  !is.na(vaf) & vaf >= low_vaf_min & vaf <= low_vaf_max
}

classify_contamination <- function(qc_table,
                                   score_col = "contamination_score",
                                   status_col = "contamination_status",
                                   warn_threshold = 0.05,
                                   fail_threshold = 0.10) {
  if (nrow(qc_table) == 0) {
    return(qc_table)
  }

  if (!score_col %in% names(qc_table)) {
    qc_table[[score_col]] <- NA_real_
  }

  score <- suppressWarnings(as.numeric(qc_table[[score_col]]))
  qc_table[[status_col]] <- ifelse(
    is.na(score),
    "not_evaluated",
    ifelse(score >= fail_threshold, "contaminated", ifelse(score >= warn_threshold, "review", "pass"))
  )
  qc_table
}

summarize_status <- function(qc_table, status_col = "contamination_status") {
  if (nrow(qc_table) == 0 || !status_col %in% names(qc_table)) {
    return(data.frame(contamination_status = character(), n = integer()))
  }

  counts <- as.data.frame(table(qc_table[[status_col]]), stringsAsFactors = FALSE)
  names(counts) <- c("contamination_status", "n")
  counts
}

plot_status_counts <- function(status_summary) {
  if (nrow(status_summary) == 0) {
    plot.new()
    text(0.5, 0.5, "No contamination status table available")
    return(invisible(NULL))
  }

  barplot(
    height = status_summary$n,
    names.arg = status_summary$contamination_status,
    las = 2,
    ylab = "Sample count",
    main = "Final contamination status"
  )
}
