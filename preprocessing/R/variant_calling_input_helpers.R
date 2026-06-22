combine_review_flags <- function(...) {
  vals <- unlist(list(...)); vals <- vals[!is.na(vals) & vals != "" & vals != "no"]
  if (length(vals)) "yes" else "no"
}
