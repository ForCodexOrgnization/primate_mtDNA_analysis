#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
})

# Intra-species contamination detection for primate mtDNA samples.
#
# The analysis is intentionally performed before coordinate liftover because
# source and recipient samples must be compared in their native, species-level
# coordinate system. Variants are therefore matched only within the same
# Species value and by CHROM/POS/REF/ALT.

parse_args <- function(args) {
  defaults <- list(
    variant_file = NULL,
    out_dir = "results/qc/intra_species_contamination",
    dp_min = 100,
    use_snv_only = TRUE,
    low_vaf_min = 0.01,
    low_vaf_max = 0.20,
    high_vaf_min = 0.99,
    mt_lower = 0.80,
    mt_upper = 0.998,
    min_low_count = 5,
    min_overlap = 3,
    candidate_overlap_fraction = 0.50,
    highconf_overlap_fraction = 0.6213636363636358,
    candidate_contamination = 0.036420574377757434,
    highconf_contamination = 0.07103935483870959
  )

  usage <- function() {
    cat(paste(
      "Usage:",
      "  Rscript qc_analysis/scripts/run_intra_species_contamination.R \\",
      "    --variant-file FILE [--out-dir DIR] [options]",
      "",
      "Required:",
      "  --variant-file FILE    Combined variant table containing Sample, Species,",
      "                         CHROM, POS, REF, ALT, FILTER, DP and VAF or AF.",
      "",
      "Options:",
      "  --out-dir DIR          Default: results/qc/intra_species_contamination",
      "  --dp-min N             Default: 100",
      "  --use-snv-only BOOL    Default: true",
      "  --low-vaf-min FLOAT    Default: 0.01",
      "  --low-vaf-max FLOAT    Default: 0.20",
      "  --high-vaf-min FLOAT   Default: 0.99",
      "  --help",
      sep = "\n"
    ))
  }

  key_map <- c(
    "variant-file" = "variant_file",
    "out-dir" = "out_dir",
    "dp-min" = "dp_min",
    "use-snv-only" = "use_snv_only",
    "low-vaf-min" = "low_vaf_min",
    "low-vaf-max" = "low_vaf_max",
    "high-vaf-min" = "high_vaf_min",
    "mt-lower" = "mt_lower",
    "mt-upper" = "mt_upper",
    "min-low-count" = "min_low_count",
    "min-overlap" = "min_overlap",
    "candidate-overlap-fraction" = "candidate_overlap_fraction",
    "highconf-overlap-fraction" = "highconf_overlap_fraction",
    "candidate-contamination" = "candidate_contamination",
    "highconf-contamination" = "highconf_contamination"
  )

  if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
    usage()
    quit(status = 0)
  }

  i <- 1
  while (i <= length(args)) {
    token <- args[[i]]
    if (!startsWith(token, "--")) stop("Unexpected argument: ", token)
    key <- sub("^--", "", token)
    if (!key %in% names(key_map)) stop("Unknown option: ", token)
    if (i == length(args)) stop("Missing value for ", token)
    field <- unname(key_map[[key]])
    value <- args[[i + 1]]
    current <- defaults[[field]]
    if (is.logical(current)) {
      value <- tolower(value)
      if (!value %in% c("true", "false", "1", "0", "yes", "no")) {
        stop("Expected true/false for ", token)
      }
      defaults[[field]] <- value %in% c("true", "1", "yes")
    } else if (is.numeric(current)) {
      parsed <- suppressWarnings(as.numeric(value))
      if (is.na(parsed)) stop("Expected numeric value for ", token)
      defaults[[field]] <- parsed
    } else {
      defaults[[field]] <- value
    }
    i <- i + 2
  }

  if (is.null(defaults$variant_file) || defaults$variant_file == "") {
    stop("--variant-file is required")
  }
  defaults
}

safe_mean <- function(x) if (length(x) == 0) NA_real_ else mean(x, na.rm = TRUE)
safe_median <- function(x) if (length(x) == 0) NA_real_ else median(x, na.rm = TRUE)

p <- parse_args(commandArgs(trailingOnly = TRUE))
list2env(p, envir = .GlobalEnv)

if (!file.exists(variant_file)) stop("Variant table not found: ", variant_file)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading variants: ", variant_file)
raw <- fread(variant_file) %>%
  mutate(
    Sample = as.character(Sample),
    Species = as.character(Species),
    CHROM = as.character(CHROM),
    POS = suppressWarnings(as.integer(POS)),
    REF = as.character(REF),
    ALT = as.character(ALT),
    FILTER = as.character(FILTER),
    Type = if ("Type" %in% names(.)) as.character(Type) else NA_character_,
    DP = suppressWarnings(as.numeric(DP)),
    VAF = if ("VAF" %in% names(.)) suppressWarnings(as.numeric(VAF))
      else if ("AF" %in% names(.)) suppressWarnings(as.numeric(AF))
      else stop("Input table has neither VAF nor AF"),
    variant_id = paste(CHROM, POS, REF, ALT, sep = ":")
  ) %>%
  filter(!is.na(Sample), Sample != "", !is.na(Species), Species != "")

required <- c("Sample", "Species", "CHROM", "POS", "REF", "ALT", "FILTER", "DP")
missing_required <- setdiff(required, names(raw))
if (length(missing_required) > 0) {
  stop("Missing required columns: ", paste(missing_required, collapse = ", "))
}

all_samples <- raw %>% distinct(Species, Sample)
species_counts <- all_samples %>% count(Species, name = "n_species_samples")
all_samples <- all_samples %>% left_join(species_counts, by = "Species")

variants <- raw %>%
  filter(
    FILTER == "PASS",
    !is.na(POS), !is.na(REF), !is.na(ALT),
    !is.na(DP), DP >= dp_min,
    !is.na(VAF), VAF >= 0, VAF <= 1
  )

if (use_snv_only) {
  variants <- variants %>%
    filter(nchar(REF) == 1, nchar(ALT) == 1, is.na(Type) | Type == "SNV")
}

low <- variants %>%
  filter(VAF >= low_vaf_min, VAF <= low_vaf_max) %>%
  select(Species, Sample, CHROM, POS, REF, ALT, variant_id, low_VAF = VAF, low_DP = DP)

high <- variants %>%
  filter(VAF >= high_vaf_min) %>%
  select(Species, Sample, CHROM, POS, REF, ALT, variant_id, high_VAF = VAF, high_DP = DP)

low_summary <- low %>%
  group_by(Species, Sample) %>%
  summarise(
    n_lowA = n_distinct(variant_id),
    mean_lowA_VAF = mean(low_VAF, na.rm = TRUE),
    median_lowA_VAF = median(low_VAF, na.rm = TRUE),
    .groups = "drop"
  )

high_summary <- high %>%
  group_by(Species, Sample) %>%
  summarise(n_high = n_distinct(variant_id), .groups = "drop")

pairwise_variant <- low %>%
  rename(Sample_A = Sample, VAF_A = low_VAF, DP_A = low_DP) %>%
  inner_join(
    high %>% rename(Sample_B = Sample, VAF_B = high_VAF, DP_B = high_DP),
    by = c("Species", "CHROM", "POS", "REF", "ALT", "variant_id")
  ) %>%
  filter(Sample_A != Sample_B)

pairwise_stats <- pairwise_variant %>%
  group_by(Species, Sample_A, Sample_B) %>%
  summarise(
    overlap = n_distinct(variant_id),
    median_overlap_lowA_VAF = median(VAF_A, na.rm = TRUE),
    median_source_highB_VAF = median(VAF_B, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(low_summary %>% select(Species, Sample_A = Sample, n_lowA),
            by = c("Species", "Sample_A")) %>%
  mutate(frac_lowA_in_highB = ifelse(n_lowA > 0, overlap / n_lowA, NA_real_)) %>%
  arrange(Species, Sample_A, desc(frac_lowA_in_highB), desc(overlap))

best_source <- pairwise_stats %>%
  group_by(Species, Sample_A) %>%
  slice_max(order_by = frac_lowA_in_highB, n = 1, with_ties = TRUE) %>%
  arrange(desc(overlap), Sample_B, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    Species,
    Sample = Sample_A,
    best_source_sample = Sample_B,
    best_overlap = overlap,
    best_frac_lowA_in_highB = frac_lowA_in_highB,
    best_median_overlap_lowA_VAF = median_overlap_lowA_VAF,
    best_median_source_highB_VAF = median_source_highB_VAF
  )

# Leave-one-out mt-high-hets estimate. A high-VAF allele observed in any other
# sample of the same species is an anchor. If that allele is depressed to
# 0.80-0.998 in sample A, 1 - mean(VAF) estimates the minor mixture fraction.
estimate_one <- function(species_i, sample_i) {
  anchors <- high %>%
    filter(Species == species_i, Sample != sample_i) %>%
    group_by(Species, CHROM, POS, REF, ALT, variant_id) %>%
    summarise(anchor_support_n_samples = n_distinct(Sample), .groups = "drop")

  if (nrow(anchors) == 0) {
    return(tibble(
      Species = species_i, Sample = sample_i,
      n_anchor_pool_excluding_A = 0L, n_anchor_tested_in_A = 0L,
      n_depressed_anchor = 0L, n_fallback_anchor = 0L,
      mean_depressed_VAF = NA_real_, mean_fallback_anchor_VAF = NA_real_,
      mt_high_hets_contamination = NA_real_, mt_high_hets_mode = "no_LOO_anchor_pool"
    ))
  }

  observed <- anchors %>%
    inner_join(
      variants %>%
        filter(Species == species_i, Sample == sample_i) %>%
        select(Species, CHROM, POS, REF, ALT, variant_id, VAF_A = VAF),
      by = c("Species", "CHROM", "POS", "REF", "ALT", "variant_id")
    )

  depressed <- observed$VAF_A[observed$VAF_A >= mt_lower & observed$VAF_A <= mt_upper]
  fallback <- observed$VAF_A[observed$VAF_A >= mt_lower & observed$VAF_A <= 1]

  contamination <- if (length(depressed) >= 3) 1 - safe_mean(depressed)
    else if (length(fallback) > 0) 1 - safe_mean(fallback)
    else NA_real_
  mode <- if (length(depressed) >= 3) "depressed_anchors"
    else if (length(fallback) > 0) "fallback_high_VAF_anchors"
    else "no_high_VAF_anchor_observed_in_A"

  tibble(
    Species = species_i,
    Sample = sample_i,
    n_anchor_pool_excluding_A = nrow(anchors),
    n_anchor_tested_in_A = nrow(observed),
    n_depressed_anchor = length(depressed),
    n_fallback_anchor = length(fallback),
    mean_depressed_VAF = safe_mean(depressed),
    mean_fallback_anchor_VAF = safe_mean(fallback),
    mt_high_hets_contamination = contamination,
    mt_high_hets_mode = mode
  )
}

mt_high_hets <- pmap_dfr(
  list(all_samples$Species, all_samples$Sample),
  function(Species, Sample) estimate_one(Species, Sample)
)

usable_summary <- variants %>%
  group_by(Species, Sample) %>%
  summarise(n_usable_variants = n_distinct(variant_id), .groups = "drop")

final <- all_samples %>%
  left_join(usable_summary, by = c("Species", "Sample")) %>%
  left_join(low_summary, by = c("Species", "Sample")) %>%
  left_join(high_summary, by = c("Species", "Sample")) %>%
  left_join(best_source, by = c("Species", "Sample")) %>%
  left_join(mt_high_hets, by = c("Species", "Sample")) %>%
  mutate(
    across(c(n_usable_variants, n_lowA, n_high, best_overlap), ~replace_na(.x, 0L)),
    best_frac_lowA_in_highB = replace_na(best_frac_lowA_in_highB, 0),
    pass_lowA_count = n_lowA >= min_low_count,
    pass_overlap_count = best_overlap >= min_overlap,
    pass_overlap_fraction_candidate = best_frac_lowA_in_highB >= candidate_overlap_fraction,
    pass_overlap_fraction_highconf = best_frac_lowA_in_highB >= highconf_overlap_fraction,
    pass_mt_high_hets_candidate = !is.na(mt_high_hets_contamination) &
      mt_high_hets_contamination >= candidate_contamination,
    pass_mt_high_hets_highconf = !is.na(mt_high_hets_contamination) &
      mt_high_hets_contamination >= highconf_contamination,
    contamination_flag_candidate = pass_lowA_count & pass_overlap_count &
      pass_overlap_fraction_candidate & pass_mt_high_hets_candidate,
    contamination_flag_highconf = pass_lowA_count & pass_overlap_count &
      pass_overlap_fraction_highconf & pass_mt_high_hets_highconf,
    contamination_flag = contamination_flag_candidate,
    contamination_status = case_when(
      contamination_flag_highconf ~ "high_confidence_contaminated",
      contamination_flag_candidate ~ "candidate_contaminated",
      n_usable_variants == 0 ~ "insufficient_variant_data",
      n_species_samples == 1 ~ "insufficient_singleton_species",
      is.na(mt_high_hets_contamination) ~ "insufficient_anchor_data",
      pass_lowA_count & pass_overlap_count & pass_overlap_fraction_candidate ~ "lowA_highB_overlap_only",
      pass_mt_high_hets_candidate ~ "mt_high_hets_only",
      TRUE ~ "no_strong_evidence"
    )
  ) %>%
  arrange(desc(contamination_flag_highconf), desc(contamination_flag_candidate),
          desc(mt_high_hets_contamination), desc(best_frac_lowA_in_highB))

write_tsv(all_samples, file.path(out_dir, "all_samples.tsv"))
write_tsv(pairwise_variant, file.path(out_dir, "pairwise_lowA_highB_overlap_variant_level.tsv"))
write_tsv(pairwise_stats, file.path(out_dir, "pairwise_lowA_highB_overlap_stats.tsv"))
write_tsv(mt_high_hets, file.path(out_dir, "mt_high_hets_all_excluding_A_summary.tsv"))
write_tsv(final, file.path(out_dir, "final_contamination_summary.tsv"))
write_tsv(final %>% filter(contamination_flag_candidate),
          file.path(out_dir, "flagged_candidate_contaminated_samples.tsv"))
write_tsv(final %>% filter(contamination_flag_highconf),
          file.path(out_dir, "flagged_high_confidence_contaminated_samples.tsv"))
write_tsv(final %>% count(contamination_status, name = "n_samples"),
          file.path(out_dir, "contamination_status_summary.tsv"))

thresholds <- tibble(
  dp_min, use_snv_only, low_vaf_min, low_vaf_max, high_vaf_min,
  mt_lower, mt_upper, min_low_count, min_overlap,
  candidate_overlap_fraction, highconf_overlap_fraction,
  candidate_contamination, highconf_contamination
)
write_tsv(thresholds, file.path(out_dir, "thresholds_used.tsv"))

message("Done. Output directory: ", out_dir)
message("Candidate contaminated samples: ", sum(final$contamination_flag_candidate, na.rm = TRUE))
message("High-confidence contaminated samples: ", sum(final$contamination_flag_highconf, na.rm = TRUE))
