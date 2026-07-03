#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
sample_file <- ifelse(length(args) >= 1, args[[1]], "data/metadata/sample_metadata.tsv")
ref_file <- ifelse(length(args) >= 2, args[[2]], "references/manifests/in_house_score_reference_inputs.tsv")
score_file <- ifelse(length(args) >= 3, args[[3]], "results/preprocessing/in_house_score/merged_in_house_score.tsv")
out_file <- ifelse(length(args) >= 4, args[[4]], "results/preprocessing/variant_calling_inputs/variant_calling_input_table.tsv")
source("preprocessing/R/read_preprocessing_inputs.R")
source("preprocessing/R/reference_manifest_helpers.R")
samples <- read_tsv_flexible(sample_file); refs <- read_tsv_flexible(ref_file)
if (!"target_species" %in% names(samples)) samples$target_species <- if ("species" %in% names(samples)) samples$species else ""
samples$join_species <- normalize_species(samples$target_species); refs$join_species <- normalize_species(refs$target_species)
vc <- merge(samples, refs, by = "join_species", all.x = TRUE, suffixes = c("", ".ref"))
if (file.exists(score_file) && file.info(score_file)$size > 0) {
  scores <- read_tsv_flexible(score_file)
  if (!"numt_mask_path" %in% names(scores)) {
    scores$numt_mask_path <- if ("REF_TYPE" %in% names(scores) && "MinimalMaskBED" %in% names(scores)) {
      eligible <- scores$REF_TYPE %in% c("#C-likely_comp", "#C-Ambiguous")
      has_final <- if ("MaskPriority" %in% names(scores)) grepl("FINAL_minimal_mask", scores$MaskPriority) else rep(FALSE, nrow(scores))
      ifelse(eligible & has_final, scores$MinimalMaskBED, "")
    } else {
      ""
    }
  }
  if (!"minimal_numt_mask_status" %in% names(scores)) {
    scores$minimal_numt_mask_status <- if ("MaskPriority" %in% names(scores)) scores$MaskPriority else "completed"
  }
  if ("sample_id" %in% names(scores) && any(scores$sample_id %in% vc$sample_id)) {
    vc <- merge(vc, scores, by = "sample_id", all.x = TRUE, suffixes = c("", ".score"))
  } else if ("Species" %in% names(scores)) {
    scores$join_species <- normalize_species(scores$Species)
    vc <- merge(vc, scores, by = "join_species", all.x = TRUE, suffixes = c("", ".score"))
  } else if ("target_species" %in% names(scores)) {
    scores$join_species <- normalize_species(scores$target_species)
    vc <- merge(vc, scores, by = "join_species", all.x = TRUE, suffixes = c("", ".score"))
  }
}

variant_ref_manifest <- "references/variant_calling/variant_calling_reference_manifest.tsv"
variant_ref_cols <- c("safe_species_id","whole_fasta","whole_fai","whole_dict","nuclear_fasta","nuclear_fai","chrM_fasta","chrM_fai","chrM_dict","chrM_shift_fasta","chrM_shift_fai","chrM_shift_dict","non_control_interval","control_region_shifted_interval","shift_back_chain","numt_mask_applied_to_whole_ref","build_status","build_message")
if (file.exists(variant_ref_manifest) && file.info(variant_ref_manifest)$size > 0) {
  variant_refs <- read_tsv_flexible(variant_ref_manifest)
  if ("target_species" %in% names(variant_refs)) {
    variant_refs$join_species <- normalize_species(variant_refs$target_species)
    keep_cols <- unique(c("join_species", variant_ref_cols[variant_ref_cols %in% names(variant_refs)]))
    variant_refs <- variant_refs[, keep_cols, drop = FALSE]
    variant_refs <- variant_refs[!duplicated(variant_refs$join_species), , drop = FALSE]
    vc <- merge(vc, variant_refs, by = "join_species", all.x = TRUE, suffixes = c("", ".variant_ref"))
  }
}
vc$target_species <- coalesce_chr(vc$target_species, vc$target_species.ref)
vc <- ensure_columns(vc, c("sample_id","cram_path","cram_index_path","final_wg_ref_species","final_wg_assembly_accession","wg_fasta_path","wg_fai_path","final_chrM_species","final_chrM_accession","chrM_fasta_path","chrM_fai_path","chrM_reference_context","reference_pairing_status","final_reference_strategy","numt_mask_path","minimal_numt_mask_status","sample_qc_status","manual_review_required","manual_review_reason", variant_ref_cols))
vc$sample_qc_status <- ifelse(vc$sample_qc_status == "", "pending", vc$sample_qc_status)
vc$minimal_numt_mask_status <- ifelse(vc$minimal_numt_mask_status == "", "pending", vc$minimal_numt_mask_status)
if ("build_status" %in% names(vc)) {
  bad_ref <- vc$build_status != "" & vc$build_status != "success"
  vc$sample_qc_status <- ifelse(bad_ref & vc$sample_qc_status %in% c("", "pending"), "reference_package_review", vc$sample_qc_status)
  vc$manual_review_required <- ifelse(bad_ref, "TRUE", vc$manual_review_required)
  vc$manual_review_reason <- ifelse(bad_ref, coalesce_chr(vc$manual_review_reason, paste0("variant_reference_package_", vc$build_status, ":", vc$build_message)), vc$manual_review_reason)
}
cols <- c("sample_id","target_species","cram_path","cram_index_path","final_wg_ref_species","final_wg_assembly_accession","wg_fasta_path","wg_fai_path","final_chrM_species","final_chrM_accession","chrM_fasta_path","chrM_fai_path","chrM_reference_context","reference_pairing_status","final_reference_strategy","numt_mask_path","minimal_numt_mask_status", variant_ref_cols, "sample_qc_status","manual_review_required","manual_review_reason")
write_tsv(vc[, cols], out_file)
message("Wrote ", out_file)
