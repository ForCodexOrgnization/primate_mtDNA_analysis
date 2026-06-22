# In-house score inputs and outputs

Input: `references/manifests/in_house_score_reference_inputs.tsv`, joined to sample metadata so each job has `target_species`, `sample_id`, `cram_path`, `wg_fasta_path`, `chrM_fasta_path`, `chrM_reference_context`, and `reference_pairing_status`.

Output: per-sample score/mask files and `results/preprocessing/in_house_score/merged_in_house_score.tsv`, including the selected minimal NUMT mask path and status.
