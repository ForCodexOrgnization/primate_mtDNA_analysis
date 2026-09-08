from qc_analysis.lib.heteroplasmy_clusters import ClusterConfig, expand_clusters, independent_seeds, max_af_coherent_count
from qc_analysis.lib.heteroplasmy_recurrence import assess_recurrence, build_species_sample_index
from qc_analysis.lib.numt_annotation import NumtInterval, best_cluster_overlap, best_species_overlap, merge_species_intervals


def test_af_coherent_seed_and_expansion():
    cfg = ClusterConfig(mt_length=1000, window_bp=100, af_span_max=0.06, min_seed_variants=3, simulations=10)
    positions = [100, 120, 150, 175, 700]
    afs = [0.11, 0.12, 0.13, 0.15, 0.40]
    assert max_af_coherent_count(positions, afs, cfg) == 4
    seeds = independent_seeds(positions, afs, threshold=3, cfg=cfg)
    assert seeds
    expanded = expand_clusters(seeds, positions, afs, cfg)
    assert len(expanded[0]) == 4


def test_sample_numt_overlap():
    rows = [
        {"source_pos": "100", "source_af": "0.11"},
        {"source_pos": "120", "source_af": "0.12"},
        {"source_pos": "150", "source_af": "0.13"},
    ]
    interval = NumtInterval(
        sample="A", species="Papio_anubis", reference_key="Papio_anubis",
        nuclear_chrom="chr1", nuclear_start=1000, nuclear_end=1200,
        chrm_start=90, chrm_end=160, tier="HIGH_CONF_NUMT",
    )
    hit = best_cluster_overlap(rows, [interval])
    assert hit["overlap_n"] == 3
    assert hit["overlap_fraction"] == 1.0


def test_species_numt_requires_other_sample_support():
    intervals = [
        NumtInterval(
            sample="B", species="Papio_anubis", reference_key="Papio_anubis",
            nuclear_chrom="chr1", nuclear_start=1000, nuclear_end=1200,
            chrm_start=90, chrm_end=160, tier="BESTHIT_ONLY_NUMT",
        )
    ]
    species = merge_species_intervals(intervals)
    rows = [
        {"source_pos": "100", "source_af": "0.11"},
        {"source_pos": "120", "source_af": "0.12"},
        {"source_pos": "150", "source_af": "0.13"},
    ]
    hit = best_species_overlap(rows, "A", "Papio_anubis", "Papio_anubis", species)
    assert hit["overlap_n"] == 3
    assert hit["other_support_samples"] == ["B"]


def test_recurrence_uses_same_species_pos_ref_alt_and_af():
    rows = [
        {"sample": "A", "source_chrom": "chrM", "source_pos": "100", "source_ref": "A", "source_alt": "G", "source_af": "0.11"},
        {"sample": "A", "source_chrom": "chrM", "source_pos": "120", "source_ref": "C", "source_alt": "T", "source_af": "0.12"},
        {"sample": "B", "source_chrom": "chrM", "source_pos": "100", "source_ref": "A", "source_alt": "G", "source_af": "0.12"},
        {"sample": "B", "source_chrom": "chrM", "source_pos": "120", "source_ref": "C", "source_alt": "T", "source_af": "0.13"},
        {"sample": "C", "source_chrom": "chrM", "source_pos": "100", "source_ref": "A", "source_alt": "G", "source_af": "0.40"},
    ]
    species = {"A": "Papio_anubis", "B": "Papio_anubis", "C": "Other_species"}
    index = build_species_sample_index(rows, species)
    result = assess_recurrence(rows[:2], "A", "Papio_anubis", index)
    assert result["recurrent"] is True
    assert result["best_recurrent_sample"] == "B"
    assert result["shared_variants"] == 2
