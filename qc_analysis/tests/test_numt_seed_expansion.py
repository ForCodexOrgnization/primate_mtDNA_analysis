import copy
import unittest
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from qc_analysis.scripts.expand_numt_seed_variants import apply_expansion


class NumtSeedExpansionTests(unittest.TestCase):
    def setUp(self):
        self.cluster_rows = [
            {
                "sample": "S1",
                "cluster_id": "S1_C1",
                "cluster_start": "1968",
                "cluster_end": "2121",
                "wraps_origin": "False",
                "median_af": "0.110",
                "af_min": "0.104",
                "af_max": "0.159",
                "af_span": "0.055",
                "sample_numt_evidence": "YES",
                "species_numt_evidence": "NO",
                "species_numt_support_n": "0",
                "numt_overlap_n": "5",
                "numt_overlap_fraction": "1.0",
                "recurrence": "NO",
                "numt_chrm_start": "1941",
                "numt_chrm_end": "2115",
                "numt_tier": "BESTHIT_ONLY_NUMT",
                "cluster_class": "NUMT_ONLY",
                "filter_action": "REMOVE",
                "filter_reason": "SAMPLE_NUMT_OVERLAP",
                "numt_scope": "SAMPLE",
            }
        ]
        self.variant_rows = [
            {
                "sample": "S1", "species": "sp", "source_chrom": "chrM",
                "source_pos": "2057", "source_ref": "C", "source_alt": "T", "source_af": "0.171",
                "filter_action": "KEEP", "filter_reason": "", "numt_scope": "NONE", "numt_tier": "NONE",
                "cluster_id": "", "cluster_class": "", "clustered": "NO", "recurrence": "NO",
            },
            {
                "sample": "S1", "species": "sp", "source_chrom": "chrM",
                "source_pos": "2071", "source_ref": "A", "source_alt": "T", "source_af": "0.166",
                "filter_action": "KEEP", "filter_reason": "", "numt_scope": "NONE", "numt_tier": "NONE",
                "cluster_id": "", "cluster_class": "", "clustered": "NO", "recurrence": "NO",
            },
            {
                "sample": "S1", "species": "sp", "source_chrom": "chrM",
                "source_pos": "2150", "source_ref": "G", "source_alt": "A", "source_af": "0.150",
                "filter_action": "KEEP", "filter_reason": "", "numt_scope": "NONE", "numt_tier": "NONE",
                "cluster_id": "", "cluster_class": "", "clustered": "NO", "recurrence": "NO",
            },
        ]

    def test_expands_threshold_edge_variants_but_not_outside_interval(self):
        clusters = copy.deepcopy(self.cluster_rows)
        variants = copy.deepcopy(self.variant_rows)
        diagnostics = apply_expansion(
            clusters, variants, max_distance_bp=250, max_delta_af=0.07, max_total_af_span=0.10
        )
        accepted = {d["source_pos"] for d in diagnostics if d["accepted"] == "YES"}
        self.assertEqual(accepted, {2057, 2071})
        by_pos = {int(v["source_pos"]): v for v in variants}
        self.assertEqual(by_pos[2057]["filter_reason"], "NUMT_SEED_EXPANSION")
        self.assertEqual(by_pos[2071]["numt_expansion_seed_id"], "S1_C1")
        self.assertEqual(by_pos[2150]["filter_action"], "KEEP")

    def test_species_numt_without_recurrence_needs_strong_evidence(self):
        clusters = copy.deepcopy(self.cluster_rows)
        clusters[0].update({
            "sample_numt_evidence": "NO",
            "species_numt_evidence": "YES",
            "species_numt_support_n": "1",
            "numt_overlap_n": "2",
            "numt_overlap_fraction": "0.75",
            "recurrence": "NO",
            "filter_action": "FLAG",
            "filter_reason": "SPECIES_NUMT_OVERLAP",
            "numt_scope": "SPECIES",
        })
        variants = copy.deepcopy(self.variant_rows)
        for v in variants:
            v["cluster_id"] = "S1_C1"
            v["clustered"] = "YES"
            v["filter_action"] = "FLAG"
            v["filter_reason"] = "SPECIES_NUMT_OVERLAP"
        diagnostics = apply_expansion(clusters, variants)
        self.assertEqual(diagnostics, [])
        self.assertTrue(all(v["filter_action"] == "FLAG" for v in variants))
        self.assertEqual(clusters[0]["strong_species_numt_evidence"], "NO")

    def test_strong_species_numt_without_recurrence_is_removed(self):
        clusters = copy.deepcopy(self.cluster_rows)
        clusters[0].update({
            "sample_numt_evidence": "NO",
            "species_numt_evidence": "YES",
            "species_numt_support_n": "2",
            "numt_overlap_n": "3",
            "numt_overlap_fraction": "0.80",
            "recurrence": "NO",
            "filter_action": "FLAG",
            "filter_reason": "SPECIES_NUMT_OVERLAP",
            "numt_scope": "SPECIES",
        })
        variants = copy.deepcopy(self.variant_rows)
        for v in variants:
            v["cluster_id"] = "S1_C1"
            v["clustered"] = "YES"
            v["filter_action"] = "FLAG"
            v["filter_reason"] = "SPECIES_NUMT_OVERLAP"
            v["numt_scope"] = "SPECIES"
        apply_expansion(clusters, variants)
        self.assertEqual(clusters[0]["strong_species_numt_evidence"], "YES")
        self.assertEqual(clusters[0]["filter_action"], "REMOVE")
        self.assertEqual(clusters[0]["filter_reason"], "STRONG_SPECIES_NUMT_OVERLAP")
        self.assertTrue(all(v["filter_action"] == "REMOVE" for v in variants))
        self.assertTrue(all(v["strong_species_numt_evidence"] == "YES" for v in variants))

    def test_total_span_cap_blocks_overexpansion(self):
        clusters = copy.deepcopy(self.cluster_rows)
        variants = copy.deepcopy(self.variant_rows[:1])
        variants[0]["source_af"] = "0.205"
        diagnostics = apply_expansion(
            clusters, variants, max_delta_af=0.10, max_total_af_span=0.10
        )
        self.assertEqual(len(diagnostics), 1)
        self.assertEqual(diagnostics[0]["accepted"], "NO")
        self.assertEqual(variants[0]["filter_action"], "KEEP")


if __name__ == "__main__":
    unittest.main()
