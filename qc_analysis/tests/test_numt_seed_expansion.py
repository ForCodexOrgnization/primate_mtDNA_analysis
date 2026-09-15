import copy
import unittest

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
                "recurrence": "NO",
                "numt_chrm_start": "1941",
                "numt_chrm_end": "2115",
                "numt_tier": "BESTHIT_ONLY_NUMT",
                "cluster_class": "NUMT_ONLY",
                "filter_action": "REMOVE",
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

    def test_species_numt_without_recurrence_cannot_drive_expansion(self):
        clusters = copy.deepcopy(self.cluster_rows)
        clusters[0]["sample_numt_evidence"] = "NO"
        clusters[0]["species_numt_evidence"] = "YES"
        clusters[0]["recurrence"] = "NO"
        variants = copy.deepcopy(self.variant_rows)
        diagnostics = apply_expansion(clusters, variants)
        self.assertEqual(diagnostics, [])
        self.assertTrue(all(v["filter_action"] == "KEEP" for v in variants))

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
