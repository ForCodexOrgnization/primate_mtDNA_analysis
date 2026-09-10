import unittest

from qc_analysis.scripts.propagate_species_numt_positions import (
    apply_species_numt_position_propagation,
    build_species_numt_position_catalogue,
)


class SpeciesNumtPositionPropagationTests(unittest.TestCase):
    def setUp(self):
        self.cluster_rows = [
            {
                "sample": "A",
                "species": "Papio_anubis",
                "reference_key": "Papio_anubis",
                "cluster_id": "A_C1",
                "any_numt_evidence": "YES",
                "numt_tier": "HIGH_CONF",
            },
            {
                "sample": "C",
                "species": "Other_species",
                "reference_key": "Other_species",
                "cluster_id": "C_C1",
                "any_numt_evidence": "NO",
                "numt_tier": "NONE",
            },
        ]
        self.variant_rows = [
            {
                "sample": "A",
                "species": "Papio_anubis",
                "source_chrom": "chrM",
                "source_pos": "1000",
                "source_ref": "A",
                "source_alt": "G",
                "source_af": "0.12",
                "source_dp": "500",
                "clustered": "YES",
                "cluster_id": "A_C1",
                "cluster_class": "NUMT_ONLY",
                "numt_scope": "SAMPLE",
                "numt_tier": "HIGH_CONF",
                "recurrence": "NO",
                "filter_action": "REMOVE",
                "filter_reason": "SAMPLE_NUMT_OVERLAP",
            },
            {
                "sample": "B",
                "species": "Papio_anubis",
                "source_chrom": "chrM",
                "source_pos": "1000",
                "source_ref": "A",
                "source_alt": "T",
                "source_af": "0.18",
                "source_dp": "700",
                "clustered": "NO",
                "cluster_id": "",
                "cluster_class": "",
                "numt_scope": "NONE",
                "numt_tier": "NONE",
                "recurrence": "NO",
                "filter_action": "KEEP",
                "filter_reason": "",
            },
            {
                "sample": "B",
                "species": "Papio_anubis",
                "source_chrom": "chrM",
                "source_pos": "1100",
                "source_ref": "C",
                "source_alt": "T",
                "source_af": "0.15",
                "source_dp": "700",
                "clustered": "NO",
                "cluster_id": "",
                "cluster_class": "",
                "numt_scope": "NONE",
                "numt_tier": "NONE",
                "recurrence": "NO",
                "filter_action": "KEEP",
                "filter_reason": "",
            },
        ]

    def test_propagates_by_position_not_ref_alt(self):
        catalogue = build_species_numt_position_catalogue(self.cluster_rows, self.variant_rows)
        out = apply_species_numt_position_propagation(self.cluster_rows, self.variant_rows, catalogue)
        propagated = out[1]
        self.assertEqual(propagated["numt_propagated"], "YES")
        self.assertEqual(propagated["filter_action"], "REMOVE")
        self.assertEqual(propagated["filter_reason"], "SAME_SPECIES_NUMT_POSITION")
        self.assertEqual(propagated["numt_source_samples"], "A")

    def test_does_not_propagate_unobserved_position(self):
        catalogue = build_species_numt_position_catalogue(self.cluster_rows, self.variant_rows)
        out = apply_species_numt_position_propagation(self.cluster_rows, self.variant_rows, catalogue)
        self.assertEqual(out[2]["numt_propagated"], "NO")
        self.assertEqual(out[2]["filter_action"], "KEEP")


if __name__ == "__main__":
    unittest.main()
