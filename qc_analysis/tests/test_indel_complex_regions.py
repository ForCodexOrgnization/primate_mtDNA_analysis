import unittest

from qc_analysis.scripts.detect_indel_complex_regions import detect_candidate_windows


class IndelComplexRegionTests(unittest.TestCase):
    def test_af_matched_region_does_not_require_preexisting_cluster(self):
        residuals = [
            {"variant_index": 0, "sample": "S1", "pos": 1000, "af": 0.20, "clustered": "NO"},
            {"variant_index": 1, "sample": "S1", "pos": 1020, "af": 0.22, "clustered": "NO"},
            {"variant_index": 2, "sample": "S1", "pos": 1040, "af": 0.24, "clustered": "NO"},
        ]
        indels = [
            {"sample": "S1", "pos": 1025, "af": 0.21, "variant_type": "INDEL"},
        ]
        hits = detect_candidate_windows(residuals, indels, [], radius_bp=100)
        self.assertEqual(len(hits), 1)
        self.assertIn("AF_MATCHED", hits[0]["patterns"])

    def test_waterfall_pattern_detected_with_wide_af_span(self):
        residuals = [
            {"variant_index": 0, "sample": "S1", "pos": 16000, "af": 0.15, "clustered": "NO"},
            {"variant_index": 1, "sample": "S1", "pos": 16020, "af": 0.35, "clustered": "NO"},
            {"variant_index": 2, "sample": "S1", "pos": 16040, "af": 0.70, "clustered": "NO"},
        ]
        indels = [
            {"sample": "S1", "pos": 16030, "af": 0.82, "variant_type": "INDEL"},
        ]
        low_snvs = [
            {"sample": "S1", "pos": 16010, "af": 0.03, "variant_type": "LOW_AF_SNV"},
            {"sample": "S1", "pos": 16050, "af": 0.06, "variant_type": "LOW_AF_SNV"},
        ]
        hits = detect_candidate_windows(
            residuals,
            indels,
            low_snvs,
            radius_bp=100,
            af_match_max_delta=0.01,
            multiple_indel_min=2,
            waterfall_min_low_snv=2,
            waterfall_min_total_snv=5,
            waterfall_min_snv_af_span=0.10,
        )
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0]["patterns"], ["WATERFALL"])

    def test_sparse_context_not_detected(self):
        residuals = [
            {"variant_index": 0, "sample": "S1", "pos": 1000, "af": 0.20, "clustered": "NO"},
            {"variant_index": 1, "sample": "S1", "pos": 1010, "af": 0.30, "clustered": "NO"},
        ]
        indels = [{"sample": "S1", "pos": 1005, "af": 0.25, "variant_type": "INDEL"}]
        hits = detect_candidate_windows(residuals, indels, [], radius_bp=100, min_residual_het=3)
        self.assertEqual(hits, [])


if __name__ == "__main__":
    unittest.main()
