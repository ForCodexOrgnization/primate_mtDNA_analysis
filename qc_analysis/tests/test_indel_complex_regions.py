import unittest

from qc_analysis.scripts.detect_indel_complex_regions import (
    build_regions,
    circular_distance,
    detect_candidate_windows,
    merge_candidate_windows,
    parse_contig_header,
)


class IndelComplexRegionTests(unittest.TestCase):
    def test_af_matched_region_does_not_require_preexisting_cluster(self):
        residuals = [
            {"variant_index": 0, "sample": "S1", "chrom": "chrM", "pos": 1000, "af": 0.20, "clustered": "NO"},
            {"variant_index": 1, "sample": "S1", "chrom": "chrM", "pos": 1020, "af": 0.22, "clustered": "NO"},
            {"variant_index": 2, "sample": "S1", "chrom": "chrM", "pos": 1040, "af": 0.24, "clustered": "NO"},
        ]
        indels = [
            {"sample": "S1", "chrom": "chrM", "pos": 1025, "af": 0.21, "variant_type": "INDEL"},
        ]
        hits = detect_candidate_windows(residuals, indels, [], mt_length=17000, radius_bp=100)
        self.assertEqual(len(hits), 1)
        self.assertIn("AF_MATCHED", hits[0]["patterns"])

        regions = build_regions(
            merge_candidate_windows(hits),
            {"S1": 17000},
            {"S1": "VCF_HEADER_SOURCE_CHROM"},
            {"S1": "chrM"},
            100,
        )
        self.assertEqual(len(regions), 1)
        self.assertFalse(regions[0]["high_confidence"])
        self.assertEqual(regions[0]["decision"], "FLAG_REGION_REMOVE_AF_MATCHED_VARIANTS")
        self.assertEqual(regions[0]["n_af_matched_het_variants"], 3)

    def test_two_independent_patterns_become_high_confidence(self):
        residuals = [
            {"variant_index": 0, "sample": "S1", "chrom": "chrM", "pos": 1000, "af": 0.20, "clustered": "NO"},
            {"variant_index": 1, "sample": "S1", "chrom": "chrM", "pos": 1020, "af": 0.22, "clustered": "NO"},
            {"variant_index": 2, "sample": "S1", "chrom": "chrM", "pos": 1040, "af": 0.24, "clustered": "NO"},
        ]
        indels = [
            {"sample": "S1", "chrom": "chrM", "pos": 1010, "af": 0.21, "variant_type": "INDEL"},
            {"sample": "S1", "chrom": "chrM", "pos": 1030, "af": 0.23, "variant_type": "INDEL"},
        ]
        hits = detect_candidate_windows(residuals, indels, [], mt_length=17000, radius_bp=100)
        regions = build_regions(
            merge_candidate_windows(hits),
            {"S1": 17000},
            {"S1": "VCF_HEADER_SOURCE_CHROM"},
            {"S1": "chrM"},
            100,
        )
        self.assertEqual(len(regions), 1)
        self.assertTrue(regions[0]["high_confidence"])
        self.assertIn("AF_MATCHED", regions[0]["patterns"])
        self.assertIn("MULTI_INDEL", regions[0]["patterns"])
        self.assertEqual(regions[0]["decision"], "REMOVE_REGION_HIGH_CONF")

    def test_waterfall_single_signal_is_flag_only(self):
        residuals = [
            {"variant_index": 0, "sample": "S1", "chrom": "chrM", "pos": 16000, "af": 0.15, "clustered": "NO"},
            {"variant_index": 1, "sample": "S1", "chrom": "chrM", "pos": 16020, "af": 0.35, "clustered": "NO"},
            {"variant_index": 2, "sample": "S1", "chrom": "chrM", "pos": 16040, "af": 0.70, "clustered": "NO"},
        ]
        indels = [
            {"sample": "S1", "chrom": "chrM", "pos": 16030, "af": 0.82, "variant_type": "INDEL"},
        ]
        low_snvs = [
            {"sample": "S1", "chrom": "chrM", "pos": 16010, "af": 0.03, "variant_type": "LOW_AF_SNV"},
            {"sample": "S1", "chrom": "chrM", "pos": 16050, "af": 0.06, "variant_type": "LOW_AF_SNV"},
        ]
        hits = detect_candidate_windows(
            residuals,
            indels,
            low_snvs,
            mt_length=17000,
            radius_bp=100,
            af_match_max_delta=0.01,
            multiple_indel_min=2,
        )
        self.assertEqual(hits[0]["patterns"], ["WATERFALL"])
        regions = build_regions(
            merge_candidate_windows(hits),
            {"S1": 17000},
            {"S1": "VCF_HEADER_SOURCE_CHROM"},
            {"S1": "chrM"},
            100,
            af_match_max_delta=0.01,
            multiple_indel_min=2,
        )
        self.assertFalse(regions[0]["high_confidence"])
        self.assertEqual(regions[0]["decision"], "FLAG_REGION_SINGLE_SIGNAL")

    def test_native_length_avoids_human_length_bug(self):
        self.assertEqual(circular_distance(17003, 403, 17100), 500)
        self.assertEqual(circular_distance(17003, 403, 16569), 16600)
        self.assertEqual(circular_distance(17003, 403, None), 16600)

    def test_parse_contig_header(self):
        contig, length = parse_contig_header("##contig=<ID=chrM,length=17123>")
        self.assertEqual(contig, "chrM")
        self.assertEqual(length, 17123)

    def test_sparse_context_not_detected(self):
        residuals = [
            {"variant_index": 0, "sample": "S1", "chrom": "chrM", "pos": 1000, "af": 0.20, "clustered": "NO"},
            {"variant_index": 1, "sample": "S1", "chrom": "chrM", "pos": 1010, "af": 0.30, "clustered": "NO"},
        ]
        indels = [{"sample": "S1", "chrom": "chrM", "pos": 1005, "af": 0.25, "variant_type": "INDEL"}]
        hits = detect_candidate_windows(residuals, indels, [], mt_length=17000, radius_bp=100, min_residual_het=3)
        self.assertEqual(hits, [])


if __name__ == "__main__":
    unittest.main()
