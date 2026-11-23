#!/usr/bin/env python3
"""Unit checks for mean-mtDNA / median-nuclear mtCN and boundary semantics."""
import unittest

from qc_analysis.lib.mtcn_qc import mtcn_from_coverage, mtcn_from_summary
from qc_analysis.scripts.run_sample_variant_filtering import evaluate


class MeanMtcnQC(unittest.TestCase):
    def test_formula(self):
        result, denominator, adjusted = mtcn_from_coverage(100, 20)
        self.assertEqual(result, 10.0)
        self.assertEqual(denominator, 20)
        self.assertFalse(adjusted)

    def test_nuclear_zero_preserves_existing_convention(self):
        result, denominator, adjusted = mtcn_from_coverage(100, 0)
        self.assertEqual(result, 200.0)
        self.assertEqual(denominator, 1)
        self.assertTrue(adjusted)

    def test_negative_or_missing_coverage_does_not_pass(self):
        self.assertIsNone(mtcn_from_coverage(None, 5)[0])
        self.assertIsNone(mtcn_from_coverage(100, -5)[0])

    def test_summary_prefers_new_formula(self):
        mtcn, _, _ = mtcn_from_summary({
            "mt_mean_coverage": "80",
            "nuclear_median_coverage": "8",
            "mtcn_median": "1",
        })
        self.assertEqual(mtcn, 20)

    def test_hetero_percent_strict_nuclear_inclusive(self):
        row = {
            "sample": "example", "species": "test",
            "mt_mean_coverage": "100",
            "nuclear_median_coverage": "5",
            "mtcn_median": "1",
            "Percent_100": "90",
        }
        thresholds = {
            "mt_median_coverage_min": 100,
            "percent_100_min": 90,
            "nuclear_median_coverage_min": 5,
            "mtcn_min": 40,
            "mad_max": 0.5,
        }
        kwargs = dict(
            thresholds=thresholds, mt_coverage_enabled=False,
            mad_enabled=False, nuclear_coverage_enabled=True,
            mtcn_enabled=True, strict_greater_than=False,
            percent_strict_greater_than=True,
        )
        at90 = evaluate(row, **kwargs)
        self.assertEqual(at90["mtcn_qc_value"], 40.0)
        self.assertTrue(at90["pass_nuclear_coverage"])
        self.assertTrue(at90["pass_mtcn"])
        self.assertFalse(at90["pass_percent_coverage"])
        row["Percent_100"] = "90.01"
        above90 = evaluate(row, **kwargs)
        self.assertEqual(above90["qc_status"], "PASS")


if __name__ == "__main__":
    unittest.main()
