import tempfile
import unittest
from pathlib import Path

from qc_analysis.lib.match_utils import (
    compare_values, lift_source_pos_to_human, load_coordinate_map,
    is_iupac_dna_base, is_iupac_rna_base, is_resolved_rna_base,
    normalize_rna_base, normalize_rna_symbol, orient_dna_base_to_rna,
    pair_effect, pair_type, rna_symbols_compatible,
)
from qc_analysis.scripts.run_rrna_match import (
    infer_species_pair_pos_from_human_pair_local, load_rrna_structure_table,
    normalize_rrna_gene,
)


class RnaMatchUtilityTests(unittest.TestCase):
    def test_rna_pair_types_and_effects_normalize_dna(self):
        self.assertEqual(normalize_rna_base("t"), "U")
        self.assertEqual(pair_type("A", "T"), "WC")
        self.assertEqual(pair_type("G", "U"), "GU_wobble")
        self.assertEqual(pair_type("A", "C"), "non_WC")
        self.assertEqual(pair_type("N", "C"), "ambiguous")
        self.assertEqual(pair_effect("WC", "non_WC"), "WC_to_non_WC")
        self.assertEqual(pair_effect("WC", "WC"), "unchanged")

    def test_iupac_symbols_are_preserved_but_not_treated_as_resolved(self):
        for symbol in "ACGTRYSWKMBDHVN":
            self.assertTrue(is_iupac_dna_base(symbol))
        for symbol in "ACGURYSWKMBDHVN":
            self.assertTrue(is_iupac_rna_base(symbol))
        self.assertEqual(normalize_rna_symbol("m"), "M")
        self.assertIsNone(normalize_rna_base("M"))
        self.assertFalse(is_resolved_rna_base("R"))

    def test_iupac_orientation_compatibility_and_conservative_pairs(self):
        self.assertEqual(orient_dna_base_to_rna("R", "+"), "R")
        expected={"R":"Y","Y":"R","M":"K","K":"M","W":"W","S":"S",
                  "B":"V","V":"B","D":"H","H":"D","A":"U","T":"A"}
        for genomic,rna in expected.items():
            self.assertEqual(orient_dna_base_to_rna(genomic, "-"), rna)
        for reference,allowed in {"R":"AGR","Y":"CUY","M":"ACM","W":"AUW"}.items():
            for observed in allowed:
                self.assertTrue(rna_symbols_compatible(reference, observed))
        for observed in "CU": self.assertFalse(rna_symbols_compatible("R", observed))
        self.assertEqual(pair_type("R", "C"), "ambiguous")
        self.assertEqual(pair_effect("ambiguous", "WC"), "NA")
        self.assertEqual(pair_effect("WC", "ambiguous"), "NA")

    def test_rrna_gene_normalization_and_pair_coordinate_inference(self):
        self.assertEqual(normalize_rrna_gene("12S"), "MT-RNR1")
        self.assertEqual(normalize_rrna_gene("rnr2"), "MT-RNR2")
        self.assertEqual(infer_species_pair_pos_from_human_pair_local({"start": "100", "end": "200", "strand": "+"}, "4"), "103")
        self.assertEqual(infer_species_pair_pos_from_human_pair_local({"start": "100", "end": "200", "strand": "-"}, "4"), "197")

    def test_structure_table_requires_columns_and_coordinate_map_lifts(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            bad = d / "bad.tsv"; bad.write_text("rrna_gene\thuman_pos\nMT-RNR1\t1\n")
            with self.assertRaises(ValueError):
                load_rrna_structure_table(bad)
            mp = d / "map.tsv"
            mp.write_text("species_pos_original\thuman_pos_canonical\n103\t203\n")
            self.assertEqual(lift_source_pos_to_human(103, load_coordinate_map(mp)), "203")
            self.assertEqual(compare_values("203", "203"), "yes")
            self.assertEqual(compare_values(".", "203"), ".")
