import tempfile
import unittest
from pathlib import Path

from qc_analysis.lib.match_utils import (
    compare_values, lift_source_pos_to_human, load_coordinate_map,
    normalize_rna_base, pair_effect, pair_type,
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
        self.assertEqual(pair_type("N", "C"), "NA")
        self.assertEqual(pair_effect("WC", "non_WC"), "WC_to_non_WC")
        self.assertEqual(pair_effect("WC", "WC"), "unchanged")

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
