import gzip
import tempfile
import unittest
from pathlib import Path

from qc_analysis.scripts.run_coordinate_liftover import Sample, SampleStats, lift_vcf


def mapped(pos):
    return {pos: {"map_status": "mapped", "human_pos_canonical": str(pos)}}


class CoordinateLiftoverTests(unittest.TestCase):
    def run_one(self, body, human="G", pos_map=None, fail=False):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            vcf = d / "in.vcf.gz"
            text = """##fileformat=VCFv4.2
##INFO=<ID=R_FIELD,Number=R,Type=Integer,Description="R">
##INFO=<ID=A_FIELD,Number=A,Type=Float,Description="A">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="AF">
##FORMAT=<ID=F1R2,Number=R,Type=Integer,Description="F1R2">
##FORMAT=<ID=F2R1,Number=R,Type=Integer,Description="F2R1">
##FORMAT=<ID=SB,Number=4,Type=Integer,Description="SB">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
""" + body
            with gzip.open(vcf, "wt") as h:
                h.write(text)
            sample = Sample("S", Path("x.fa"), vcf, Path("x.cov"))
            stats = SampleStats("S")
            out = d / "out.vcf"
            unresolved = d / "unresolved.tsv"
            lift_vcf(sample, pos_map or mapped(1), out, human, "chrM", True, fail, stats, unresolved, True)
            return out.read_text(), unresolved.read_text(), stats

    def records(self, text):
        return [l for l in text.splitlines() if l and not l.startswith("#")]

    def test_ref_match_unchanged(self):
        out, unr, stats = self.run_one("chrS\t1\t.\tG\tA\t.\tPASS\t.\tGT:AD:AF\t0/1:90,10:0.1\n", "G")
        rec = self.records(out)[0].split("\t")
        self.assertEqual(rec[3:5], ["G", "A"])
        self.assertEqual(rec[9], "0/1:90,10:0.1")
        self.assertIn("LIFTOVER_ALLELE_STATUS=REF_MATCH", rec[7])
        self.assertEqual(stats.ref_match_count, 1)

    def test_alt_ref_flip_fields(self):
        out, unr, stats = self.run_one("chrS\t1\t.\tG\tA\t.\tPASS\t.\tGT:AD:AF:F1R2:F2R1:SB\t0/1:90,10:0.1:80,8:10,2:50,40,6,4\n", "A")
        rec = self.records(out)[0].split("\t")
        self.assertEqual(rec[3:5], ["A", "G"])
        self.assertEqual(rec[9], "1/0:10,90:0.9:8,80:2,10:6,4,50,40")
        self.assertIn("LIFTOVER_ALLELE_STATUS=ALT_REF_FLIP", rec[7])
        self.assertEqual(stats.alt_ref_flip_count, 1)

    def test_af_fallback_without_ad(self):
        out, _, _ = self.run_one("chrS\t1\t.\tG\tA\t.\tPASS\t.\tGT:AF\t0/1:0.1\n", "A")
        self.assertIn("1/0:0.9", self.records(out)[0])

    def test_multiallelic_unresolved(self):
        out, unr, stats = self.run_one("chrS\t1\t.\tG\tA,C\t.\tPASS\t.\tGT\t0/1\n", "A")
        self.assertEqual(self.records(out), [])
        self.assertIn("MULTIALLELIC", unr)

    def test_indel_flip_candidate_unsupported(self):
        out, unr, stats = self.run_one("chrS\t1\t.\tC\tGA\t.\tPASS\t.\tGT\t0/1\n", "GA")
        self.assertEqual(self.records(out), [])
        self.assertIn("FLIP_UNSUPPORTED_VARIANT_TYPE", unr)
        self.assertEqual(stats.unsupported_flip_count, 1)

    def test_target_ref_neither_continue_and_fail(self):
        out, unr, stats = self.run_one("chrS\t1\t.\tG\tA\t.\tPASS\t.\tGT\t0/1\n", "C", fail=False)
        self.assertIn("TARGET_REF_NOT_SOURCE_REF_OR_ALT", unr)
        with self.assertRaises(ValueError):
            self.run_one("chrS\t1\t.\tG\tA\t.\tPASS\t.\tGT\t0/1\n", "C", fail=True)

    def test_number_r_and_number_a_info(self):
        out, _, _ = self.run_one("chrS\t1\t.\tG\tA\t.\tPASS\tR_FIELD=90,10;A_FIELD=5\tGT\t0/1\n", "A")
        info = self.records(out)[0].split("\t")[7]
        self.assertIn("R_FIELD=10,90", info)
        self.assertNotIn("A_FIELD=5", info)
        self.assertIn("LIFTOVER_DROPPED_INFO_FIELDS=A_FIELD", info)

    def test_final_ref_check_blocks_bad_ref(self):
        out, unr, stats = self.run_one("chrS\t1\t.\tGG\tA\t.\tPASS\t.\tGT\t0/1\n", "G")
        self.assertEqual(self.records(out), [])
        self.assertIn("TARGET_REF_NOT_SOURCE_REF_OR_ALT", unr)

    def test_stats_mixed(self):
        body = "".join([
            "chrS\t1\t.\tG\tA\t.\tPASS\t.\tGT\t0/1\n",
            "chrS\t2\t.\tG\tA\t.\tPASS\t.\tGT:AD:AF\t0/1:90,10:0.1\n",
            "chrS\t3\t.\tG\tA,C\t.\tPASS\t.\tGT\t0/1\n",
            "chrS\t4\t.\tG\tA\t.\tPASS\t.\tGT\t0/1\n",
        ])
        out, unr, stats = self.run_one(body, "GACC", pos_map=mapped(1)|mapped(2)|mapped(3)|mapped(4))
        self.assertEqual(stats.vcf_variants_lifted, 2)
        self.assertEqual(stats.vcf_variants_failed_liftover, 2)
        self.assertEqual(stats.ref_match_count, 1)
        self.assertEqual(stats.ref_mismatch_count, 2)
        self.assertEqual(stats.alt_ref_flip_count, 1)
        self.assertEqual(stats.unresolved_ref_mismatch_count, 1)
        self.assertEqual(stats.unsupported_flip_count, 0)


if __name__ == "__main__":
    unittest.main()
