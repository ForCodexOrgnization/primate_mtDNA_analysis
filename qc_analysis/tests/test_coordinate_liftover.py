import configparser
import gzip
import tempfile
import unittest
from pathlib import Path

from qc_analysis.scripts.run_coordinate_liftover import (
    Sample, SampleStats, _validate_input_file, build_map, find_sample_file, lift_vcf,
    load_anchor_positions, read_fasta, restore_human_pos, sample_from_row,
)
from qc_analysis.lib.mt_anchor_utils import IUPAC_DNA, mask_ambiguity_for_alignment, sequence_sha256


def mapped(pos):
    return {pos: {"map_status": "mapped", "human_pos_canonical": str(pos)}}


class CoordinateLiftoverTests(unittest.TestCase):
    def run_one(self, body, human="G", pos_map=None, fail=False, source_reference_seq=""):
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
            lift_vcf(sample, pos_map or mapped(1), out, human, "chrM", True, fail, stats, unresolved, True, source_reference_seq)
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

    def test_ambiguous_source_reference_variant_is_unresolved_but_other_record_lifts(self):
        body = "chrS\t1\t.\tY\tA\t.\tPASS\t.\tGT\t0/1\nchrS\t2\t.\tG\tA\t.\tPASS\t.\tGT\t0/1\n"
        out, unresolved, stats = self.run_one(body, "GG", mapped(1) | mapped(2), source_reference_seq="YG")
        self.assertEqual(len(self.records(out)), 1)
        self.assertIn("SOURCE_REFERENCE_AMBIGUOUS", unresolved)
        self.assertIn("\tY\tSOURCE_REFERENCE_AMBIGUOUS", unresolved)
        self.assertEqual(stats.variants_overlapping_ambiguous_reference, 1)


class HumanCoordinateRestorationTests(unittest.TestCase):
    def test_rotated_human_boundary_positions_restore_from_anchor(self):
        self.assertEqual(restore_human_pos(1, 3059, 16568, 0), 3059)
        self.assertEqual(restore_human_pos(2, 3059, 16568, 0), 3060)
        self.assertEqual(restore_human_pos(13511, 3059, 16568, 0), 1)

    def test_rotated_human_alignment_matches_unrotated_canonical_map(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            sample = Sample("synthetic", Path("species.fa"), Path("in.vcf"), Path("in.cov"))
            canonical = "ACGTACGT"
            # The pairwise alignment is deliberately ungapped: this isolates
            # coordinate restoration from alignment scoring.
            unrotated = d / "unrotated.fa"
            unrotated.write_text(f">species\n{canonical}\n>human\n{canonical}\n")
            rotated = d / "rotated.fa"
            rotated_human = canonical[2:] + canonical[:2]
            rotated.write_text(f">species\n{rotated_human}\n>human\n{rotated_human}\n")
            canonical_map, _ = build_map(sample, unrotated, d / "canonical.tsv", 1, 1, len(canonical), 0)
            rotated_map, _ = build_map(sample, rotated, d / "rotated.tsv", 3, 3, len(canonical), 0)
            self.assertEqual(
                {pos: row["human_pos_canonical"] for pos, row in canonical_map.items()},
                {pos: row["human_pos_canonical"] for pos, row in rotated_map.items()},
            )


class IupacFastaTests(unittest.TestCase):
    def test_all_standard_iupac_codes_are_accepted_and_masked(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "iupac.fa"; path.write_text(">chrM\nACGTRYSWKMBDHVN\n")
            rec = read_fasta(path)
        self.assertEqual(set(rec.seq), IUPAC_DNA)
        self.assertEqual(mask_ambiguity_for_alignment(rec.seq), "ACGT" + "N" * 11)
        self.assertEqual(len(mask_ambiguity_for_alignment(rec.seq)), len(rec.seq))

    def test_y_and_r_are_accepted_and_hash_original_sequence(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "Cercopithecus_mitis.fa"; path.write_text(">chrM\nacgtyr\n")
            rec = read_fasta(path)
        self.assertEqual(rec.seq, "ACGTYR")
        self.assertNotEqual(sequence_sha256(rec.seq), sequence_sha256(mask_ambiguity_for_alignment(rec.seq)))

    def test_non_iupac_fasta_base_is_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "bad.fa"; path.write_text(">chrM\nACGTX\n")
            with self.assertRaisesRegex(ValueError, "X"):
                read_fasta(path)


class CoordinateLiftoverInputTests(unittest.TestCase):
    sample = "ERS12091931"
    vcf_stem = "ERS12091931.round2.original_coords.clean.final.split"
    cov_name = "ERS12091931.merged.max_depth.per_base_coverage.tsv"

    def config(self, directory: Path) -> configparser.ConfigParser:
        cfg = configparser.ConfigParser()
        cfg["paths"] = {
            "vcf_dir": str(directory / "vcf"),
            "vcf_pattern": "{sample}.round2.original_coords.clean.final.split.vcf.gz,{sample}.round2.original_coords.clean.final.split.vcf",
            "cov_dir": str(directory / "cov"),
            "cov_pattern": "{sample}.merged.max_depth.per_base_coverage.tsv",
            "species_fasta_dir": str(directory / "fasta"),
        }
        return cfg

    def touch(self, path: Path, text="x") -> Path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return path

    def vcf(self, directory: Path, suffix: str) -> Path:
        return self.touch(directory / "vcf" / f"{self.vcf_stem}{suffix}")

    def test_only_compressed_vcf_is_selected(self):
        with tempfile.TemporaryDirectory() as td:
            directory = Path(td); cfg = self.config(directory)
            compressed = self.vcf(directory, ".vcf.gz")
            self.assertEqual(find_sample_file(self.sample, cfg, "vcf_dir", "vcf_pattern", "VCF"), compressed)

    def test_only_uncompressed_vcf_is_selected_and_accepted(self):
        with tempfile.TemporaryDirectory() as td:
            directory = Path(td); cfg = self.config(directory)
            uncompressed = self.vcf(directory, ".vcf")
            self.assertEqual(find_sample_file(self.sample, cfg, "vcf_dir", "vcf_pattern", "VCF"), uncompressed)
            self.assertEqual(_validate_input_file(uncompressed, "VCF"), [])

    def test_compressed_vcf_has_priority_when_both_exist(self):
        with tempfile.TemporaryDirectory() as td:
            directory = Path(td); cfg = self.config(directory)
            compressed = self.vcf(directory, ".vcf.gz")
            self.vcf(directory, ".vcf")
            self.assertEqual(find_sample_file(self.sample, cfg, "vcf_dir", "vcf_pattern", "VCF"), compressed)

    def test_broken_compressed_symlink_uses_uncompressed_fallback(self):
        with tempfile.TemporaryDirectory() as td:
            directory = Path(td); cfg = self.config(directory)
            broken = directory / "vcf" / f"{self.vcf_stem}.vcf.gz"
            broken.parent.mkdir(parents=True)
            broken.symlink_to(directory / "missing.vcf.gz")
            uncompressed = self.vcf(directory, ".vcf")
            selected = find_sample_file(self.sample, cfg, "vcf_dir", "vcf_pattern", "VCF")
            self.assertEqual(selected, uncompressed)
            self.assertEqual(_validate_input_file(selected, "VCF"), [])

    def test_broken_compressed_symlink_without_fallback_is_diagnostic(self):
        with tempfile.TemporaryDirectory() as td:
            directory = Path(td); cfg = self.config(directory)
            broken = directory / "vcf" / f"{self.vcf_stem}.vcf.gz"
            broken.parent.mkdir(parents=True)
            target = directory / "old" / "missing.vcf.gz"
            broken.symlink_to(target)
            with self.assertRaisesRegex(FileNotFoundError, r"BROKEN_SYMLINK: .*\.vcf\.gz -> .*missing\.vcf\.gz") as raised:
                find_sample_file(self.sample, cfg, "vcf_dir", "vcf_pattern", "VCF")
            self.assertIn(f"{self.vcf_stem}.vcf", str(raised.exception))

    def test_valid_compressed_symlink_is_accepted(self):
        with tempfile.TemporaryDirectory() as td:
            directory = Path(td); cfg = self.config(directory)
            target = self.touch(directory / "target.vcf.gz")
            link = directory / "vcf" / f"{self.vcf_stem}.vcf.gz"
            link.parent.mkdir(parents=True)
            link.symlink_to(target)
            self.assertEqual(find_sample_file(self.sample, cfg, "vcf_dir", "vcf_pattern", "VCF"), link)

    def test_only_merged_coverage_is_selected(self):
        with tempfile.TemporaryDirectory() as td:
            directory = Path(td); cfg = self.config(directory)
            merged = self.touch(directory / "cov" / self.cov_name)
            self.touch(directory / "cov" / f"{self.sample}.round2.original_coords.per_base_coverage.tsv")
            self.assertEqual(find_sample_file(self.sample, cfg, "cov_dir", "cov_pattern", "COV"), merged)

    def test_original_coverage_is_not_a_fallback(self):
        with tempfile.TemporaryDirectory() as td:
            directory = Path(td); cfg = self.config(directory)
            self.touch(directory / "cov" / f"{self.sample}.round2.original_coords.per_base_coverage.tsv")
            with self.assertRaisesRegex(FileNotFoundError, "No merged coverage file found"):
                find_sample_file(self.sample, cfg, "cov_dir", "cov_pattern", "COV")

    def test_missing_vcf_reports_both_candidates(self):
        with tempfile.TemporaryDirectory() as td:
            directory = Path(td); cfg = self.config(directory)
            with self.assertRaisesRegex(FileNotFoundError, r"\.vcf\.gz.*\.vcf"):
                find_sample_file(self.sample, cfg, "vcf_dir", "vcf_pattern", "VCF")

    def test_sample_loading_accepts_uncompressed_vcf_and_merged_coverage(self):
        with tempfile.TemporaryDirectory() as td:
            directory = Path(td); cfg = self.config(directory)
            self.touch(directory / "fasta" / "species.fa", ">chrM\nACGT\n")
            vcf = self.vcf(directory, ".vcf")
            cov = self.touch(directory / "cov" / self.cov_name)
            sample = sample_from_row({"sample": self.sample, "species": "species"}, cfg)
            self.assertEqual(sample.vcf, vcf)
            self.assertEqual(sample.cov, cov)
            self.assertEqual(sample.missing_files, [])


if __name__ == "__main__":
    unittest.main()

class AnchorValidationTests(unittest.TestCase):
    def test_sequence_sha256_alias_reuses_validated_anchor(self):
        from qc_analysis.scripts.run_coordinate_liftover import select_runtime_anchor, read_workflow_config
        with tempfile.TemporaryDirectory() as td:
            d = Path(td); cfgp = d / 'c.yaml'; af = d / 'anchors.tsv'
            seq = 'ACGT'
            sha = sequence_sha256(seq)
            af.write_text(
                'reference_id\tsequence_sha256\tsequence_length\tanchor_original_position\thuman_anchor_original_position\tanchor_method\tanchor_qc_status\n'
                f'ref_primary\t{sha}\t4\t2\t3\tGLOBAL_MSA_ANCHOR\tPASS\n'
            )
            cfgp.write_text(
                f'coordinate_liftover:\n  coordinates:\n    anchor_positions_file: {af}\n'
                '  anchor:\n    require_validated_anchor: true\n    verify_sequence_sha256: true\n'
                '    allow_pairwise_anchor_fallback: false\n    allow_anchor_position_one_fallback: false\n'
            )
            cfg = read_workflow_config(cfgp)
            by_ref, by_sha = load_anchor_positions(cfg)
            for ref_id, expected_method in [('ref_primary', 'REFERENCE_ID'), ('ref_alias', 'SEQUENCE_SHA256_ALIAS')]:
                sample = Sample('S', Path('ref.fa'), Path('x.vcf.gz'), Path('x.cov'), reference_id=ref_id)
                got = select_runtime_anchor(sample, seq, seq, cfg, by_ref, '', by_sha)
                self.assertEqual(got[0:2], (2, 3))
                self.assertEqual(got[-1], expected_method)

    def test_sequence_hash_mismatch_rejected(self):
        import tempfile, csv
        from qc_analysis.scripts.run_coordinate_liftover import select_runtime_anchor, read_workflow_config
        with tempfile.TemporaryDirectory() as td:
            d=Path(td); cfgp=d/'c.yaml'; af=d/'anchors.tsv'
            af.write_text('reference_id\tsequence_sha256\tsequence_length\tanchor_original_position\thuman_anchor_original_position\tanchor_method\tanchor_qc_status\nref1\tbad\t4\t2\t1\tGLOBAL_MSA_ANCHOR\tPASS\n')
            cfgp.write_text(f'coordinate_liftover:\n  coordinates:\n    anchor_positions_file: {af}\n  anchor:\n    require_validated_anchor: true\n    verify_sequence_sha256: true\n    allow_pairwise_anchor_fallback: false\n    allow_anchor_position_one_fallback: false\n')
            cfg=read_workflow_config(cfgp); anchors={'ref1': next(csv.DictReader(af.open(), delimiter='\t'))}
            s=Sample('S', Path('ref.fa'), Path('x.vcf.gz'), Path('x.cov'), reference_id='ref1')
            with self.assertRaisesRegex(ValueError, 'ANCHOR_REFERENCE_HASH_MISMATCH'):
                select_runtime_anchor(s, 'ACGT', 'ACGT', cfg, anchors, '')

    def test_sample_override_backward_compatible(self):
        from qc_analysis.scripts.run_coordinate_liftover import select_runtime_anchor, read_workflow_config
        with tempfile.TemporaryDirectory() as td:
            cfgp=Path(td)/'c.yaml'; cfgp.write_text('coordinate_liftover:\n  anchor:\n    require_validated_anchor: true\n    allow_pairwise_anchor_fallback: false\n')
            cfg=read_workflow_config(cfgp)
            s=Sample('S', Path('ref.fa'), Path('x.vcf.gz'), Path('x.cov'), rotate_anchor=3)
            got=select_runtime_anchor(s, 'ACGT', 'ACGT', cfg, {}, '')
            self.assertEqual(got[0],3); self.assertEqual(got[2], 'SAMPLE_OVERRIDE')
