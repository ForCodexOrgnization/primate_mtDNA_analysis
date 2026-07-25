import csv
import gzip
import tempfile
import unittest
from pathlib import Path

from qc_analysis.lib.reference_utils import normalized_fasta_sequence_sha256
from qc_analysis.scripts.compare_genbank_mitos2_reference_annotations import compare


class CompareReferenceAnnotationsTests(unittest.TestCase):
    def rows(self, key='', sha='', species='', accession='A.1', end=3, base='ATG', **extra):
        rows = []
        for position in range(1, end + 1):
            row = dict(reference_key=key, coordinate_reference_accession=accession,
                       gene='MT-ND1', pos=str(position),
                       ref_base_genome=base[position - 1], strand='+', codon_index='1',
                       codon_pos_in_triplet=str(position), codon_seq=base,
                       target_species=species)
            if sha:
                row['genbank_record_sequence_sha256'] = sha
            row.update(extra)
            rows.append(row)
        return rows

    def write(self, path, rows):
        with path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=sorted({k for r in rows for k in r}), delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)

    def run_compare(self, directory, genbank_rows, mitos_rows, fail=True):
        g, m = directory / 'g.tsv', directory / 'm.tsv'
        self.write(g, genbank_rows)
        self.write(m, mitos_rows)
        paths = directory / 'genes.tsv', directory / 'summary.tsv', directory / 'diagnostics.tsv'
        compare(g, m, *paths, fail_no_shared=fail)
        result = []
        for path in paths:
            with path.open() as handle:
                result.append(list(csv.DictReader(handle, delimiter='\t')))
        return result

    def test_hash_pairing_does_not_require_shared_reference_keys(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            genes, summary, diagnostics = self.run_compare(
                d, self.rows('genbank-key', 'a' * 64),
                self.rows('mitos-key', species='Species one', mitos2_input_sequence_sha256='a' * 64))
            self.assertEqual(len(genes), 13)
            self.assertEqual(len(summary), 1)
            self.assertFalse(diagnostics)
            self.assertEqual(summary[0]['comparison_sequence_sha256'], 'a' * 64)

    def test_legacy_mitos_table_hashes_input_fasta_and_ignores_header(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            genbank_fasta = d / 'genbank.fa'
            genbank_fasta.write_text('>GenBank header\na t\ng\n')
            mitos_fasta = d / 'mitos.fa.gz'
            with gzip.open(mitos_fasta, 'wt') as handle:
                handle.write('>chrM completely different header\nATG\n')
            sha = normalized_fasta_sequence_sha256(genbank_fasta)['sequence_sha256']
            legacy = self.rows('', species='Legacy species', mitos2_input_fasta=str(mitos_fasta))
            genes, summary, _ = self.run_compare(d, self.rows('gb', sha), legacy)
            self.assertEqual(len(genes), 13)
            self.assertEqual(summary[0]['mitos2_reference_key'], 'Legacy species')
            self.assertEqual(summary[0]['mitos2_input_sequence_sha256'], sha)

    def test_legacy_falls_back_to_coordinate_fasta(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            fasta = d / 'coordinate.fasta'
            fasta.write_text('>old\nATG\n')
            sha = normalized_fasta_sequence_sha256(fasta)['sequence_sha256']
            _, summary, _ = self.run_compare(
                d, self.rows('gb', sha),
                self.rows('', species='Legacy', coordinate_reference_fasta=str(fasta)))
            self.assertEqual(summary[0]['comparison_sequence_sha256'], sha)

    def test_multiple_targets_are_compared_once_and_all_are_recorded(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            sha = 'b' * 64
            mitos = (self.rows('m1', species='Species one', mitos2_input_sequence_sha256=sha, accession='M1') +
                     self.rows('m2', species='Species two', mitos2_input_sequence_sha256=sha, accession='M2'))
            genes, summary, diagnostics = self.run_compare(d, self.rows('g', sha), mitos)
            self.assertEqual(len(genes), 13)
            self.assertEqual(len(summary), 1)
            self.assertFalse(diagnostics)
            self.assertEqual(summary[0]['sequence_usage_category'], 'shared_sequence_multiple_targets')
            self.assertEqual(summary[0]['mitos2_target_species'], 'Species one;Species two')
            self.assertEqual(summary[0]['mitos2_accessions'], 'M1;M2')
            self.assertEqual(summary[0]['n_mitos2_groups'], '2')

    def test_inconsistent_same_sequence_annotations_are_diagnostic(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            sha = 'c' * 64
            genbank = self.rows('g1', sha) + self.rows('g2', sha, base='ACG')
            genes, summary, diagnostics = self.run_compare(
                d, genbank,
                self.rows('m', species='Species', mitos2_input_sequence_sha256=sha), fail=False)
            self.assertFalse(genes)
            self.assertFalse(summary)
            self.assertEqual(diagnostics[0]['reason_comparison_skipped'],
                             'same_sequence_inconsistent_genbank_annotation')

    def test_inconsistent_mitos_annotations_are_diagnostic(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            sha = 'd' * 64
            mitos = (self.rows('m1', species='One', mitos2_input_sequence_sha256=sha) +
                     self.rows('m2', species='Two', mitos2_input_sequence_sha256=sha, base='ACG'))
            _, _, diagnostics = self.run_compare(d, self.rows('g', sha), mitos, fail=False)
            self.assertEqual(diagnostics[0]['sequence_compatibility_category'],
                             'same_sequence_inconsistent_mitos2_annotation')

    def test_unmatched_references_do_not_block_valid_comparison(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            genbank = self.rows('shared', 'e' * 64) + self.rows('g-only', 'f' * 64)
            mitos = (self.rows('shared-m', species='Shared', mitos2_input_sequence_sha256='e' * 64) +
                     self.rows('m-only', species='Only M', mitos2_input_sequence_sha256='0' * 64))
            genes, summary, diagnostics = self.run_compare(d, genbank, mitos)
            self.assertEqual((len(genes), len(summary)), (13, 1))
            self.assertEqual({r['reason_comparison_skipped'] for r in diagnostics}, {
                'no_matching_mitos2_sequence_hash', 'no_matching_genbank_sequence_hash'})


if __name__ == '__main__':
    unittest.main()
