import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

ROOT = Path(__file__).resolve().parents[2]


@unittest.skipUnless(BIOPYTHON_AVAILABLE, 'Biopython is required for GenBank fixture tests')
class BuildPrimateCodonTableTests(unittest.TestCase):
    def make_record(self, path):
        # Plus CDS at 1..6 is ATGAAA. Minus CDS positions 7..12 encode ATGCCC.
        record = SeqRecord(Seq('ATGAAAGGGCAT'), id='TEST.1', name='TEST')
        record.annotations['molecule_type'] = 'DNA'
        record.features = [
            SeqFeature(FeatureLocation(0, 6, strand=1), type='CDS', qualifiers={'gene':['ND1'], 'codon_start':['1'], 'transl_table':['2']}),
            SeqFeature(CompoundLocation([FeatureLocation(6, 9, strand=-1), FeatureLocation(9, 12, strand=-1)]), type='CDS', qualifiers={'gene':['COI'], 'codon_start':['1']}),
        ]
        SeqIO.write(record, path, 'genbank')

    def test_local_genbank_builds_coding_orientation_and_duplicates_samples(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td); gbdir = d / 'gb'; gbdir.mkdir(); self.make_record(gbdir / 'TEST.1.gb')
            refs = d / 'refs.tsv'; refs.write_text('sample\tspecies\taccession\tfamily\nS1\tSpecies one\tTEST.1\tFam\nS2\tSpecies one\tTEST.1\tFam\n')
            config = d / 'config.yaml'; output = d / 'table.tsv'; failures = d / 'failed.tsv'; summary = d / 'summary.tsv'
            config.write_text(f'''build_primate_codon_table:
  paths:
    sample_ref_file: {refs}
    genbank_dir: {gbdir}
    output_table: {output}
    failed_downloads_table: {failures}
    summary_table: {summary}
  settings:
    accession_columns: accession,reference_id
    skip_existing_genbank: true
''')
            result = subprocess.run([sys.executable, str(ROOT / 'qc_analysis/scripts/build_primate_codon_table.py'), '--config', str(config)], cwd=ROOT, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            rows = list(csv.DictReader(output.open(), delimiter='\t'))
            self.assertEqual({r['sample'] for r in rows}, {'S1', 'S2'})
            s1 = [r for r in rows if r['sample'] == 'S1']
            self.assertEqual(len(s1), 12)
            plus = next(r for r in s1 if r['pos'] == '1')
            self.assertEqual((plus['codon_seq'], plus['codon_pos_in_triplet'], plus['codon_pos1_genomic']), ('ATG', '1', '1'))
            minus = next(r for r in s1 if r['pos'] == '12')
            self.assertEqual((minus['codon_seq'], minus['codon_pos_in_triplet']), ('ATG', '1'))
            self.assertEqual((minus['codon_pos1_genomic'], minus['codon_pos2_genomic'], minus['codon_pos3_genomic']), ('12', '11', '10'))
            self.assertEqual(minus['gene'], 'MT-CO1')

    def test_missing_accession_is_reported_without_stopping(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td); refs = d / 'refs.tsv'; refs.write_text('sample\tspecies\nmissing\tSpecies\n')
            config = d / 'config.yaml'; failures = d / 'failed.tsv'
            config.write_text(f'''build_primate_codon_table:
  paths:
    sample_ref_file: {refs}
    genbank_dir: {d / 'gb'}
    output_table: {d / 'table.tsv'}
    failed_downloads_table: {failures}
    summary_table: {d / 'summary.tsv'}
  settings: {{accession_columns: accession, sample_column: sample, species_column: species}}
''')
            result = subprocess.run([sys.executable, str(ROOT / 'qc_analysis/scripts/build_primate_codon_table.py'), '--config', str(config)], cwd=ROOT, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            failure = next(csv.DictReader(failures.open(), delimiter='\t'))
            self.assertEqual(failure['sample'], 'missing')
            self.assertIn('No accession', failure['reason'])


if __name__ == '__main__':
    unittest.main()
