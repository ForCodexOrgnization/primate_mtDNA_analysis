import csv, tempfile, unittest
from pathlib import Path
from qc_analysis.lib.mt_anchor_utils import rotate_sequence, sequence_sha256
from qc_analysis.scripts.discover_global_liftover_anchor import main as discover_main


def write_fa(p, name, seq): p.write_text(f'>{name}\n{seq}\n')
def rows(p):
    with p.open() as h: return list(csv.DictReader(h, delimiter='\t'))

class GlobalAnchorDiscoveryTests(unittest.TestCase):
    def run_discovery(self, manifest_rows, refs, min_len=20):
        td=tempfile.TemporaryDirectory(); d=Path(td.name)
        human_seq='AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT'
        write_fa(d/'human.fa','chrM',human_seq)
        fasta_paths={}
        for name,seq in refs.items():
            fp=d/f'{name}.fa'; write_fa(fp,'chrM',seq); fasta_paths[name]=fp
        man=d/'samples.tsv'
        with man.open('w') as h:
            h.write('sample\tspecies\treference_id\tspecies_fasta\n')
            for sample,species,rid,refname in manifest_rows:
                h.write(f'{sample}\t{species}\t{rid}\t{fasta_paths[refname]}\n')
        cfg=d/'cfg.yaml'; out=d/'ga'
        cfg.write_text(f'''global_anchor_discovery:\n  output_dir: {out}\n  sample_ref_file: {man}\n  human_fasta: {d/'human.fa'}\n  min_reference_length: {min_len}\n  max_reference_length: 100\n  max_n_fraction: 0.2\n  candidate_window_size: 5\n  min_anchor_column_occupancy: 0.66\n  min_window_mean_occupancy: 0.66\n  max_window_gap_fraction: 0.34\n  max_homopolymer_run: 20\n  aligner: mafft\n  aligner_options: "--auto --quiet"\n  allow_simple_alignment_fallback: true\n''')
        discover_main(['--config', str(cfg)])
        return td,out

    def test_identical_circular_sequence_different_starts(self):
        base='AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT'
        refs={'r1':rotate_sequence(base,1),'r2':rotate_sequence(base,9)}
        td,out=self.run_discovery([('A','sp','r1','r1'),('B','sp','r2','r2')], refs)
        pos=rows(out/'reference_anchor_positions.tsv')
        self.assertEqual(len(pos),2)
        self.assertNotEqual(pos[0]['anchor_original_position'], pos[1]['anchor_original_position'])
        starts={rotate_sequence(refs[r['reference_id']], int(r['anchor_original_position']))[:8] for r in pos}
        self.assertEqual(len(starts),1)
        td.cleanup()

    def test_duplicate_reference_sequence_deduplicates(self):
        refs={'r1':'AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT'}
        td,out=self.run_discovery([('A','sp','r1','r1'),('B','sp','r1','r1')], refs)
        man=rows(out/'unique_reference_manifest.tsv')
        self.assertEqual(len(man),1); self.assertEqual(man[0]['sample_count'],'2')
        td.cleanup()

    def test_low_quality_truncated_excluded(self):
        refs={'short':'ACGTACGT'}
        td,out=self.run_discovery([('A','sp','short','short')], refs, min_len=20)
        exc=rows(out/'excluded_references.tsv')
        self.assertIn('REFERENCE_TOO_SHORT', exc[0]['exclusion_reason'])
        td.cleanup()

    def test_deterministic_input_ordering(self):
        refs={'a':'AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT','b':'AAAACCCCGGGGTTTTAAAACCCCGGGGTTTA'}
        td1,out1=self.run_discovery([('A','sp','a','a'),('B','sp','b','b')], refs)
        td2,out2=self.run_discovery([('B','sp','b','b'),('A','sp','a','a')], refs)
        self.assertEqual((out1/'global_anchor_selection.tsv').read_text(), (out2/'global_anchor_selection.tsv').read_text())
        td1.cleanup(); td2.cleanup()

if __name__=='__main__': unittest.main()
