#!/usr/bin/env python3
"""Reference-discovery entry point for primate mtDNA preprocessing.

This wrapper preserves the expected CLI/output contract for the discovery step. If
an older full discovery implementation is available as find_primate_wg_chrM_refs.v5.py
in the repository root, this script delegates to it. Otherwise it creates a
schema-valid manual-review manifest from the supplied species table so downstream
preprocessing can be developed and tested without network discovery.
"""
import argparse, csv, os, subprocess, sys
from pathlib import Path

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--species', required=True); p.add_argument('--mito-fasta', required=True)
    p.add_argument('--tree', required=True); p.add_argument('--outdir', required=True)
    p.add_argument('--email', required=True); p.add_argument('--max-nearest', default='200')
    p.add_argument('--delay', default='0.34')
    a=p.parse_args()
    legacy=Path('find_primate_wg_chrM_refs.v5.py')
    if legacy.exists():
        os.execvp('python3', ['python3', str(legacy), '--species', a.species, '--mito-fasta', a.mito_fasta, '--tree', a.tree, '--outdir', a.outdir, '--email', a.email, '--max-nearest', str(a.max_nearest), '--delay', str(a.delay)])
    outdir=Path(a.outdir); (outdir/'cache').mkdir(parents=True, exist_ok=True)
    with open(a.species, newline='') as fh: rows=list(csv.DictReader(fh, delimiter='\t'))
    cols=['target_species','sample_count','preprint_REFERENCE_SPECIES','final_wg_ref_species','final_wg_ref_source','final_wg_assembly_accession','final_wg_assembly_level','final_wg_ftp_path','final_chrM_species','final_chrM_source','final_chrM_accession','final_chrM_contig_name','final_chrM_length','final_chrM_assembly_accession','final_chrM_assembly_source','final_chrM_genbank_accn','final_chrM_refseq_accn','final_reference_strategy']
    with open(outdir/'species_reference_chrM_summary.tsv','w',newline='') as fh:
        w=csv.DictWriter(fh, fieldnames=cols, delimiter='\t'); w.writeheader()
        for r in rows:
            w.writerow({'target_species':r.get('species') or r.get('target_species',''), 'sample_count':r.get('sample_count',''), 'preprint_REFERENCE_SPECIES':r.get('preprint_REFERENCE_SPECIES') or r.get('REFERENCE_SPECIES',''), 'final_reference_strategy':'manual_discovery_required'})
    for name in ['all_candidate_wg_refs.tsv','nuccore_mito_hits.tsv']:
        (outdir/name).write_text('status\tmessage\nmanual_discovery_required\tNo legacy discovery implementation found.\n')
    (outdir/'species_reference_chrM_summary.status_counts.tsv').write_text(f'status\tn\nmanual_discovery_required\t{len(rows)}\n')
if __name__=='__main__': main()
